const std = @import("std");
const png = @import("png.zig");
const kd_tree = @import("kd_tree.zig");
const vp_tree = @import("vp_tree.zig");
const colours_lib = @import("colours.zig");
const flags = @import("flags.zig");
const Pixel = @import("colours.zig").Pixel;

const seed_flag = flags.ValueFlag(
    u32,
    "seed",
    "e",
    0,
    "The seed for the RNG",
);

const starts_flag = flags.ValueFlag(
    u16,
    "starts",
    "n",
    1,
    "The number of pixels to place before attempting nearest-neighbour",
);
const wrap_flag = flags.BooleanFlag(
    "wrap",
    false,
    "If nearest-neighbour checks should wrap around the image boundary",
);
const output_file_flag = flags.ValueFlag(
    []const u8,
    "output",
    "o",
    "out.png",
    "The output PNG to write to",
);
const sort_flag = flags.EnumFlag(
    ColourSort,
    "sort_type",
    .hue,
    "How to sort the initial set of colours before attempting nearest-neighbour",
);
const verify_flag = flags.BooleanFlag(
    "verify",
    true,
    "Run verification after generating the output",
);
const depth_flag = flags.ValueFlag(
    u8,
    "depth",
    "d",
    8,
    "The bit depth of each RGB channel - valid range is [1, 8] (inclusive)",
);
const approximate_flag = flags.BooleanFlag(
    "approx",
    false,
    "Approximate nearest-neighbour",
);
const tree_type_flag = flags.EnumFlag(
    TreeType,
    "tree_type",
    .kd,
    "The type of tree to use",
);
const stride_flag = flags.ValueFlag(
    u8,
    "stride",
    "r",
    1,
    "The step size between successive pixels when iterating through for nearest-neighbbour",
);
const source_flag = flags.ValueFlag(
    ?[]const u8,
    "source",
    "s",
    null,
    "The source image to use for colours",
);

const args_parser = flags.ArgParser(.{
    output_file_flag,
    source_flag,
    depth_flag,
    seed_flag,
    sort_flag,
    tree_type_flag,
    starts_flag,
    wrap_flag,
    approximate_flag,
    stride_flag,
    verify_flag,
    flags.helpFlag,
});

const ImageData = struct {
    buffer: []Pixel,
    size_x: usize,
    size_y: usize,

    fn at(self: *const ImageData, coord: ImageCoord) Pixel {
        const idx = self.toIndex(coord);
        return self.buffer[idx];
    }

    fn put(self: *ImageData, coord: ImageCoord, value: Pixel) void {
        const idx = self.toIndex(coord);
        self.buffer[idx] = value;
    }

    fn toIndex(self: *const ImageData, coord: ImageCoord) usize {
        const bounded = self.boundCoord(coord);
        return self.size_x * bounded.y + bounded.x;
    }

    // Ensures the coordinate is within the bounds of the image data.
    // Wraps with the image bounds.
    fn boundCoord(self: *const ImageData, coord: ImageCoord) ImageCoord {
        return .{
            .x = coord.x % self.size_x,
            .y = coord.y % self.size_y,
        };
    }
};

const ImageCoord = struct {
    x: usize,
    y: usize,
};

fn StridedIterator(comptime T: type) type {
    return struct {
        data: []const T,
        start: usize = 0,
        count: usize = 0,
        strides: usize = 1,
        current_stride: usize = 0,

        pub fn next(it: *StridedIterator(T)) ?T {
            if (it.count >= it.data.len) {
                it.current_stride += 1;
                if (it.current_stride >= it.strides) {
                    return null;
                }
                it.count = it.current_stride;
            }
            const current = (it.start + it.count) % it.data.len;
            it.count += it.strides;
            return it.data[current];
        }

        pub fn iterate(data: []const T, start: usize, stride: usize) StridedIterator(T) {
            const clamped_stride = if (stride > 0) stride else 1;
            return .{
                .data = data,
                .start = start,
                .count = 0,
                .strides = clamped_stride,
                .current_stride = 0,
            };
        }
    };
}

const ChannelIterator = union(enum) {
    repeat: RepeatingIterator,
    reversing: ReversingIterator,

    pub fn next(it: *ChannelIterator) usize {
        return switch (it.*) {
            inline else => |*impl| impl.next(),
        };
    }
};

const RepeatingIterator = struct {
    count: usize = 256,
    idx: usize = 0,
    step: usize = 1,

    pub fn next(it: *RepeatingIterator) usize {
        const v = it.idx * it.step;
        it.idx += 1;
        it.idx %= it.count;
        return v;
    }

    pub fn iter(self: *RepeatingIterator) ChannelIterator {
        return .{
            .self = self,
            .next = RepeatingIterator.next,
        };
    }
};

const ReversingIterator = struct {
    forwards: bool = true,
    count: usize = 256,
    idx: usize = 0,
    step: usize = 1,

    pub fn next(it: *ReversingIterator) usize {
        const max = it.count - 1;
        const v = it.idx * it.step;
        if (it.forwards and it.idx == max) {
            it.forwards = false;
        } else if (!it.forwards and it.idx == 0) {
            it.forwards = true;
        } else if (it.forwards) {
            it.idx += 1;
        } else {
            it.idx -= 1;
        }
        return v;
    }

    pub fn iter(self: *ReversingIterator) ChannelIterator {
        return .{
            .self = self,
            .next = &ReversingIterator.next,
        };
    }
};

const TreeType = enum {
    kd,
    vp,
    vp2,
    vpinf,

    fn treeType(comptime tree_type: TreeType) type {
        return treeTypeValue(tree_type, ImageCoord);
    }

    fn treeTypeValue(comptime tree_type: TreeType, V: type) type {
        return switch (tree_type) {
            .kd => kd_tree.KdTree(Pixel, V),
            .vp => vp_tree.VpTree(Pixel, V, Pixel.l1Distance),
            .vp2 => vp_tree.VpTree(Pixel, V, Pixel.l2Distance),
            .vpinf => vp_tree.VpTree(Pixel, V, Pixel.lInfDistance),
        };
    }
};

fn TreeStore(comptime tree_type: TreeType) type {
    return struct {
        tree: TreeType.treeType(tree_type),

        pub fn init() @This() {
            return .{ .tree = .{} };
        }

        pub fn getNearest(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            return self.tree.getNearest(key);
        }

        pub fn getNear(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            return self.tree.getNear(key);
        }

        pub fn get(self: *@This(), key: Pixel) ?ImageCoord {
            return self.tree.get(key);
        }

        pub fn getMutable(self: *@This(), key: Pixel) ?*ImageCoord {
            return self.tree.get(key);
        }

        pub fn add(
            self: *@This(),
            allocator: std.mem.Allocator,
            key: Pixel,
            value: ImageCoord,
        ) !void {
            return self.tree.add(allocator, key, value);
        }

        pub fn remove(self: *@This(), key: Pixel, _: ImageCoord) !void {
            return self.tree.remove(key);
        }

        pub fn leaves(self: *@This()) usize {
            return self.tree.leaf_count;
        }

        pub fn emptyLeaves(self: *@This()) usize {
            return self.tree.empty_leaf_count;
        }
    };
}

fn TreeMultiStore(comptime tree_type: TreeType) type {
    return struct {
        tree: TreeType.treeTypeValue(tree_type, std.ArrayList(ImageCoord)),
        keys_count: usize = 0,
        values_count: usize = 0,

        fn init() @This() {
            return .{ .tree = .{}, .keys_count = 0, .values_count = 0 };
        }

        fn getNearest(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const data = self.tree.getNearest(key) orelse return null;
            if (data.items.len == 0) {
                std.debug.print("SHOULD NEVER HAPPEN\n", .{});
                return null;
            }
            return .{ found, data.items[0] };
        }

        fn getNear(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const data = self.tree.getNear(key) orelse return null;
            if (data.items.len == 0) {
                return null;
            }
            return .{ found, data.items[0] };
        }

        fn add(
            self: *@This(),
            allocator: std.mem.Allocator,
            key: Pixel,
            value: ImageCoord,
        ) !void {
            const maybe_data = self.tree.getMutable(key);
            if (maybe_data == null) {
                var data: std.ArrayList(ImageCoord) = .empty;
                try data.append(allocator, value);
                try self.tree.add(allocator, key, data);
                self.keys_count += 1;
            } else {
                try maybe_data.?.append(allocator, value);
            }
            self.values_count += 1;
        }

        fn remove(self: *@This(), key: Pixel, value: ImageCoord) !void {
            const data = self.tree.getMutable(key) orelse return error.MissingKey;
            for (0..data.items.len) |idx| {
                if (std.meta.eql(value, data.items[idx])) {
                    _ = data.swapRemove(idx);
                    break;
                }
            }
            self.values_count -= 1;
            if (data.items.len == 0) {
                self.tree.remove(key) catch unreachable;
                self.keys_count -= 1;
            }
        }

        fn leaves(self: *@This()) usize {
            return self.tree.leaf_count;
        }

        fn emptyLeaves(self: *@This()) usize {
            return self.tree.empty_leaf_count;
        }
    };
}

fn createColours(allocator: std.mem.Allocator, channel_depth: u8, zigzag: bool) ![]Pixel {
    const channel_size: usize = @as(usize, 1) << @as(u4, @intCast(channel_depth));
    const data_size = channel_size * channel_size * channel_size;
    var colours: []Pixel = try allocator.alloc(Pixel, data_size);
    var idx: usize = 0;
    const step = 256 / channel_size;

    var red_it: ChannelIterator = undefined;
    var green_it: ChannelIterator = undefined;
    var blue_it: ChannelIterator = undefined;
    if (zigzag) {
        red_it = .{ .reversing = .{ .count = channel_size, .step = step } };
        green_it = .{ .reversing = .{ .count = channel_size, .step = step } };
        blue_it = .{ .reversing = .{ .count = channel_size, .step = step } };
    } else {
        red_it = .{ .repeat = .{ .count = channel_size, .step = step } };
        green_it = .{ .repeat = .{ .count = channel_size, .step = step } };
        blue_it = .{ .repeat = .{ .count = channel_size, .step = step } };
    }

    for (0..channel_size) |_| {
        const red = red_it.next();
        for (0..channel_size) |_| {
            const green = green_it.next();
            for (0..channel_size) |_| {
                const blue = blue_it.next();
                const pixel: Pixel = .{
                    .red = @truncate(red),
                    .green = @truncate(green),
                    .blue = @truncate(blue),
                    .alpha = 0xFF,
                };
                colours[idx] = pixel;
                idx += 1;
            }
        }
    }
    return colours;
}

// Returns a subslice of the provided buffer populate with open neighbours the
// given point. This accounts for wrapping at the image boundaries.
fn getAvailableNeighboursWrapped(
    image: ImageData,
    point: ImageCoord,
    buffer: []ImageCoord,
) []ImageCoord {
    var candidates = std.mem.zeroes([8]ImageCoord);
    var c_idx: usize = 0;
    for (0..3) |y_idx| {
        const dy = @as(i32, @intCast(y_idx)) - 1;
        const y: i32 = @mod(
            @as(i32, @intCast(point.y)) - dy,
            @as(i32, @intCast(image.size_y)),
        );
        for (0..3) |x_idx| {
            const dx = @as(i32, @intCast(x_idx)) - 1;
            const x: i32 = @mod(
                @as(i32, @intCast(point.x)) - dx,
                @as(i32, @intCast(image.size_x)),
            );

            if (x == point.x and y == point.y) {
                continue;
            }
            candidates[c_idx] = .{
                .x = @as(usize, @intCast(x)),
                .y = @as(usize, @intCast(y)),
            };
            c_idx += 1;
        }
    }

    var out_idx: usize = 0;
    for (0..c_idx) |idx| {
        const test_point = candidates[idx];
        if (image.at(test_point).alpha == 0) {
            buffer[out_idx] = test_point;
            out_idx += 1;
        }
    }
    return buffer[0..out_idx];
}

/// Returns a subslice of the provided buffer populated with the open neighbours
/// of the given point.
fn getAvailableNeighboursClamped(
    image: ImageData,
    point: ImageCoord,
    buffer: []ImageCoord,
) []ImageCoord {
    const bounded = image.boundCoord(point);
    var idx: usize = 0;

    const left_edge = bounded.x == 0;
    const right_edge = bounded.x == image.size_x - 1;
    const top_edge = bounded.y == 0;
    const bottom_edge = bounded.y == image.size_y - 1;

    var candidates = std.mem.zeroes([8]ImageCoord);
    if (!top_edge and !left_edge) {
        candidates[idx] = .{ .x = bounded.x - 1, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!top_edge) {
        candidates[idx] = .{ .x = bounded.x, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!top_edge and !right_edge) {
        candidates[idx] = .{ .x = bounded.x + 1, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!left_edge) {
        candidates[idx] = .{ .x = bounded.x - 1, .y = bounded.y };
        idx += 1;
    }
    if (!right_edge) {
        candidates[idx] = .{ .x = bounded.x + 1, .y = bounded.y };
        idx += 1;
    }
    if (!left_edge and !bottom_edge) {
        candidates[idx] = .{ .x = bounded.x - 1, .y = bounded.y + 1 };
        idx += 1;
    }
    if (!bottom_edge) {
        candidates[idx] = .{ .x = bounded.x, .y = bounded.y + 1 };
        idx += 1;
    }
    if (!right_edge and !bottom_edge) {
        candidates[idx] = .{ .x = bounded.x + 1, .y = bounded.y + 1 };
        idx += 1;
    }
    var out_idx: usize = 0;
    for (0..idx) |candidate_idx| {
        const test_point = candidates[candidate_idx];
        if (image.at(test_point).alpha == 0) {
            buffer[out_idx] = test_point;
            out_idx += 1;
        }
    }
    return buffer[0..out_idx];
}

fn getAvailableNeighbours(image: ImageData, point: ImageCoord, buffer: []ImageCoord, wrap: bool) []ImageCoord {
    if (wrap) {
        return getAvailableNeighboursWrapped(image, point, buffer);
    } else {
        return getAvailableNeighboursClamped(image, point, buffer);
    }
}

fn getFirstNeighbour(image: ImageData, point: ImageCoord, wrap: bool) ?ImageCoord {
    var buffer: [8]ImageCoord = std.mem.zeroes([8]ImageCoord);
    const available = getAvailableNeighbours(image, point, &buffer, wrap);
    if (available.len > 0) {
        return available[0];
    }
    return null;
}

const FillParameters = struct {
    seed: u32,
    starts: u16,
    stride: u8,
    approximate: bool,
    wrap: bool,
};

const ProgressMonitor = struct {
    writer: *std.Io.Writer,
    used: bool = false,

    pub fn update(self: *@This(), comptime format: []const u8, args: anytype) !void {
        const csi = "\x1B[";
        const clear_line = csi ++ "1K";
        const text = clear_line ++ "\r" ++ format;
        try self.writer.print(text, args);
        self.used = true;
    }

    pub fn finish(self: *@This()) void {
        self.writer.print("\n", .{}) catch return;
    }
};
/// Creates a new tree using the partial image.
fn populateNewTree(
    comptime tree_type: type,
    allocator: std.mem.Allocator,
    image: ImageData,
    parameters: FillParameters,
) !tree_type {
    var tree: tree_type = tree_type.init();
    for (0..image.size_y) |y| {
        for (0..image.size_x) |x| {
            const slot: ImageCoord = .{ .x = x, .y = y };
            const pixel = image.at(slot);
            if (pixel.alpha == 0) {
                continue;
            }
            if (getFirstNeighbour(image, slot, parameters.wrap) == null) {
                continue;
            }
            try tree.add(allocator, pixel, slot);
        }
    }
    return tree;
}

fn fillImage(
    comptime tree_type: type,
    allocator: std.mem.Allocator,
    colours: []const Pixel,
    image: *ImageData,
    progress_monitor: *ProgressMonitor,
    parameters: FillParameters,
) !void {
    var prng = std.Random.DefaultPrng.init(parameters.seed);
    const rng = prng.random();
    for (0..image.buffer.len) |idx| {
        image.buffer[idx] = .{};
    }
    defer progress_monitor.finish();
    var tree_alloc = std.heap.ArenaAllocator.init(allocator);
    defer tree_alloc.deinit();
    const tree_allocator = tree_alloc.allocator();
    var tree: tree_type = try populateNewTree(tree_type, tree_allocator, image.*, parameters);

    const start_idx = rng.intRangeLessThan(usize, 0, colours.len);
    var colours_it = StridedIterator(Pixel).iterate(
        colours,
        start_idx,
        parameters.stride,
    );

    // place the initial pixels and seed the tree
    for (0..parameters.starts) |_| {
        const initial_x = rng.intRangeLessThan(usize, 0, image.size_x);
        const initial_y = rng.intRangeLessThan(usize, 0, image.size_y);
        const initial_pixel = colours_it.next() orelse return error.NoColours;
        const initial_coord = ImageCoord{ .x = initial_x, .y = initial_y };
        try tree.add(tree_allocator, initial_pixel, initial_coord);
        image.put(initial_coord, initial_pixel);
    }
    var c_count: usize = parameters.starts;
    var since_rebuild: usize = 0;
    var rebuilds: usize = 0;
    while (colours_it.next()) |colour| {
        since_rebuild += 1;
        const empties = tree.emptyLeaves();
        const ratio: f32 = @as(f32, @floatFromInt(empties)) / @as(f32, @floatFromInt(tree.leaves()));
        const rebuild = (empties > 1 and ratio >= 0.1) or (since_rebuild > 1024 * 1024);
        if (rebuild) {
            _ = tree_alloc.reset(.retain_capacity);
            tree = try populateNewTree(tree_type, tree_allocator, image.*, parameters);
            since_rebuild = 0;
            rebuilds += 1;
        }

        var slot: ImageCoord = .{ .x = 0, .y = 0 };
        var pixel: Pixel = .{};
        const search_function = if (parameters.approximate)
            &tree_type.getNear
        else
            &tree_type.getNearest;
        while (true) {
            const closest_pixel, const closest_idx = search_function(
                &tree,
                colour,
            ) orelse return error.InvalidStateEmptyTree;
            var buffer = std.mem.zeroes([8]ImageCoord);
            const available = getAvailableNeighbours(image.*, closest_idx, &buffer, parameters.wrap);
            if (available.len == 0) {
                try tree.remove(closest_pixel, closest_idx);
                continue;
            }
            // we're about to remove this pixel's last available neighbour
            if (available.len == 1) {
                try tree.remove(closest_pixel, closest_idx);
            }
            const pick_idx = rng.intRangeLessThan(usize, 0, available.len);
            slot = available[pick_idx];
            pixel = closest_pixel;
            break;
        }
        image.put(slot, colour);
        if (getFirstNeighbour(image.*, slot, parameters.wrap) != null) {
            try tree.add(tree_allocator, colour, slot);
        }
        c_count += 1;
        if (c_count % 1024 * 16 == 0) {
            const percentage: f64 = 100 * @as(f64, @floatFromInt(c_count)) / @as(f64, @floatFromInt(colours.len));
            try progress_monitor.update(
                "{d:5.1}% - {d:6}/{d} \tTree leaves: {d:4} - Empty: {d:3} - Rebuilds {d:3}",
                .{
                    percentage,
                    c_count,
                    image.buffer.len,
                    tree.leaves(),
                    tree.emptyLeaves(),
                    rebuilds,
                },
            );
        }
    }
}

fn verifyFullImagePopulated(image: ImageData) bool {
    for (0..image.buffer.len) |idx| {
        const pixel = image.buffer[idx];
        if (pixel.alpha == 0) {
            const x = idx % image.size_x;
            const y = idx / image.size_x;
            std.debug.print("Pixel {any} at {d},{d} was not set\n", .{ pixel, x, y });
            return false;
        }
    }
    return true;
}

fn verifyAllPixels(
    allocator: std.mem.Allocator,
    image: ImageData,
    colours: []colours_lib.Pixel,
) !bool {
    var found: std.AutoHashMapUnmanaged(Pixel, usize) = .{};
    defer found.deinit(allocator);
    for (0..image.buffer.len) |idx| {
        const pixel = image.buffer[idx];
        const count = found.get(pixel) orelse 0;
        try found.put(allocator, pixel, count + 1);
    }

    for (0..colours.len) |idx| {
        const pixel = colours[idx];
        const count = found.get(pixel) orelse return false;
        if (count == 0) {
            return false;
        } else if (count == 0) {
            _ = found.remove(pixel);
        } else {
            try found.put(allocator, pixel, count - 1);
        }
    }
    var it = found.valueIterator();
    while (it.next()) |count| {
        if (count.* != 0) {
            return false;
        }
    }
    return true;
}

const ColourSort = enum {
    hue,
    hsp,
    zigzag,
    zorder,
    none,
};

const Parameters = struct {
    seed: u32 = 0,
    starts: u16 = 0,
    sort_type: ColourSort = ColourSort.hue,
    output_file: []const u8 = "",
    verify: bool,
    channel_depth: u8 = 0,
    stride: u8 = 1,
    tree_type: TreeType = .kd,
    approximate: bool = false,
    wrap: bool = false,
    source: ?[]const u8 = null,
    help: bool = false,
};

fn parseArguments(args: std.process.Args) !Parameters {
    const params = args_parser.parse(args) catch return error.InvaldArguments;

    if (params.depth == 0 or params.depth > 8) {
        return error.InvalidArguments;
    }

    return .{
        .seed = params.seed,
        .starts = params.starts,
        .output_file = params.output,
        .sort_type = params.sort_type,
        .verify = params.verify,
        .channel_depth = params.depth,
        .stride = params.stride,
        .approximate = params.approx,
        .wrap = params.wrap,
        .tree_type = params.tree_type,
        .source = params.source,
        .help = params.help,
    };
}

fn imageSizeFromBitDepth(channel_depth: u8) struct { u32, u32 } {
    var size_x: u32 = 4096;
    var size_y: u32 = 4096;
    if (channel_depth >= 8) {
        return .{ size_x, size_y };
    }
    const steps = 8 - channel_depth;
    // every drop of 1 bit means our image shrinks by a factor of 8
    for (0..steps) |_| {
        // shrink the larger dimension twice
        if (size_x <= size_y) {
            size_x >>= 1;
            size_y >>= 2;
        } else {
            size_y >>= 1;
            size_x >>= 2;
        }
    }
    return .{ size_x, size_y };
}

pub fn main(init: std.process.Init) !void {
    var gpa = std.heap.DebugAllocator(.{}){};
    const gen_alloc = gpa.allocator();
    var arena = std.heap.ArenaAllocator.init(gen_alloc);
    defer arena.deinit();
    const allocator = arena.allocator();

    const stdout = std.Io.File.stdout();
    const stderr = std.Io.File.stderr();
    // unbuffered - don't plan many writes but the inputs may be large.
    var stdout_writer = stdout.writer(init.io, &.{});
    var stderr_writer = stderr.writer(init.io, &.{});

    var size_x: u32 = 0;
    var size_y: u32 = 0;
    const parameters = parseArguments(init.minimal.args) catch |err| {
        try stdout_writer.interface.print("{s}\n", .{args_parser.helpText()});
        return err;
    };
    if (parameters.help) {
        try stdout_writer.interface.print("{s}\n", .{args_parser.helpText()});
        return;
    }

    var colours: []colours_lib.Pixel = undefined;
    if (parameters.source != null) {
        try stderr_writer.interface.print(
            "Reading from {s}\n",
            .{parameters.source.?},
        );
        const source_data = try png.readPng(
            allocator,
            init.io,
            parameters.source.?,
        );
        colours = source_data.data;
        size_x = source_data.size_x;
        size_y = source_data.size_y;
    } else {
        colours = try createColours(
            allocator,
            parameters.channel_depth,
            parameters.sort_type == ColourSort.zigzag,
        );
        size_x, size_y = imageSizeFromBitDepth(parameters.channel_depth);
    }

    switch (parameters.sort_type) {
        .hue => std.mem.sort(Pixel, colours, {}, colours_lib.hueCompare),
        .hsp => std.mem.sort(Pixel, colours, {}, colours_lib.hspCompare),
        .zigzag => {}, // this is really for generation order
        .zorder => std.mem.sort(Pixel, colours, {}, colours_lib.zOrderCompare),
        .none => {},
    }

    const buffer = try allocator.alloc(Pixel, colours.len);
    var image = ImageData{
        .buffer = buffer,
        .size_x = size_x,
        .size_y = size_y,
    };
    try stderr_writer.interface.print(
        "Producing a {d}x{d} image - {s}\n",
        .{ size_x, size_y, parameters.output_file },
    );

    const fill_params = FillParameters{
        .seed = parameters.seed,
        .approximate = parameters.approximate,
        .stride = parameters.stride,
        .starts = parameters.starts,
        .wrap = parameters.wrap,
    };
    var progress_monitor = ProgressMonitor{ .writer = &stdout_writer.interface };
    switch (parameters.tree_type) {
        inline else => |t| {
            if (parameters.source == null) {
                const tree = TreeStore(t);
                fillImage(
                    tree,
                    gen_alloc,
                    colours,
                    &image,
                    &progress_monitor,
                    fill_params,
                ) catch |err| {
                    stderr_writer.interface.print("Error filling image: {any}\n", .{err}) catch {};
                };
            } else {
                const tree = TreeMultiStore(t);
                fillImage(
                    tree,
                    gen_alloc,
                    colours,
                    &image,
                    &progress_monitor,
                    fill_params,
                ) catch |err| {
                    stderr_writer.interface.print("Error filling image: {any}\n", .{err}) catch {};
                };
            }
        },
    }
    if (parameters.verify) {
        var failed = false;
        if (!verifyFullImagePopulated(image)) {
            try stderr_writer.interface.print("Not all pixels were written to.\n", .{});
            failed = true;
        } else {}
        if (!try verifyAllPixels(allocator, image, colours)) {
            failed = true;
            try stderr_writer.interface.print("Written pixels do not match source pixels\n", .{});
        }
        if (!failed) {
            try stderr_writer.interface.print("Image verified\n", .{});
        }
    }

    try png.writePng(
        allocator,
        init.io,
        buffer,
        parameters.output_file,
        size_x,
        size_y,
    );
}

test "Strided iterator unstrided" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const data = &[_]i32{ 1, 2, 3, 4, 5, 6, 7, 8 };
    var it = StridedIterator(i32).iterate(data, 0, 1);

    var it_out: std.ArrayList(i32) = .empty;
    while (it.next()) |v| {
        try it_out.append(allocator, v);
    }

    try std.testing.expectEqualSlices(i32, data, it_out.items);
}

test "Strided iterator with offset" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const data = &[_]i32{ 1, 2, 3, 4, 5, 6, 7, 8 };
    var it = StridedIterator(i32).iterate(data, 3, 1);

    var it_out: std.ArrayList(i32) = .empty;
    while (it.next()) |v| {
        try it_out.append(allocator, v);
    }

    const expected = &[_]i32{ 4, 5, 6, 7, 8, 1, 2, 3 };
    try std.testing.expectEqualSlices(i32, expected, it_out.items);
}

test "Strided iterator stride and offset" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const data = &[_]i32{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
    var it = StridedIterator(i32).iterate(data, 5, 3);

    var it_out: std.ArrayList(i32) = .empty;
    while (it.next()) |v| {
        try it_out.append(allocator, v);
    }

    const expected = &[_]i32{ 6, 9, 12, 1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3 };
    try std.testing.expectEqual(data.len, expected.len);
    try std.testing.expectEqualSlices(i32, expected, it_out.items);
}
