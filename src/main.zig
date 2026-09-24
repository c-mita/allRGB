const std = @import("std");
const png = @import("png.zig");
const kd_tree = @import("kd_tree.zig");
const vp_tree = @import("vp_tree.zig");
const colours_lib = @import("colours.zig");
const flags = @import("flags.zig");
const fills = @import("fills.zig");

const Pixel = @import("colours.zig").Pixel;
const ImageCoord = @import("image.zig").ImageCoord;
const ImageData = @import("image.zig").ImageData;

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

const target_flag = flags.ValueFlag(
    ?[]const u8,
    "target",
    "t",
    null,
    "A target image to match to",
);

const args_parser = flags.ArgParser(.{
    output_file_flag,
    source_flag,
    target_flag,
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

        pub fn init() @This() {
            return .{ .tree = .{}, .keys_count = 0, .values_count = 0 };
        }

        pub fn getNearest(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const data = self.tree.getNearest(key) orelse return null;
            if (data.items.len == 0) {
                std.debug.print("SHOULD NEVER HAPPEN\n", .{});
                return null;
            }
            return .{ found, data.items[0] };
        }

        pub fn getNear(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const data = self.tree.getNear(key) orelse return null;
            if (data.items.len == 0) {
                return null;
            }
            return .{ found, data.items[0] };
        }

        pub fn add(
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

        pub fn remove(self: *@This(), key: Pixel, value: ImageCoord) !void {
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

        pub fn leaves(self: *@This()) usize {
            return self.tree.leaf_count;
        }

        pub fn emptyLeaves(self: *@This()) usize {
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

/// If necessary, reprocess the colours array to match the target size
fn resizeColours(allocator: std.mem.Allocator, colours: []Pixel, size: usize) ![]Pixel {
    if (colours.len == size) {
        return colours;
    }
    const new_colours: []Pixel = try allocator.alloc(Pixel, size);
    const repeats = size / colours.len;
    const remainder = size - colours.len * repeats;
    for (0..colours.len) |idx| {
        for (0..repeats) |r_idx| {
            new_colours[idx + r_idx] = colours[idx];
        }
    }

    // Evenly distribute the remainder over the input colours
    var new_colours_idx = colours.len * repeats;
    if (remainder > 0) {
        const step = colours.len / remainder;
        for (0..remainder) |r_idx| {
            const idx = step * r_idx;
            new_colours[new_colours_idx] = colours[idx];
            new_colours_idx += 1;
        }
    }

    return new_colours;
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

fn fillImage(
    colours: []const Pixel,
    filler: *fills.ImageFill,
    image: *ImageData,
    progress_monitor: *ProgressMonitor,
    parameters: FillParameters,
) !void {
    for (0..image.buffer.len) |idx| {
        image.buffer[idx] = .{};
    }
    defer progress_monitor.finish();

    const start_idx = filler.rng.intRangeLessThan(usize, 0, colours.len);
    var colours_it = StridedIterator(Pixel).iterate(
        colours,
        start_idx,
        parameters.stride,
    );

    try filler.repopulate();

    // place the initial pixels and seed the tree
    for (0..parameters.starts) |_| {
        const initial_pixel = colours_it.next() orelse return error.NoColours;
        try filler.placeRandomly(initial_pixel);
    }
    var c_count: usize = parameters.starts;
    var since_rebuild: usize = 0;
    var rebuilds: usize = 0;
    while (colours_it.next()) |colour| {
        since_rebuild += 1;
        const empties = filler.emptyLeaves();
        const ratio: f32 = @as(f32, @floatFromInt(empties)) / @as(f32, @floatFromInt(filler.leaves()));
        const rebuild = (empties > 1 and ratio >= 0.1) or (since_rebuild > 1024 * 1024);
        if (rebuild) {
            try filler.repopulate();
            since_rebuild = 0;
            rebuilds += 1;
        }
        try filler.matchAndPlace(colour);
        c_count += 1;
        if (c_count % 1024 * 16 == 0) {
            const percentage: f64 = 100 * @as(f64, @floatFromInt(c_count)) / @as(f64, @floatFromInt(colours.len));
            try progress_monitor.update(
                "{d:5.1}% - {d:6}/{d} \tTree leaves: {d:4} - Empty: {d:3} - Rebuilds {d:3}",
                .{
                    percentage,
                    c_count,
                    image.buffer.len,
                    filler.leaves(),
                    filler.emptyLeaves(),
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
    target: ?[]const u8 = null,
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
        .target = params.target,
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

fn minFillOf(
    allocator: std.mem.Allocator,
    tree_allocator: std.mem.Allocator,
    comptime tree_type: TreeType,
    comptime multi: bool,
    image: *ImageData,
    wrap: bool,
    approximate: bool,
    rng: std.Random,
) !fills.ImageFill {
    const wrapped_tree = if (multi) TreeMultiStore(tree_type) else TreeStore(tree_type);

    var filler = try allocator.create(fills.MinFill(wrapped_tree));
    filler.* = .init(
        tree_allocator,
        image,
        wrap,
        approximate,
        rng,
    );
    return filler.filler();
}

fn targetFillOf(
    allocator: std.mem.Allocator,
    tree_allocator: std.mem.Allocator,
    comptime tree_type: TreeType,
    image: *ImageData,
    target: *ImageData,
    approximate: bool,
    rng: std.Random,
) !fills.ImageFill {
    const wrapped_tree = TreeMultiStore(tree_type);

    var filler = try allocator.create(fills.TargetFill(wrapped_tree));
    filler.* = .init(
        tree_allocator,
        image,
        target,
        approximate,
        rng,
    );
    return filler.filler();
}

/// Create the appropriate filler strategy.
/// The initial allocator is used to create the filler struct
/// The tree_allocator is passed to the created filler to manage
/// its data structures.
fn createFiller(
    allocator: std.mem.Allocator,
    tree_allocator: std.mem.Allocator,
    image: *ImageData,
    tree_type: TreeType,
    fill_type: fills.Fills,
    target: *?ImageData,
    multi: bool,
    wrap: bool,
    approximate: bool,
    rng: std.Random,
) !fills.ImageFill {
    switch (multi) {
        inline else => |is_multi| {
            switch (fill_type) {
                .min => {
                    return switch (tree_type) {
                        inline else => |tree| minFillOf(
                            allocator,
                            tree_allocator,
                            tree,
                            is_multi,
                            image,
                            wrap,
                            approximate,
                            rng,
                        ),
                    };
                },
                .target => {
                    return switch (tree_type) {
                        inline else => |tree| targetFillOf(
                            allocator,
                            tree_allocator,
                            tree,
                            image,
                            &(target.*.?),
                            approximate,
                            rng,
                        ),
                    };
                },
            }
        },
    }
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
    var colours_contains_duplicates: bool = false;
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
        colours_contains_duplicates = true;
    } else {
        colours = try createColours(
            allocator,
            parameters.channel_depth,
            parameters.sort_type == ColourSort.zigzag,
        );
        size_x, size_y = imageSizeFromBitDepth(parameters.channel_depth);
    }

    var target: ?ImageData = null;
    if (parameters.target != null) {
        try stderr_writer.interface.print(
            "Reading from {s}\n",
            .{parameters.target.?},
        );
        const target_data = try png.readPng(
            allocator,
            init.io,
            parameters.target.?,
        );
        // If a target image is set we should match its size
        size_x = target_data.size_x;
        size_y = target_data.size_y;
        target = ImageData{
            .buffer = target_data.data,
            .size_x = target_data.size_x,
            .size_y = target_data.size_y,
        };

        if (target_data.data.len > colours.len) {
            colours_contains_duplicates = true;
        }
        colours = try resizeColours(allocator, colours, target_data.data.len);
    }

    switch (parameters.sort_type) {
        .hue => std.mem.sort(Pixel, colours, {}, colours_lib.hueCompare),
        .hsp => std.mem.sort(Pixel, colours, {}, colours_lib.hspCompare),
        .zigzag => {}, // this is really for generation order
        .zorder => std.mem.sort(Pixel, colours, {}, colours_lib.zOrderCompare),
        .none => {},
    }

    const output_size = size_x * size_y;
    const buffer = try allocator.alloc(Pixel, output_size);
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

    var prng = std.Random.DefaultPrng.init(parameters.seed);
    const rng = prng.random();
    var progress_monitor = ProgressMonitor{ .writer = &stdout_writer.interface };

    var filler_buffer: [128]u8 = undefined;
    var filler_obj_alloc = std.heap.FixedBufferAllocator.init(&filler_buffer);

    var fill_type = fills.Fills.min;
    if (target != null) {
        fill_type = fills.Fills.target;
    }
    var filler = try createFiller(
        filler_obj_alloc.allocator(),
        gen_alloc,
        &image,
        parameters.tree_type,
        fill_type,
        &target,
        colours_contains_duplicates,
        parameters.wrap,
        parameters.approximate,
        rng,
    );

    fillImage(
        colours,
        &filler,
        &image,
        &progress_monitor,
        fill_params,
    ) catch |err| {
        stderr_writer.interface.print("Error filling image: {any}\n", .{err}) catch {};
    };

    try png.writePng(
        allocator,
        init.io,
        buffer,
        parameters.output_file,
        size_x,
        size_y,
    );

    if (parameters.verify) {
        var failed = false;
        if (!verifyFullImagePopulated(image)) {
            try stderr_writer.interface.print("Not all pixels were written to.\n", .{});
            failed = true;
        }
        if (!try verifyAllPixels(allocator, image, colours)) {
            failed = true;
            try stderr_writer.interface.print("Written pixels do not match source pixels\n", .{});
        }
        if (!failed) {
            try stderr_writer.interface.print("Image verified\n", .{});
        }
    }
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
