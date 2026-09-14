const std = @import("std");
const png = @import("png.zig");
const kd_tree = @import("kd_tree.zig");
const vp_tree = @import("vp_tree.zig");
const colours_lib = @import("colours.zig");
const Pixel = @import("colours.zig").Pixel;

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
};

const TreeStore = union(enum) {
    kd: kd_tree.KdTree(Pixel, ImageCoord),
    vp: vp_tree.VpTree(Pixel, ImageCoord, Pixel.l1Distance),

    fn getNearest(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
        return switch (self.*) {
            inline else => |*impl| impl.getNearest(key),
        };
    }

    fn getNear(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
        return switch (self.*) {
            inline else => |*impl| impl.getNear(key),
        };
    }

    fn add(self: *@This(), allocator: std.mem.Allocator, key: Pixel, value: ImageCoord) !void {
        return switch (self.*) {
            inline else => |*impl| impl.add(allocator, key, value),
        };
    }

    fn remove(self: *@This(), key: Pixel) !void {
        return switch (self.*) {
            inline else => |*impl| impl.remove(key),
        };
    }

    fn leaves(self: *@This()) usize {
        return switch (self.*) {
            inline else => |*impl| impl.leaf_count,
        };
    }

    fn emptyLeaves(self: *@This()) usize {
        return switch (self.*) {
            inline else => |*impl| impl.empty_leaf_count,
        };
    }
};

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

/// Returns a subslice of the provided buffer populated with the open neighbours
/// of the given point.
fn getAvailableNeighbours(image: ImageData, point: ImageCoord, buffer: []ImageCoord) []ImageCoord {
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

fn getFirstNeighbour(image: ImageData, point: ImageCoord) ?ImageCoord {
    var buffer: [8]ImageCoord = std.mem.zeroes([8]ImageCoord);
    const available = getAvailableNeighbours(image, point, &buffer);
    if (available.len > 0) {
        return available[0];
    }
    return null;
}

/// Creates a new tree using the partial image.
fn populateNewTree(
    allocator: std.mem.Allocator,
    image: ImageData,
    tree_type: TreeType,
) !TreeStore {
    var tree: TreeStore = switch (tree_type) {
        .kd => .{ .kd = .{} },
        .vp => .{ .vp = .{} },
    };
    for (0..image.size_y) |y| {
        for (0..image.size_x) |x| {
            const slot: ImageCoord = .{ .x = x, .y = y };
            const pixel = image.at(slot);
            if (pixel.alpha == 0) {
                continue;
            }
            if (getFirstNeighbour(image, slot) == null) {
                continue;
            }
            try tree.add(allocator, pixel, slot);
        }
    }
    return tree;
}

fn fillImage(
    allocator: std.mem.Allocator,
    colours: []const Pixel,
    image: *ImageData,
    starts: u16,
    stride: u8,
    tree_type: TreeType,
    approximate: bool,
    seed: u32,
) !void {
    var prng = std.Random.DefaultPrng.init(seed);
    const rng = prng.random();
    for (0..image.buffer.len) |idx| {
        image.buffer[idx] = .{};
    }
    var tree_alloc = std.heap.ArenaAllocator.init(allocator);
    defer tree_alloc.deinit();
    const tree_allocator = tree_alloc.allocator();
    var tree: TreeStore = try populateNewTree(tree_allocator, image.*, tree_type);

    const start_idx = rng.intRangeLessThan(usize, 0, colours.len);
    var colours_it = StridedIterator(Pixel).iterate(
        colours,
        start_idx,
        stride,
    );

    // place the initial pixels and seed the tree
    for (0..starts) |_| {
        const initial_x = rng.intRangeLessThan(usize, 0, image.size_x);
        const initial_y = rng.intRangeLessThan(usize, 0, image.size_y);
        const initial_pixel = colours_it.next() orelse return error.NoColours;
        const initial_coord = ImageCoord{ .x = initial_x, .y = initial_y };
        try tree.add(tree_allocator, initial_pixel, initial_coord);
        image.put(initial_coord, initial_pixel);
    }
    var percentage: i32 = 1;
    var c_count: usize = starts;
    var in_tree: usize = starts;
    var since_rebuild: usize = 0;
    const percent_mod = @max(1, colours.len / 100);
    while (colours_it.next()) |colour| {
        since_rebuild += 1;
        if (c_count % percent_mod == 0) {
            std.debug.print("Progress: {d} - Tree leaves: {d} - Empty: {d} - Elements: {d} - Placed: {d}\n", .{
                percentage,
                tree.leaves(),
                tree.emptyLeaves(),
                in_tree,
                c_count,
            });
            percentage += 1;
        }
        c_count += 1;
        const empties = tree.emptyLeaves();
        const ratio: f32 = @as(f32, @floatFromInt(empties)) / @as(f32, @floatFromInt(tree.leaves()));
        const rebuild = (empties > 1 and ratio >= 0.1) or (since_rebuild > 1024 * 1024);
        if (rebuild) {
            std.debug.print("Rebuilding tree - leaves: {any} - empty {any}\n", .{ tree.leaves(), tree.emptyLeaves() });
            _ = tree_alloc.reset(.retain_capacity);
            tree = try populateNewTree(tree_allocator, image.*, tree_type);
            since_rebuild = 0;
        }

        var slot: ImageCoord = .{ .x = 0, .y = 0 };
        var pixel: Pixel = .{};
        const search_function = if (approximate)
            &TreeStore.getNear
        else
            &TreeStore.getNearest;
        while (true) {
            const closest_pixel, const closest_idx = search_function(
                &tree,
                colour,
            ) orelse return error.InvalidStateEmptyTree;
            var buffer = std.mem.zeroes([8]ImageCoord);
            const available = getAvailableNeighbours(image.*, closest_idx, &buffer);
            if (available.len == 0) {
                try tree.remove(closest_pixel);
                in_tree -= 1;
                continue;
            }
            // we're about to remove this pixel's last available neighbour
            if (available.len == 1) {
                try tree.remove(closest_pixel);
                in_tree -= 1;
            }
            const pick_idx = rng.intRangeLessThan(usize, 0, available.len);
            slot = available[pick_idx];
            pixel = closest_pixel;
            break;
        }
        image.put(slot, colour);
        if (getFirstNeighbour(image.*, slot) != null) {
            try tree.add(tree_allocator, colour, slot);
            in_tree += 1;
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

fn verifyUniquePixels(allocator: std.mem.Allocator, image: ImageData) !bool {
    var found: std.AutoHashMapUnmanaged(Pixel, bool) = .{};
    defer found.deinit(allocator);
    for (0..image.buffer.len) |idx| {
        const pixel = image.buffer[idx];
        if (found.contains(pixel)) {
            return false;
        }
        try found.put(allocator, pixel, true);
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
};

fn parseArguments(args: std.process.Args) !Parameters {
    var seed: u32 = 0;
    var starts: u16 = 1;
    var output_file: []const u8 = "out.png";
    var sort_type = ColourSort.hue;
    var verify = false;
    var channel_depth: u8 = 8;
    var stride: u8 = 1;
    var approximate: bool = false;
    var tree_type: TreeType = .kd;
    var it = args.iterate();

    while (it.next()) |arg| {
        if (std.mem.eql(u8, "--seed", arg)) {
            const seed_str = it.next() orelse return error.InvalidArguments;
            seed = try std.fmt.parseInt(u32, seed_str, 10);
        } else if (std.mem.eql(u8, "--starts", arg)) {
            const starts_str = it.next() orelse return error.InvalidArguments;
            starts = try std.fmt.parseInt(u16, starts_str, 10);
        } else if (std.mem.eql(u8, "--out", arg)) {
            output_file = it.next() orelse return error.InvalidArguments;
        } else if (std.mem.eql(u8, "--hue", arg)) {
            sort_type = ColourSort.hue;
        } else if (std.mem.eql(u8, "--hsp", arg)) {
            sort_type = ColourSort.hsp;
        } else if (std.mem.eql(u8, "--zigzag", arg)) {
            sort_type = ColourSort.zigzag;
        } else if (std.mem.eql(u8, "--zorder", arg)) {
            sort_type = ColourSort.zorder;
        } else if (std.mem.eql(u8, "--none", arg)) {
            sort_type = ColourSort.none;
        } else if (std.mem.eql(u8, "--verify", arg)) {
            verify = true;
        } else if (std.mem.eql(u8, "--depth", arg)) {
            const depth_str = it.next() orelse return error.InvalidArguments;
            channel_depth = std.fmt.parseInt(u8, depth_str, 10) catch return error.InvalidArguments;
            if (channel_depth > 8) {
                std.debug.print("Maximum channel depth is 8\n", .{});
                return error.InvalidArguments;
            } else if (channel_depth == 0) {
                std.debug.print("Minimum channel depth is 1\n", .{});
                return error.InvalidArguments;
            }
        } else if (std.mem.eql(u8, "--stride", arg)) {
            const stride_str = it.next() orelse return error.InvalidArguments;
            stride = std.fmt.parseInt(u8, stride_str, 10) catch return error.InvalidArguments;
        } else if (std.mem.eql(u8, "--kd", arg)) {
            tree_type = .kd;
        } else if (std.mem.eql(u8, "--vp", arg)) {
            tree_type = .vp;
        } else if (std.mem.eql(u8, "--approx", arg)) {
            approximate = true;
        }
    }

    return .{
        .seed = seed,
        .starts = starts,
        .output_file = output_file,
        .sort_type = sort_type,
        .verify = verify,
        .channel_depth = channel_depth,
        .stride = stride,
        .approximate = approximate,
        .tree_type = tree_type,
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

    const parameters = try parseArguments(init.minimal.args);

    const colours = try createColours(allocator, parameters.channel_depth, parameters.sort_type == ColourSort.zigzag);

    switch (parameters.sort_type) {
        .hue => std.mem.sort(Pixel, colours, {}, colours_lib.hueCompare),
        .hsp => std.mem.sort(Pixel, colours, {}, colours_lib.hspCompare),
        .zigzag => {}, // this is really for generation order
        .zorder => std.mem.sort(Pixel, colours, {}, colours_lib.zOrderCompare),
        .none => {},
    }

    const size_x, const size_y = imageSizeFromBitDepth(parameters.channel_depth);

    const buffer = try allocator.alloc(Pixel, colours.len);
    var image = ImageData{
        .buffer = buffer,
        .size_x = size_x,
        .size_y = size_y,
    };
    std.debug.print("Producing a {d}x{d} image\n", .{ size_x, size_y });
    fillImage(
        gen_alloc,
        colours,
        &image,
        parameters.starts,
        parameters.stride,
        parameters.tree_type,
        parameters.approximate,
        parameters.seed,
    ) catch |err| {
        std.debug.print("Error filling image: {any}\n", .{err});
    };

    if (parameters.verify) {
        if (!verifyFullImagePopulated(image)) {
            std.debug.print("Image verified\n", .{});
        }
        if (!try verifyUniquePixels(allocator, image)) {
            std.debug.print("Duplicate pixels detected\n", .{});
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
