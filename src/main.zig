const std = @import("std");
const png = @import("png.zig");
const kd_tree = @import("kd_tree.zig");
const Pixel = @import("colours.zig").Pixel;
const hspCompare = @import("colours.zig").hspCompare;
const hueCompare = @import("colours.zig").hueCompare;

const BIT_DEPTH = 8;
const CHANNEL_SIZE = 1 << BIT_DEPTH;
const DATA_SIZE = CHANNEL_SIZE * CHANNEL_SIZE * CHANNEL_SIZE;
const IMAGE_X = 4096;
const IMAGE_Y = 4096;

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

const ColoursIterator = struct {
    buffer: []const Pixel,
    start: usize = 0,
    count: usize = 0,
    len: usize = 0,

    pub fn next(it: *ColoursIterator) ?Pixel {
        if (it.count >= it.len) {
            return null;
        }
        const current = (it.start + it.count) % it.buffer.len;
        it.count += 1;
        return it.buffer[current];
    }

    pub fn iterateColours(colours: []const Pixel, start: usize) ColoursIterator {
        return .{
            .buffer = colours,
            .start = start,
            .count = 0,
            .len = colours.len,
        };
    }
};

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

fn createColours(allocator: std.mem.Allocator, zigzag: bool) ![]Pixel {
    var colours: []Pixel = try allocator.alloc(Pixel, DATA_SIZE);
    var idx: usize = 0;
    const step = 256 / CHANNEL_SIZE;

    var red_it: ChannelIterator = undefined;
    var green_it: ChannelIterator = undefined;
    var blue_it: ChannelIterator = undefined;
    if (zigzag) {
        red_it = .{ .reversing = .{ .count = CHANNEL_SIZE, .step = step } };
        green_it = .{ .reversing = .{ .count = CHANNEL_SIZE, .step = step } };
        blue_it = .{ .reversing = .{ .count = CHANNEL_SIZE, .step = step } };
    } else {
        red_it = .{ .repeat = .{ .count = CHANNEL_SIZE, .step = step } };
        green_it = .{ .repeat = .{ .count = CHANNEL_SIZE, .step = step } };
        blue_it = .{ .repeat = .{ .count = CHANNEL_SIZE, .step = step } };
    }

    for (0..CHANNEL_SIZE) |_| {
        const red = red_it.next();
        for (0..CHANNEL_SIZE) |_| {
            const green = green_it.next();
            for (0..CHANNEL_SIZE) |_| {
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
) !kd_tree.KdTree(Pixel, ImageCoord) {
    var tree: kd_tree.KdTree(Pixel, ImageCoord) = .{};
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
    var tree = kd_tree.KdTree(Pixel, ImageCoord){};

    const start_idx = rng.intRangeLessThan(usize, 0, colours.len);
    var colours_it = ColoursIterator.iterateColours(
        colours,
        start_idx,
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
    var c_count: usize = 0;
    while (colours_it.next()) |colour| {
        c_count += 1;
        if (c_count % (colours.len / 100) == 0) {
            std.debug.print("Progress: {d} - Tree leaves: {d}\n", .{ percentage, tree.leaf_count });
            percentage += 1;
        }
        const empties = tree.empty_leaf_count;
        const ratio: f32 = @as(f32, @floatFromInt(empties)) / @as(f32, @floatFromInt(tree.leaf_count));
        const rebuild = empties > 1 and ratio >= 0.05;
        if (rebuild) {
            std.debug.print("Rebuilding tree - leaves: {any} - empty {any}\n", .{ tree.leaf_count, tree.empty_leaf_count });
            _ = tree_alloc.reset(.retain_capacity);
            tree = try populateNewTree(tree_allocator, image.*);
        }

        var slot: ImageCoord = .{ .x = 0, .y = 0 };
        var pixel: Pixel = .{};
        while (true) {
            const closest_pixel, const closest_idx = tree.getNearest(
                colour,
            ) orelse return error.InvalidStateEmptyTree;
            var buffer = std.mem.zeroes([8]ImageCoord);
            const available = getAvailableNeighbours(image.*, closest_idx, &buffer);
            if (available.len == 0) {
                try tree.remove(closest_pixel);
                continue;
            }
            // we're about to remove this pixel's last available neighbour
            if (available.len == 1) {
                try tree.remove(closest_pixel);
            }
            const pick_idx = rng.intRangeLessThan(usize, 0, available.len);
            slot = available[pick_idx];
            pixel = closest_pixel;
            break;
        }
        image.put(slot, colour);
        if (getFirstNeighbour(image.*, slot) != null) {
            try tree.add(tree_allocator, colour, slot);
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

const ColourSort = enum {
    hue,
    hsp,
    zigzag,
    none,
};

const Parameters = struct {
    seed: u32 = 0,
    starts: u16 = 0,
    sort_type: ColourSort = ColourSort.hue,
    output_file: []const u8 = "",
};

fn parseArguments(args: std.process.Args) !Parameters {
    var seed: u32 = 0;
    var starts: u16 = 1;
    var output_file: []const u8 = "out.png";
    var sort_type = ColourSort.hue;
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
        } else if (std.mem.eql(u8, "--none", arg)) {
            sort_type = ColourSort.none;
        }
    }

    return .{
        .seed = seed,
        .starts = starts,
        .output_file = output_file,
        .sort_type = sort_type,
    };
}

pub fn main(init: std.process.Init) !void {
    var gpa = std.heap.DebugAllocator(.{}){};
    const gen_alloc = gpa.allocator();
    var arena = std.heap.ArenaAllocator.init(gen_alloc);
    defer arena.deinit();
    const allocator = arena.allocator();

    const parameters = try parseArguments(init.minimal.args);

    const colours = try createColours(allocator, parameters.sort_type == ColourSort.zigzag);

    switch (parameters.sort_type) {
        .hue => std.mem.sort(Pixel, colours, {}, hueCompare),
        .hsp => std.mem.sort(Pixel, colours, {}, hspCompare),
        .zigzag => {}, // this is really for generation order
        .none => {},
    }

    const buffer = try allocator.alloc(Pixel, DATA_SIZE);
    var image = ImageData{
        .buffer = buffer,
        .size_x = IMAGE_X,
        .size_y = IMAGE_Y,
    };
    fillImage(gen_alloc, colours, &image, parameters.starts, parameters.seed) catch |err| {
        std.debug.print("Error filling image: {any}\n", .{err});
    };

    if (!verifyFullImagePopulated(image)) {
        std.debug.print("Image verified\n", .{});
    }

    try png.writePng(
        allocator,
        init.io,
        buffer,
        parameters.output_file,
        IMAGE_X,
        IMAGE_Y,
    );
}
