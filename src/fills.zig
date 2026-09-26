const std = @import("std");

const ImageCoord = @import("image.zig").ImageCoord;
const ImageData = @import("image.zig").ImageData;
const Pixel = @import("colours.zig").Pixel;

pub const Error = error{
    MissingKey,
    TreeError,
    InvalidState,
};

pub const Fills = enum {
    min,
    target,
    mean,
};

pub const ImageFill = struct {
    ptr: *anyopaque,
    rng: std.Random,
    vtable: *const VTable,

    const VTable = struct {
        place: *const fn (*anyopaque, colour: Pixel, coord: ImageCoord) Error!void,
        placeRandomly: *const fn (*anyopaque, colour: Pixel) Error!void,
        matchAndPlace: *const fn (*anyopaque, colour: Pixel) Error!void,
        repopulate: *const fn (*anyopaque) Error!void,
        emptyLeaves: *const fn (*anyopaque) usize,
        leaves: *const fn (*anyopaque) usize,
    };

    pub fn place(self: *const @This(), colour: Pixel, coord: ImageCoord) Error!void {
        return self.vtable.place(self.ptr, colour, coord);
    }

    pub fn placeRandomly(self: *const @This(), colour: Pixel) Error!void {
        return self.vtable.placeRandomly(self.ptr, colour);
    }

    pub fn matchAndPlace(self: *@This(), colour: Pixel) Error!void {
        return self.vtable.matchAndPlace(self.ptr, colour);
    }

    pub fn repopulate(self: *@This()) Error!void {
        return self.vtable.repopulate(self.ptr);
    }

    pub fn leaves(self: *@This()) usize {
        return self.vtable.leaves(self.ptr);
    }

    pub fn emptyLeaves(self: *@This()) usize {
        return self.vtable.emptyLeaves(self.ptr);
    }
};

pub fn MinFill(comptime tree_type: type) type {
    return struct {
        allocator: std.heap.ArenaAllocator,
        image: *ImageData,
        tree: tree_type,
        wrap: bool,
        approximate: bool,
        rng: std.Random,

        pub fn init(
            allocator: std.mem.Allocator,
            image: *ImageData,
            wrap: bool,
            approximate: bool,
            rng: std.Random,
        ) @This() {
            const tree_arena = std.heap.ArenaAllocator.init(allocator);
            const tree: tree_type = tree_type.init();
            return .{
                .allocator = tree_arena,
                .image = image,
                .tree = tree,
                .wrap = wrap,
                .approximate = approximate,
                .rng = rng,
            };
        }

        pub fn place(
            ptr: *anyopaque,
            colour: Pixel,
            coord: ImageCoord,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            self.image.put(coord, colour);
            self.tree.add(
                self.allocator.allocator(),
                colour,
                coord,
            ) catch return Error.TreeError;
        }

        pub fn placeRandomly(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            const initial_x = self.rng.intRangeLessThan(usize, 0, self.image.size_x);
            const initial_y = self.rng.intRangeLessThan(usize, 0, self.image.size_y);
            try place(
                self,
                colour,
                .{ .x = initial_x, .y = initial_y },
            );
        }

        // Wipes and reconstructs the held tree using the current image state
        pub fn repopulate(
            ptr: *anyopaque,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            _ = self.allocator.reset(.retain_capacity);
            self.tree = tree_type.init();
            for (0..self.image.size_y) |y| {
                for (0..self.image.size_x) |x| {
                    const slot: ImageCoord = .{ .x = x, .y = y };
                    const pixel = self.image.at(slot);
                    if (pixel.alpha == 0) {
                        continue;
                    }
                    if (hasOpenNeighbours(self.image.*, slot, self.wrap)) {
                        self.tree.add(
                            self.allocator.allocator(),
                            pixel,
                            slot,
                        ) catch return Error.TreeError;
                    }
                }
            }
        }

        pub fn matchAndPlace(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            const search_fn = if (self.approximate)
                &tree_type.getNear
            else
                &tree_type.getNearest;
            while (true) {
                const closest_pixel, const closest_idx = search_fn(
                    &self.tree,
                    colour,
                ) orelse return Error.InvalidState;
                var buffer = std.mem.zeroes([8]ImageCoord);
                const available = getOpenNeighbours(
                    self.image.*,
                    closest_idx,
                    &buffer,
                    self.wrap,
                );

                // we found a pixel that should have been removed but never was
                if (available.len == 0) {
                    self.tree.remove(closest_pixel, closest_idx) catch return Error.MissingKey;
                    continue;
                }
                // we're about to remove this pixel's last available neighbour
                if (available.len == 1) {
                    self.tree.remove(closest_pixel, closest_idx) catch return Error.MissingKey;
                }

                const pick_idx = self.rng.intRangeLessThan(
                    usize,
                    0,
                    available.len,
                );
                const slot = available[pick_idx];
                self.image.put(slot, colour);
                if (hasOpenNeighbours(self.image.*, slot, self.wrap)) {
                    self.tree.add(self.allocator.allocator(), colour, slot) catch return Error.TreeError;
                }
                break;
            }
        }

        pub fn leaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.leaves();
        }

        pub fn emptyLeaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.emptyLeaves();
        }

        pub fn filler(self: *@This()) ImageFill {
            return .{
                .ptr = self,
                .rng = self.rng,
                .vtable = &.{
                    .place = &@This().place,
                    .placeRandomly = &@This().placeRandomly,
                    .matchAndPlace = &@This().matchAndPlace,
                    .repopulate = &@This().repopulate,
                    .leaves = &@This().leaves,
                    .emptyLeaves = &@This().emptyLeaves,
                },
            };
        }
    };
}

/// A fill strategy that tries to place pixels according
/// to a reference image.
pub fn TargetFill(comptime tree_type: type) type {
    return struct {
        allocator: std.heap.ArenaAllocator,
        image: *ImageData,
        reference: *ImageData,
        tree: tree_type,
        approximate: bool,
        rng: std.Random,

        pub fn init(
            allocator: std.mem.Allocator,
            image: *ImageData,
            reference: *ImageData,
            approximate: bool,
            rng: std.Random,
        ) @This() {
            const tree_arena = std.heap.ArenaAllocator.init(allocator);
            const tree: tree_type = tree_type.init();
            return .{
                .allocator = tree_arena,
                .image = image,
                .reference = reference,
                .tree = tree,
                .approximate = approximate,
                .rng = rng,
            };
        }

        pub fn place(
            ptr: *anyopaque,
            colour: Pixel,
            coord: ImageCoord,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            self.image.put(coord, colour);
            const value = self.reference.at(coord);
            try self.tree.remove(
                value,
                coord,
            );
        }

        pub fn placeRandomly(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            return matchAndPlace(
                ptr,
                colour,
            );
        }

        pub fn matchAndPlace(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            const search_fn = if (self.approximate)
                &tree_type.getNear
            else
                &tree_type.getNearest;

            const closest_pixel, const closest_idx = search_fn(
                &self.tree,
                colour,
            ) orelse return Error.InvalidState;

            self.image.put(closest_idx, colour);
            self.tree.remove(closest_pixel, closest_idx) catch return Error.MissingKey;
        }

        pub fn repopulate(
            ptr: *anyopaque,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            _ = self.allocator.reset(.retain_capacity);
            self.tree = tree_type.init();

            for (0..self.reference.size_y) |y| {
                for (0..self.reference.size_x) |x| {
                    const slot: ImageCoord = .{ .x = x, .y = y };
                    if (self.image.at(slot).alpha == 0) {
                        self.tree.add(
                            self.allocator.allocator(),
                            self.reference.at(slot),
                            slot,
                        ) catch return Error.TreeError;
                    }
                }
            }
        }

        pub fn leaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.leaves();
        }

        pub fn emptyLeaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.emptyLeaves();
        }

        pub fn filler(self: *@This()) ImageFill {
            return .{
                .ptr = self,
                .rng = self.rng,
                .vtable = &.{
                    .place = &@This().place,
                    .placeRandomly = &@This().placeRandomly,
                    .matchAndPlace = &@This().matchAndPlace,
                    .repopulate = &@This().repopulate,
                    .leaves = &@This().leaves,
                    .emptyLeaves = &@This().emptyLeaves,
                },
            };
        }
    };
}

const MeanVal = struct {
    red: usize = 0,
    green: usize = 0,
    blue: usize = 0,
    count: usize = 0,

    pub fn add(self: *const MeanVal, other: Pixel) MeanVal {
        return .{
            .red = self.red + other.red,
            .green = self.green + other.green,
            .blue = self.blue + other.blue,
            .count = self.count + 1,
        };
    }

    pub fn toPixel(self: *const MeanVal) Pixel {
        const red = self.red / self.count;
        const green = self.green / self.count;
        const blue = self.blue / self.count;
        return .{
            .red = @intCast(red),
            .green = @intCast(green),
            .blue = @intCast(blue),
            .alpha = 0xFF,
        };
    }
};

/// Match according the mean of already populate neighbours
pub fn MeanFill(comptime tree_type: type) type {
    return struct {
        allocator: std.heap.ArenaAllocator,
        image: *ImageData,
        means: []MeanVal,
        tree: tree_type,
        wrap: bool,
        approximate: bool,
        rng: std.Random,

        pub fn init(
            allocator: std.mem.Allocator,
            image: *ImageData,
            wrap: bool,
            approximate: bool,
            rng: std.Random,
        ) !@This() {
            var means = try allocator.alloc(MeanVal, image.buffer.len);
            for (0..means.len) |idx| {
                means[idx] = .{};
            }
            const tree_arena = std.heap.ArenaAllocator.init(allocator);
            const tree: tree_type = tree_type.init();
            return .{
                .allocator = tree_arena,
                .image = image,
                .means = means,
                .tree = tree,
                .wrap = wrap,
                .approximate = approximate,
                .rng = rng,
            };
        }

        pub fn repopulate(
            ptr: *anyopaque,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            _ = self.allocator.reset(.retain_capacity);
            self.tree = tree_type.init();

            for (0..self.means.len) |mean_idx| {
                const mean_val = self.means[mean_idx];
                if (mean_val.count == 0) {
                    continue;
                }
                const mean_coord = self.image.fromIndex(mean_idx);
                self.tree.add(
                    self.allocator.allocator(),
                    mean_val.toPixel(),
                    mean_coord,
                ) catch return Error.TreeError;
            }
        }

        pub fn place(
            ptr: *anyopaque,
            colour: Pixel,
            coord: ImageCoord,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            self.image.put(coord, colour);

            // remove the mean value from the tree (if it's set)
            const to_remove = self.means[self.image.toIndex(coord)];
            if (to_remove.count > 0) {
                try self.tree.remove(to_remove.toPixel(), coord);
                // Erase it from our mean set by setting its count to 0
                self.means[self.image.toIndex(coord)] = .{};
            }

            // Update the means of all our neighbours
            // The maths should work out for neighbours with a mean count of 0
            var buffer = std.mem.zeroes([8]ImageCoord);
            const neighbours = getOpenNeighbours(self.image.*, coord, &buffer, self.wrap);
            for (0..neighbours.len) |idx| {
                const neighbour = neighbours[idx];
                const linear_idx = self.image.toIndex(neighbour);
                const old_val = self.means[linear_idx];
                const new_val = old_val.add(colour);
                self.means[linear_idx] = new_val;
                if (old_val.count != 0) {
                    try self.tree.remove(old_val.toPixel(), neighbour);
                }
                self.tree.add(
                    self.allocator.allocator(),
                    new_val.toPixel(),
                    neighbour,
                ) catch return Error.TreeError;
            }
        }

        pub fn matchAndPlace(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            // Mean filling requires a multi tree and we want to randomly
            // select a matching coordinate to avoid biasing the pixel
            // placement in a particular direction.
            const search_fn = if (self.approximate)
                &tree_type.getAllNear
            else
                &tree_type.getAllNearest;

            _, const indices = search_fn(
                &self.tree,
                colour,
            ) orelse return Error.InvalidState;

            const selected_index = self.rng.intRangeLessThan(usize, 0, indices.len);
            const closest_idx = indices[selected_index];

            try place(self, colour, closest_idx);
        }

        pub fn placeRandomly(
            ptr: *anyopaque,
            colour: Pixel,
        ) Error!void {
            const self: *@This() = @ptrCast(@alignCast(ptr));

            const initial_x = self.rng.intRangeLessThan(usize, 0, self.image.size_x);
            const initial_y = self.rng.intRangeLessThan(usize, 0, self.image.size_y);
            try place(
                self,
                colour,
                .{ .x = initial_x, .y = initial_y },
            );
        }

        pub fn leaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.leaves();
        }

        pub fn emptyLeaves(ptr: *anyopaque) usize {
            const self: *@This() = @ptrCast(@alignCast(ptr));
            return self.tree.emptyLeaves();
        }

        pub fn filler(self: *@This()) ImageFill {
            return .{
                .ptr = self,
                .rng = self.rng,
                .vtable = &.{
                    .place = &@This().place,
                    .placeRandomly = &@This().placeRandomly,
                    .matchAndPlace = &@This().matchAndPlace,
                    .repopulate = &@This().repopulate,
                    .leaves = &@This().leaves,
                    .emptyLeaves = &@This().emptyLeaves,
                },
            };
        }
    };
}

fn getNeighboursWrapped(
    image: ImageData,
    point: ImageCoord,
    buffer: []ImageCoord,
) []ImageCoord {
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
            buffer[c_idx] = .{
                .x = @as(usize, @intCast(x)),
                .y = @as(usize, @intCast(y)),
            };
            c_idx += 1;
        }
    }
    return buffer[0..c_idx];
}

fn getNeighboursClamped(
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

    if (!top_edge and !left_edge) {
        buffer[idx] = .{ .x = bounded.x - 1, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!top_edge) {
        buffer[idx] = .{ .x = bounded.x, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!top_edge and !right_edge) {
        buffer[idx] = .{ .x = bounded.x + 1, .y = bounded.y - 1 };
        idx += 1;
    }
    if (!left_edge) {
        buffer[idx] = .{ .x = bounded.x - 1, .y = bounded.y };
        idx += 1;
    }
    if (!right_edge) {
        buffer[idx] = .{ .x = bounded.x + 1, .y = bounded.y };
        idx += 1;
    }
    if (!left_edge and !bottom_edge) {
        buffer[idx] = .{ .x = bounded.x - 1, .y = bounded.y + 1 };
        idx += 1;
    }
    if (!bottom_edge) {
        buffer[idx] = .{ .x = bounded.x, .y = bounded.y + 1 };
        idx += 1;
    }
    if (!right_edge and !bottom_edge) {
        buffer[idx] = .{ .x = bounded.x + 1, .y = bounded.y + 1 };
        idx += 1;
    }
    return buffer[0..idx];
}

// Fills the buffer with the unfilled neighbours of the given pixel.
fn getOpenNeighbours(
    image: ImageData,
    point: ImageCoord,
    buffer: []ImageCoord,
    wrap: bool,
) []ImageCoord {
    var to_fill = std.mem.zeroes([8]ImageCoord);
    var candidates: []ImageCoord = undefined;
    if (wrap) {
        candidates = getNeighboursWrapped(image, point, &to_fill);
    } else {
        candidates = getNeighboursClamped(image, point, &to_fill);
    }

    var out_idx: usize = 0;
    for (0..candidates.len) |candidate_idx| {
        const test_point = candidates[candidate_idx];
        if (image.at(test_point).alpha == 0) {
            buffer[out_idx] = test_point;
            out_idx += 1;
        }
    }
    return buffer[0..out_idx];
}

fn hasOpenNeighbours(image: ImageData, point: ImageCoord, wrap: bool) bool {
    var buffer: [8]ImageCoord = std.mem.zeroes([8]ImageCoord);
    const available = getOpenNeighbours(image, point, &buffer, wrap);
    return available.len > 0;
}
