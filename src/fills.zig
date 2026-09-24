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

// Returns a subslice of the provided buffer populate with open neighbours the
// given point. This accounts for wrapping at the image boundaries.
fn getOpenNeighboursWrapped(
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
fn getOpenNeighboursClamped(
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

// Fills the buffer with the unfilled neighbours of the given pixel.
fn getOpenNeighbours(
    image: ImageData,
    point: ImageCoord,
    buffer: []ImageCoord,
    wrap: bool,
) []ImageCoord {
    if (wrap) {
        return getOpenNeighboursWrapped(image, point, buffer);
    } else {
        return getOpenNeighboursClamped(image, point, buffer);
    }
}

fn hasOpenNeighbours(image: ImageData, point: ImageCoord, wrap: bool) bool {
    var buffer: [8]ImageCoord = std.mem.zeroes([8]ImageCoord);
    const available = getOpenNeighbours(image, point, &buffer, wrap);
    return available.len > 0;
}
