const std = @import("std");

const Bucket = @import("bucket.zig").Bucket;

fn VpTreeLeaf(
    comptime K: type,
    comptime V: type,
    comptime distanceFunc: fn (K, K) usize,
) type {
    return Bucket(K, V, distanceFunc);
}

fn VpTreeNode(
    comptime K: type,
    comptime V: type,
    comptime distanceFunc: fn (K, K) usize,
) type {
    return struct {
        leaf: ?*VpTreeLeaf(K, V, distanceFunc) = null,
        left: ?*@This() = null,
        right: ?*@This() = null,
        centre: K = undefined,
        radius: usize = 0,

        fn initLeaf(allocator: std.mem.Allocator) !*@This() {
            const leaf = try allocator.create(VpTreeLeaf(K, V, distanceFunc));
            leaf.* = .{};
            const node = try allocator.create(@This());
            node.* = .{ .leaf = leaf };
            return node;
        }

        fn getNodeForKey(self: *@This(), key: K) *@This() {
            var node = self;
            while (node.leaf == null) {
                const distance = distanceFunc(key, node.centre);
                if (distance <= node.radius) {
                    node = node.left orelse unreachable;
                } else {
                    node = node.right orelse unreachable;
                }
            }
            return node;
        }

        fn add(self: *@This(), allocator: std.mem.Allocator, key: K, value: V) !void {
            var node = self.getNodeForKey(key);
            var leaf = node.leaf.?;
            if (leaf.count == leaf.keys.len) {
                try node.split(allocator);
                node = node.getNodeForKey(key);
                leaf = node.leaf.?;
            }
            leaf.putValue(key, value) catch unreachable;
        }

        fn remove(self: *@This(), key: K) !void {
            const leaf = self.getNodeForKey(key).leaf orelse unreachable;
            try leaf.removeKey(key);
        }

        fn getNearest(self: *@This(), key: K) ?struct { K, V } {
            if (self.leaf != null) {
                return self.leaf.?.getNearest(key);
            }

            const key_distance = distanceFunc(self.centre, key);
            const inside = key_distance <= self.radius;
            const first = if (inside) self.left.? else self.right.?;
            const second = if (inside) self.right.? else self.left.?;

            const first_candidate, const first_value = first.getNearest(key) orelse return second.getNearest(key);

            const found_distance = distanceFunc(first_candidate, key);
            // if we're inside the ball then we want to to check that
            // d(key, centre) + d(key, found) <= ball_radius
            //
            // otherwise we want to check that
            // ball_radius + d(key, found) <= d(key, centre)
            //
            // note the <= in the latter case is ok because even though it may
            // mean we technically intersect the closed ball, any point we may
            // find there will have the same distance
            const allowed_distance = if (inside) self.radius else key_distance;
            const actual_distance = if (inside) key_distance + found_distance else self.radius + found_distance;

            if (actual_distance <= allowed_distance) {
                return .{ first_candidate, first_value };
            }

            const second_result = second.getNearest(key);
            if (second_result == null) {
                return .{ first_candidate, first_value };
            }
            const second_candidate, const second_value = second_result.?;
            const second_found_distance = distanceFunc(second_candidate, key);
            if (second_found_distance < found_distance) {
                return .{ second_candidate, second_value };
            }
            return .{ first_candidate, first_value };
        }

        fn getNear(self: *@This(), key: K) ?struct { K, V } {
            if (self.leaf != null) {
                return self.leaf.?.getNearest(key);
            }

            const distance = distanceFunc(self.centre, key);
            const first = if (distance < self.radius) self.left.? else self.right.?;
            const second = if (distance < self.radius) self.right.? else self.left.?;

            return first.getNear(key) orelse second.getNear(key);
        }

        fn distanceSort(origin: K, lhs: K, rhs: K) bool {
            const left = distanceFunc(origin, lhs);
            const right = distanceFunc(origin, rhs);
            return left < right;
        }

        fn split(self: *@This(), allocator: std.mem.Allocator) !void {
            if (self.left != null or self.right != null) {
                return error.VpTreeUnsplittable;
            }

            const leaf = self.leaf orelse return error.VpTreeUnsplittable;
            if (leaf.count < 2) {
                return error.VpTreeUnsplittable;
            }

            const leaf_size = leaf.keys.len;

            var keys: [leaf_size]K = undefined;
            @memcpy(&keys, &leaf.keys);
            const origin = keys[0];
            std.mem.sortUnstable(K, keys[0..leaf.count], origin, @This().distanceSort);

            const left_leaf_ptr = try allocator.create(
                VpTreeLeaf(K, V, distanceFunc),
            );
            errdefer allocator.destroy(left_leaf_ptr);
            const right_leaf_ptr = try allocator.create(
                VpTreeLeaf(K, V, distanceFunc),
            );
            errdefer allocator.destroy(right_leaf_ptr);
            left_leaf_ptr.* = .{};
            right_leaf_ptr.* = .{};

            var mid_point_offset: usize = 0;
            var valid_split = false;
            var radius: usize = 0;
            // There is a chance that the "median" distance doesn't split anything
            // For instance, over half the other keys are exactly the same distance
            // away from our selected origin. So we need to trim back the radius
            // until we actually split something.
            while (!valid_split) {
                // the worse case (excluding the case of this leaf being
                // stuffed with identical keys) is that we split so the left
                // side contains only the origin point and the right contains
                // everything else.
                const mid_point = keys[leaf.count / 2 - mid_point_offset];
                radius = distanceFunc(origin, mid_point);
                for (0..leaf.count) |idx| {
                    const key = leaf.keys[idx];
                    const value = leaf.values[idx];
                    const distance = distanceFunc(origin, key);
                    if (distance <= radius) {
                        left_leaf_ptr.*.simplePut(key, value) catch unreachable;
                    } else {
                        right_leaf_ptr.*.simplePut(key, value) catch unreachable;
                    }
                }
                if (left_leaf_ptr.*.count == leaf_size or right_leaf_ptr.*.count == leaf_size) {
                    left_leaf_ptr.*.count = 0;
                    right_leaf_ptr.*.count = 0;
                    mid_point_offset += 1;
                } else {
                    valid_split = true;
                }
            }

            const left = try allocator.create(@This());
            errdefer allocator.destroy(left);
            const right = try allocator.create(@This());
            errdefer allocator.destroy(right);
            left.* = .{
                .leaf = left_leaf_ptr,
            };
            right.* = .{
                .leaf = right_leaf_ptr,
            };

            allocator.destroy(self.leaf.?);
            self.leaf = null;
            self.left = left;
            self.right = right;
            self.centre = origin;
            self.radius = radius;
        }
    };
}

/// A vantage-point tree to store key-value pairs.
/// Permits a nearest-neighbour look using the specified metric.
/// Note that the provided distance function must satisfy the
/// triangle inequality.
///
/// More general than a kd-tree (or any BSP tree, which require something
/// like an inner-product vector space) since this can work on any metric
/// space. However the efficiency is not quite as good.
pub fn VpTree(
    comptime K: type,
    comptime V: type,
    comptime distanceFunc: fn (K, K) usize,
) type {
    return struct {
        root: ?*VpTreeNode(K, V, distanceFunc) = null,
        leaf_count: usize = 0,
        empty_leaf_count: usize = 0,

        pub fn getNearest(self: *@This(), key: K) ?struct { K, V } {
            const root = self.root orelse return null;
            return root.getNearest(key);
        }

        pub fn getNear(self: *@This(), key: K) ?struct { K, V } {
            const root = self.root orelse return null;
            return root.getNear(key);
        }

        pub fn add(self: *@This(), allocator: std.mem.Allocator, key: K, value: V) !void {
            if (self.root == null) {
                self.root = try VpTreeNode(K, V, distanceFunc).initLeaf(allocator);
                self.leaf_count = 1;
                self.empty_leaf_count = 1;
            }
            const node = self.root.?.getNodeForKey(key);
            const leaf = node.leaf orelse unreachable;
            if (leaf.count == 0) {
                self.empty_leaf_count -= 1;
            }
            try node.add(allocator, key, value);
            if (node.leaf == null) {
                // the node had to be split, giving us a new leaf
                self.leaf_count += 1;
            }
        }

        pub fn remove(self: *@This(), key: K) !void {
            if (self.root == null) {
                return error.VpTreeKeyNotFound;
            }
            const node = self.root.?.getNodeForKey(key);
            try node.leaf.?.removeKey(key);
            if (node.leaf.?.count == 0) {
                self.empty_leaf_count += 1;
            }
        }
    };
}

fn intDistance(x: i32, y: i32) usize {
    return @abs(y - x);
}

test "VpTree add" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    var tree = VpTree(i32, usize, intDistance){};

    try tree.add(allocator, 0, 0);
    try tree.add(allocator, 16, 1024);
    try tree.add(allocator, 12, 2048);

    try std.testing.expectEqual(1, tree.leaf_count);
    try std.testing.expectEqual(0, tree.empty_leaf_count);

    const leaf = tree.root.?.leaf.?;
    try std.testing.expectEqual(3, leaf.count);
    try std.testing.expectEqualSlices(
        i32,
        &.{ 0, 16, 12 },
        leaf.keys[0..leaf.count],
    );
    try std.testing.expectEqualSlices(
        usize,
        &.{ 0, 1024, 2048 },
        leaf.values[0..leaf.count],
    );
}

test "VpTree lookup" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    var tree = VpTree(i32, usize, intDistance){};
    for (0..1024) |idx| {
        const k: i32 = @as(i32, @intCast(idx));
        const key = @mod(k * 8 - 128, 2048);
        try tree.add(allocator, key, idx);
    }
    try std.testing.expectEqual(0, tree.empty_leaf_count);
    try std.testing.expectEqual(true, tree.leaf_count > 0);

    // (1000 * 8 - 128) % 2048 == 1728
    try std.testing.expectEqual(.{ 1728, 1000 }, tree.getNearest(1728));
    try std.testing.expectEqual(.{ 1728, 1000 }, tree.getNearest(1730));
    try std.testing.expectEqual(.{ 1728, 1000 }, tree.getNearest(1732));
    try std.testing.expectEqual(.{ 1736, 1001 }, tree.getNearest(1733));
}
