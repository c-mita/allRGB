const std = @import("std");
const colours = @import("colours.zig");

const Bucket = @import("bucket.zig").Bucket;

fn KdTreeLeaf(comptime K: type, comptime V: type) type {
    return Bucket(K, V, K.distance);
}

fn KdTreeNode(comptime K: type, comptime V: type) type {
    return struct {
        left: ?*KdTreeNode(K, V) = null,
        right: ?*KdTreeNode(K, V) = null,
        leaf: ?*KdTreeLeaf(K, V) = null,
        // TODO: The split information should be a normal vector
        split_axis: u32 = 0,
        split_value: u8 = 0,

        /// Create a new tree node that contains an empty leaf.
        fn initLeafNode(allocator: std.mem.Allocator) !*KdTreeNode(K, V) {
            const leaf = try allocator.create(KdTreeLeaf(K, V));
            leaf.* = .{};
            const node = try allocator.create(KdTreeNode(K, V));
            node.* = .{ .leaf = leaf };
            return node;
        }

        fn getNodeForKey(self: *KdTreeNode(K, V), key: K) *KdTreeNode(K, V) {
            var node = self;
            while (node.leaf == null) {
                if (key.lessThan(node.split_axis, node.split_value)) {
                    node = node.left orelse unreachable;
                } else {
                    node = node.right orelse unreachable;
                }
            }
            return node;
        }

        /// Puts a value into the tree.
        /// Overwrites the previous value if the key was already added.
        fn add(self: *KdTreeNode(K, V), allocator: std.mem.Allocator, key: K, value: V) !void {
            var node = self.getNodeForKey(key);
            var leaf = node.leaf.?;
            if (leaf.count == leaf.keys.len) {
                try node.split(allocator);
                node = node.getNodeForKey(key);
                leaf = node.leaf.?;
            }
            leaf.putValue(key, value) catch {
                unreachable;
            };
        }

        /// Removes a key from the tree.
        fn remove(self: *KdTreeNode(K, V), key: K) !void {
            const leaf = self.getNodeForKey(key);
            try leaf.leaf.?.removeKey(key);
        }

        fn split(self: *KdTreeNode(K, V), allocator: std.mem.Allocator) !void {
            const leaf = self.leaf orelse return error.KdTreeUnsplittable;
            if (self.left != null or self.right != null) {
                return error.KdTreeUnsplittable;
            }
            const axis, const split_point = K.splitElements(leaf.keys[0..leaf.count]);
            const left_leaf_ptr = try allocator.create(KdTreeLeaf(K, V));
            errdefer allocator.destroy(left_leaf_ptr);
            const right_leaf_ptr = try allocator.create(KdTreeLeaf(K, V));
            errdefer allocator.destroy(right_leaf_ptr);
            left_leaf_ptr.* = .{};
            right_leaf_ptr.* = .{};
            for (0..leaf.count) |idx| {
                const key = leaf.keys[idx];
                const value = leaf.values[idx];
                if (key.lessThan(axis, split_point)) {
                    try left_leaf_ptr.*.simplePut(key, value);
                } else {
                    try right_leaf_ptr.*.simplePut(key, value);
                }
            }
            if (left_leaf_ptr.*.count == 0 or right_leaf_ptr.*.count == 0) {
                std.debug.print(
                    "Bad split on axis {any} for value {any}\n",
                    .{ axis, split_point },
                );
                std.debug.print("To split: {any}\n", .{leaf.keys});
                return error.KdTreeBadSplit;
            }
            const left = try allocator.create(KdTreeNode(K, V));
            errdefer allocator.destroy(left);
            const right = try allocator.create(KdTreeNode(K, V));
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
            self.split_axis = axis;
            self.split_value = split_point;
        }

        fn getNearest(self: *const KdTreeNode(K, V), key: K) ?struct { K, V } {
            if (self.leaf != null) {
                return self.leaf.?.getNearest(key);
            }
            const distance_to_split = key.distanceToSplit(self.split_axis, self.split_value);
            const dsquared = distance_to_split * distance_to_split;
            const first = if (distance_to_split < 0) self.left.? else self.right.?;
            const second = if (distance_to_split < 0) self.right.? else self.left.?;

            const first_result = first.getNearest(key);
            if (first_result != null) {
                const first_candidate, const first_value = first_result.?;
                const first_distance = key.distanceSquared(first_candidate);
                if (first_distance <= dsquared) {
                    return .{ first_candidate, first_value };
                }
            }
            const second_result = second.getNearest(key);
            if (second_result == null) {
                return first_result;
            } else if (first_result == null) {
                return second_result;
            }
            const first_candidate, const first_value = first_result.?;
            const first_distance = key.distanceSquared(first_candidate);
            const second_candidate, const second_value = second_result.?;
            const second_distance = key.distanceSquared(second_candidate);
            if (second_distance < first_distance) {
                return .{ second_candidate, second_value };
            } else {
                return .{ first_candidate, first_value };
            }
        }

        /// Only every queries the "other" side of the spit if there are no
        /// keys on the same side as the split as the query node.
        fn getNear(self: *const KdTreeNode(K, V), key: K) ?struct { K, V } {
            if (self.leaf != null) {
                return self.leaf.?.getNearest(key);
            }
            const distance_to_split = key.distanceToSplit(self.split_axis, self.split_value);
            const first = if (distance_to_split < 0) self.left.? else self.right.?;
            const second = if (distance_to_split < 0) self.right.? else self.left.?;

            return first.getNear(key) orelse second.getNear(key);
        }

        fn clear(self: *KdTreeNode(K, V), allocator: std.mem.Allocator) void {
            if (self.leaf != null) {
                allocator.destroy(self.leaf.?);
            }
            if (self.left != null) {
                self.left.?.clear(allocator);
                allocator.destroy(self.left.?);
            }
            if (self.right != null) {
                self.right.?.clear(allocator);
                allocator.destroy(self.right.?);
            }
        }
    };
}

pub fn KdTree(comptime K: type, comptime V: type) type {
    return struct {
        root: ?*KdTreeNode(K, V) = null,
        leaf_count: usize = 0,
        empty_leaf_count: usize = 0,

        /// Returns the closest match to the input key in the tree.
        pub fn getNearest(self: *KdTree(K, V), key: K) ?struct { K, V } {
            return if (self.root != null) self.root.?.getNearest(key) else null;
        }

        /// Returns a potentially "near" match to the input key. If a match
        /// is found in the same "subspace" as the key then return it (only if
        /// there is no match on the same side of a splitting plane as the key
        /// will the other side be searched).
        pub fn getNear(self: *KdTree(K, V), key: K) ?struct { K, V } {
            return if (self.root != null) self.root.?.getNear(key) else null;
        }

        /// Adds a key value pair to the tree
        pub fn add(self: *KdTree(K, V), allocator: std.mem.Allocator, key: K, value: V) !void {
            if (self.root == null) {
                self.root = try KdTreeNode(K, V).initLeafNode(allocator);
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
                // this means the node split and we have an extra leaf
                self.leaf_count += 1;
            }
        }

        pub fn remove(self: *KdTree(K, V), key: K) !void {
            if (self.root == null) {
                return error.KdTreeKeyNotFound;
            }
            const node = self.root.?.getNodeForKey(key);
            try node.leaf.?.removeKey(key);
            if (node.leaf.?.count == 0) {
                self.empty_leaf_count += 1;
            }
        }

        pub fn clear(self: *KdTree(K, V), allocator: std.mem.Allocator) void {
            self.leaf_count = 0;
            self.empty_leaf_count = 0;
            if (self.root == null) {
                return;
            }
            self.root.?.clear(allocator);
        }
    };
}

test "KdTreeNode lookup found on left" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const left_leaf = try KdTreeNode(colours.Pixel, i32).initLeafNode(allocator);
    _ = try left_leaf.*.add(allocator, .{ .green = 10 }, 1);
    _ = try left_leaf.*.add(allocator, .{ .red = 20 }, 2);
    const right_leaf = try KdTreeNode(colours.Pixel, i32).initLeafNode(allocator);
    _ = try right_leaf.*.add(allocator, .{ .blue = 100 }, 3);

    const node = KdTreeNode(colours.Pixel, i32){
        .left = left_leaf,
        .right = right_leaf,
        .leaf = null,
        .split_axis = 2,
        .split_value = 30,
    };

    try std.testing.expectEqual(
        .{ colours.Pixel{ .red = 20 }, 2 },
        node.getNearest(.{ .red = 30, .green = 11, .blue = 1 }),
    );
    try std.testing.expectEqual(
        .{ colours.Pixel{ .green = 10 }, 1 },
        node.getNearest(.{ .red = 1, .green = 10, .blue = 29 }),
    );
}

test "KdTreeNode lookup found on right" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const left_leaf = try KdTreeNode(colours.Pixel, i32).initLeafNode(allocator);
    _ = try left_leaf.*.add(allocator, .{ .green = 10 }, 1);
    _ = try left_leaf.*.add(allocator, .{ .red = 20 }, 2);
    const right_leaf = try KdTreeNode(colours.Pixel, i32).initLeafNode(allocator);
    _ = try right_leaf.*.add(allocator, .{ .blue = 100 }, 3);

    const node = KdTreeNode(colours.Pixel, i32){
        .left = left_leaf,
        .right = right_leaf,
        .leaf = null,
        .split_axis = 2,
        .split_value = 30,
    };

    try std.testing.expectEqual(
        .{ colours.Pixel{ .blue = 100 }, 3 },
        node.getNearest(.{ .red = 30, .green = 11, .blue = 190 }),
    );
}

test "KdNode split" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    const leaf = try allocator.create(KdTreeLeaf(colours.Pixel, i32));
    leaf.* = .{};
    for (0..16) |idx| {
        const v: u8 = @intCast(idx);
        try leaf.putValue(.{ .green = v * 2 }, @intCast(idx));
    }
    try std.testing.expectEqual(16, leaf.count);

    var node = KdTreeNode(colours.Pixel, i32){
        .left = null,
        .right = null,
        .leaf = leaf,
    };
    try node.split(allocator);

    try std.testing.expectEqual(
        8,
        node.left.?.leaf.?.count,
    );
    try std.testing.expectEqual(
        8,
        node.right.?.leaf.?.count,
    );
}

test "KdNode split on full" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    var node = try KdTreeNode(colours.Pixel, usize).initLeafNode(allocator);

    const leaf_size = node.leaf.?.values.len;
    for (0..leaf_size * 2) |idx| {
        const v: u8 = @intCast(idx);
        _ = try node.add(allocator, .{ .green = v }, @intCast(idx));
    }
    // Because our pixels are added in ascending order we expect a tree
    // biased to the right
    const half_leaf = leaf_size / 2;
    try std.testing.expectEqual(half_leaf, node.left.?.leaf.?.count);
    try std.testing.expectEqual(0, node.left.?.leaf.?.values[0]);
    try std.testing.expectEqual(half_leaf, node.right.?.left.?.leaf.?.count);
    try std.testing.expectEqual(half_leaf, node.right.?.left.?.leaf.?.values[0]);
    try std.testing.expectEqual(leaf_size - 1, node.right.?.left.?.leaf.?.values[half_leaf - 1]);
    try std.testing.expectEqual(leaf_size, node.right.?.right.?.leaf.?.count);
    try std.testing.expectEqual(leaf_size, node.right.?.right.?.leaf.?.values[0]);
    try std.testing.expectEqual(leaf_size * 2 - 1, node.right.?.right.?.leaf.?.values[leaf_size - 1]);
}

test "KdTree lookup" {
    const test_alloc = std.testing.allocator;
    var arena_alloc = std.heap.ArenaAllocator.init(test_alloc);
    defer arena_alloc.deinit();
    const allocator = arena_alloc.allocator();

    var tree = KdTree(colours.Pixel, i32){};
    for (0..1000) |idx| {
        const red: u8 = @intCast(idx % 256);
        const green: u8 = @intCast((idx * 3) % 256);
        const blue: u8 = @intCast((idx * 5) % 256);
        const px = colours.Pixel{ .red = red, .green = green, .blue = blue };
        try tree.add(allocator, px, @intCast(idx));
    }

    try std.testing.expectEqual(
        // 772 % 256 == 4 so this should element should match the fourth thing we add
        .{ colours.Pixel{ .red = 4, .green = 12, .blue = 20 }, 772 },
        tree.getNearest(.{ .red = 4, .green = 12, .blue = 20 }),
    );
    try std.testing.expectEqual(
        .{ colours.Pixel{ .red = 4, .green = 12, .blue = 20 }, 772 },
        tree.getNearest(.{ .red = 5, .green = 11, .blue = 18 }),
    );
}
