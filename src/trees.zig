const std = @import("std");
const kd_tree = @import("kd_tree.zig");
const vp_tree = @import("vp_tree.zig");
const colours_lib = @import("colours.zig");

const Pixel = @import("colours.zig").Pixel;
const ImageCoord = @import("image.zig").ImageCoord;
const ImageData = @import("image.zig").ImageData;

pub const TreeType = enum {
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

pub fn TreeStore(comptime tree_type: TreeType) type {
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

/// Wrap a tree to handle the caes of duplicated keys.
pub fn TreeMultiStore(comptime tree_type: TreeType) type {
    return struct {
        tree: TreeType.treeTypeValue(tree_type, std.ArrayList(ImageCoord)),
        keys_count: usize = 0,
        values_count: usize = 0,

        pub fn init() @This() {
            return .{ .tree = .{}, .keys_count = 0, .values_count = 0 };
        }

        pub fn getNearest(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const items = self.getAllNearest(key) orelse return null;
            return .{ found, items[0] };
        }

        pub fn getNear(self: *@This(), key: Pixel) ?struct { Pixel, ImageCoord } {
            const found, const items = self.getAllNear(key) orelse return null;
            return .{ found, items[0] };
        }

        pub fn getAllNearest(self: *@This(), key: Pixel) ?struct { Pixel, []ImageCoord } {
            const found, const data = self.tree.getNearest(key) orelse return null;
            if (data.items.len == 0) {
                std.debug.print("SHOULD NEVER HAPPEN\n", .{});
                return null;
            }
            return .{ found, data.items };
        }

        pub fn getAllNear(self: *@This(), key: Pixel) ?struct { Pixel, []ImageCoord } {
            const found, const data = self.tree.getNear(key) orelse return null;
            if (data.items.len == 0) {
                return null;
            }
            return .{ found, data.items };
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
