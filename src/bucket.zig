const std = @import("std");

const LEAF_SIZE = 128;

/// A simple collection of keys and values backed by static arrays.
pub fn Bucket(comptime K: type, comptime V: type, comptime distance_func: fn (K, K) usize) type {
    return struct {
        keys: [LEAF_SIZE]K = std.mem.zeroes([LEAF_SIZE]K),
        values: [LEAF_SIZE]V = std.mem.zeroes([LEAF_SIZE]V),
        count: u8 = 0,

        /// Sets the value for the given key within this bucket.
        /// Returns an error if the bucket is full.
        pub fn putValue(self: *@This(), key: K, value: V) !void {
            if (self.count == self.keys.len) {
                return error.BucketFull;
            }
            const idx = self.findForInsertion(key);
            self.keys[idx] = key;
            self.values[idx] = value;
            if (idx == self.count) {
                self.count += 1;
            }
        }

        /// Force puts a new key value without checking for duplicates.
        /// Only use when certain there are no duplicates.
        pub fn simplePut(self: *@This(), key: K, value: V) !void {
            if (self.count == self.keys.len) {
                return error.BucketFull;
            }
            const idx = self.count;
            self.keys[idx] = key;
            self.values[idx] = value;
            self.count += 1;
        }

        /// Returns the value associate with this key.
        pub fn getValue(self: *@This(), key: K) ?V {
            for (0..self.count) |idx| {
                if (std.meta.eql(key, self.keys[idx])) {
                    return self.values[idx];
                }
            }
            return null;
        }

        /// Returns a pointer to the value in the store
        pub fn getValuePtr(self: *@This(), key: K) ?*V {
            for (0..self.count) |idx| {
                if (std.meta.eql(key, self.keys[idx])) {
                    return &self.values[idx];
                }
            }
            return null;
        }

        /// Returns the key closest to the given key in this leaf.
        pub fn getNearest(self: *@This(), key: K) ?struct { K, V } {
            var nearest: ?struct { K, V } = null;
            var d_min: usize = 0xFFFFFFFF;
            for (0..self.count) |idx| {
                const candidate = self.keys[idx];
                const distance = distance_func(key, candidate);
                if (distance < d_min) {
                    d_min = distance;
                    nearest = .{ candidate, self.values[idx] };
                }
            }
            return nearest;
        }

        /// Returns the index into this leaf for the key to insert.
        /// If the key already exists it returns the occupied slot.
        fn findForInsertion(self: *@This(), key: K) usize {
            for (0..self.count) |idx| {
                if (std.meta.eql(self.keys[idx], key)) {
                    return idx;
                }
            }
            return self.count;
        }

        /// Removes the first instance of the passed in value
        pub fn removeKey(self: *@This(), key: K) !void {
            for (0..self.count) |idx| {
                if (std.meta.eql(key, self.keys[idx])) {
                    @memmove(self.keys[idx .. self.keys.len - 1], self.keys[idx + 1 ..]);
                    @memmove(self.values[idx .. self.values.len - 1], self.values[idx + 1 ..]);
                    self.count -= 1;
                    return;
                }
            } else {
                return error.BucketKeyNotFound;
            }
        }
    };
}

fn intDistance(left: i32, right: i32) usize {
    const v: i32 = right - left;
    return @intCast(v * v);
}

test "Bucket add" {
    var bucket: Bucket(i32, i32, intDistance) = .{};

    try bucket.putValue(17, 1024);
    try bucket.putValue(19, 96);

    try std.testing.expectEqual(2, bucket.count);
    try std.testing.expectEqual(1024, bucket.getValue(17));
    try std.testing.expectEqual(96, bucket.getValue(19));
}

test "Bucket remove" {
    var bucket = Bucket(i32, i32, intDistance){};

    try bucket.putValue(12, 1);
    try bucket.putValue(19, 2);
    try bucket.putValue(7, 3);
    try bucket.removeKey(19);

    try std.testing.expectEqual(2, bucket.count);
    try std.testing.expectEqual(1, bucket.getValue(12));
    try std.testing.expectEqual(null, bucket.getValue(19));
    try std.testing.expectEqual(3, bucket.getValue(7));
}

test "Bucket nearest" {
    var bucket = Bucket(i32, i32, intDistance){};

    try bucket.putValue(10, 1);
    try bucket.putValue(20, 2);
    try bucket.putValue(30, 3);
    try bucket.putValue(40, 4);

    try std.testing.expectEqual(
        .{ 20, 2 },
        bucket.getNearest(16),
    );
    try std.testing.expectEqual(
        .{ 30, 3 },
        bucket.getNearest(32),
    );
    try std.testing.expectEqual(
        .{ 10, 1 },
        bucket.getNearest(-1000),
    );
}

test "Bucket duplicate" {
    var bucket = Bucket(i32, i32, intDistance){};
    const value: i32 = 3142;

    try bucket.putValue(value, 16);
    try bucket.putValue(value, 1);
    try bucket.putValue(value, 17);

    try std.testing.expectEqual(1, bucket.count);
    try std.testing.expectEqual(
        .{ value, 17 },
        bucket.getNearest(0),
    );
}

test "Bucket add returns error on full" {
    var bucket = Bucket(i32, i32, intDistance){};
    for (0..LEAF_SIZE) |idx| {
        const k: i32 = @intCast(idx);
        try bucket.putValue(k, k);
    }

    try std.testing.expectError(
        error.BucketFull,
        bucket.putValue(-1, -1),
    );
}
