const Pixel = @import("colours.zig").Pixel;

pub const ImageData = struct {
    buffer: []Pixel,
    size_x: usize,
    size_y: usize,

    pub fn at(self: *const ImageData, coord: ImageCoord) Pixel {
        const idx = self.toIndex(coord);
        return self.buffer[idx];
    }

    pub fn put(self: *ImageData, coord: ImageCoord, value: Pixel) void {
        const idx = self.toIndex(coord);
        self.buffer[idx] = value;
    }

    pub fn toIndex(self: *const ImageData, coord: ImageCoord) usize {
        const bounded = self.boundCoord(coord);
        return self.size_x * bounded.y + bounded.x;
    }

    pub fn fromIndex(self: *const ImageData, idx: usize) ImageCoord {
        const bounded = idx % (self.size_x * self.size_y);
        const y = bounded / self.size_x;
        const x = bounded % self.size_x;
        return .{
            .x = x,
            .y = y,
        };
    }

    // Ensures the coordinate is within the bounds of the image data.
    // Wraps with the image bounds.
    pub fn boundCoord(self: *const ImageData, coord: ImageCoord) ImageCoord {
        return .{
            .x = coord.x % self.size_x,
            .y = coord.y % self.size_y,
        };
    }
};

pub const ImageCoord = struct {
    x: usize,
    y: usize,
};
