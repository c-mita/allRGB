const std = @import("std");

pub const Pixel = struct {
    red: u8 = 0,
    green: u8 = 0,
    blue: u8 = 0,
    alpha: u8 = 0,

    /// The square of the euclidean distance of the RGB values
    pub fn distanceSquared(self: *const Pixel, other: Pixel) i32 {
        const rdiff = @as(i32, self.red) - other.red;
        const gdiff = @as(i32, self.green) - other.green;
        const bdiff = @as(i32, self.blue) - other.blue;
        return rdiff * rdiff + gdiff * gdiff + bdiff * bdiff;
    }

    pub fn distance(lhs: Pixel, rhs: Pixel) usize {
        return @intCast(lhs.distanceSquared(rhs));
    }

    pub fn l1Distance(lhs: Pixel, rhs: Pixel) usize {
        const rdiff = @as(i32, @max(lhs.red, rhs.red)) - @min(lhs.red, rhs.red);
        const gdiff = @as(i32, @max(lhs.green, rhs.green)) - @min(lhs.green, rhs.green);
        const bdiff = @as(i32, @max(lhs.blue, rhs.blue)) - @min(lhs.blue, rhs.blue);
        return @intCast(rdiff + gdiff + bdiff);
    }

    pub fn l2Distance(lhs: Pixel, rhs: Pixel) f64 {
        const rdiff = @as(i32, @max(lhs.red, rhs.red)) - @min(lhs.red, rhs.red);
        const gdiff = @as(i32, @max(lhs.green, rhs.green)) - @min(lhs.green, rhs.green);
        const bdiff = @as(i32, @max(lhs.blue, rhs.blue)) - @min(lhs.blue, rhs.blue);
        const dsq = rdiff * rdiff + gdiff * gdiff + bdiff * bdiff;
        return @sqrt(@floatFromInt(dsq));
    }

    pub fn lInfDistance(lhs: Pixel, rhs: Pixel) i32 {
        const rdiff = @as(i32, @max(lhs.red, rhs.red)) - @min(lhs.red, rhs.red);
        const gdiff = @as(i32, @max(lhs.green, rhs.green)) - @min(lhs.green, rhs.green);
        const bdiff = @as(i32, @max(lhs.blue, rhs.blue)) - @min(lhs.blue, rhs.blue);
        return @max(rdiff, gdiff, bdiff);
    }

    /// The signed distance between this point and the splitting plane.
    /// Result is simply self.[red|green|blue] - point
    pub fn distanceToSplit(self: *const Pixel, axis: u32, point: u8) i32 {
        const value = self.selectAxis(axis);
        return @as(i32, value) - point;
    }

    pub fn selectAxis(self: *const Pixel, axis: u32) u8 {
        return switch (axis % 3) {
            0 => self.red,
            1 => self.green,
            2 => self.blue,
            else => unreachable,
        };
    }

    pub fn lessThan(self: *const Pixel, axis: u32, value: u8) bool {
        const self_value = self.selectAxis(axis);
        return self_value < value;
    }

    pub fn splitElements(pixels: []const Pixel) struct { u32, u8 } {
        // find the best splitting point
        var red: usize = 0;
        var green: usize = 0;
        var blue: usize = 0;
        for (pixels) |pixel| {
            red += pixel.red;
            green += pixel.green;
            blue += pixel.blue;
        }
        const red_extra: u8 = if (red % pixels.len != 0) 1 else 0;
        const green_extra: u8 = if (green % pixels.len != 0) 1 else 0;
        const blue_extra: u8 = if (blue % pixels.len != 0) 1 else 0;
        red = (red / pixels.len) + red_extra;
        green = (green / pixels.len) + green_extra;
        blue = (blue / pixels.len) + blue_extra;

        var red_min: u8 = 0xFF;
        var green_min: u8 = 0xFF;
        var blue_min: u8 = 0xFF;
        var red_max: u8 = 0;
        var green_max: u8 = 0;
        var blue_max: u8 = 0;
        for (pixels) |pixel| {
            red_min = @min(red_min, pixel.red);
            green_min = @min(green_min, pixel.green);
            blue_min = @min(blue_min, pixel.blue);
            red_max = @max(red_max, pixel.red);
            green_max = @max(green_max, pixel.green);
            blue_max = @max(blue_max, pixel.blue);
        }

        const red_diff = red_max - red_min;
        const green_diff = green_max - green_min;
        const blue_diff = blue_max - blue_min;

        var split_axis: u32 = 0;
        var split_value: u8 = @intCast(red);
        split_axis = if (red_diff > green_diff) 0 else 1;
        const tmp = if (red_diff > green_diff) red_diff else green_diff;
        split_axis = if (tmp >= blue_diff) split_axis else 2;
        split_value = switch (split_axis) {
            0 => @intCast(red),
            1 => @intCast(green),
            2 => @intCast(blue),
            else => unreachable,
        };
        return .{ split_axis, split_value };
    }
};

pub fn hspCompare(_: void, lhs: Pixel, rhs: Pixel) bool {
    const lh_red = @as(f64, @floatFromInt(lhs.red));
    const lh_green = @as(f64, @floatFromInt(lhs.green));
    const lh_blue = @as(f64, @floatFromInt(lhs.blue));
    const rh_red = @as(f64, @floatFromInt(rhs.red));
    const rh_green = @as(f64, @floatFromInt(rhs.green));
    const rh_blue = @as(f64, @floatFromInt(rhs.blue));
    const bx = 0.299 * lh_red * lh_red + 0.587 * lh_green * lh_green + 0.144 * lh_blue * lh_blue;
    const by = 0.299 * rh_red * rh_red + 0.587 * rh_green * rh_green + 0.144 * rh_blue * rh_blue;

    return if (bx < by) true else false;
}

const HueAngle = struct {
    // hue = atan2(sqrt(3) * (G-B), 2 * R - G - B)
    // Mathematically, atan(sqrt(3) * (G-B) / (2 * R - G - B))
    // Since we only need to "compare" hue values, we can avoid normalization
    // and the computation of atan (or atan2) since atan is strictly
    // increasing within a given quadrant.

    quadrant: i32,
    numerator: i32,
    denominator: i32,

    pub fn of(pixel: Pixel) HueAngle {
        const numerator: i32 = @as(i32, pixel.green) - pixel.blue;
        var denominator: i32 = 2 * @as(i32, pixel.red) - pixel.green - pixel.blue;
        if (numerator == 0 and denominator == 0) {
            denominator = 1;
        }

        const pos_num = numerator >= 0;
        const pos_den = denominator >= 0;
        // odd quadrant ordering but this preserves the ordering of older code
        const quadrant: i32 = if (pos_num and pos_den)
            2
        else if (pos_num and !pos_den)
            3
        else if (!pos_num and !pos_den)
            0
        else
            1;

        return .{
            .quadrant = quadrant,
            .numerator = numerator,
            .denominator = denominator,
        };
    }

    pub fn lessThan(self: *const HueAngle, rhs: HueAngle) bool {
        // Can quickly check quadrants first.
        const lhs = self.*;
        if (lhs.quadrant < rhs.quadrant) {
            return true;
        } else if (rhs.quadrant < lhs.quadrant) {
            return false;
        }

        // Within a quadrant so the signs of the numerator and denominator line up.
        // In this case (n1 / d1 < n2 / d2 ==> n1 * d2 < n2 * d1).
        // And atan is monotonic within a quadrant so we can just perform that check.
        const left = lhs.numerator * rhs.denominator;
        const right = rhs.numerator * lhs.denominator;
        return left < right;
    }
};

/// Returns true if lhs < rhs according to the hue angle
pub fn hueCompare(_: void, lhs: Pixel, rhs: Pixel) bool {
    const left = HueAngle.of(lhs);
    const right = HueAngle.of(rhs);
    return left.lessThan(right);
}

pub fn zOrderCompare(_: void, lhs: Pixel, rhs: Pixel) bool {
    // The Morton Z-Order Curve is the result of interleaving the bits of the
    // RGB values.
    // So [r1..r8][g1..g8][b1..b8] becomes
    // [r1g1b1, r2g2b2, ... r8g8b8]
    // But we're just doing a comparison so we only need to find the most
    // significant bit of the difference between the values for each channel.
    const lhs_data = &[_]u8{ lhs.red, lhs.green, lhs.blue };
    const rhs_data = &[_]u8{ rhs.red, rhs.green, rhs.blue };
    var channel: usize = 0;
    for (1..lhs_data.len) |channel_idx| {
        const lhs_v = lhs_data[channel_idx];
        const rhs_v = rhs_data[channel_idx];
        const current_diff = lhs_v ^ rhs_v;
        const best_diff = lhs_data[channel] ^ rhs_data[channel];
        if ((8 - @clz(current_diff)) > (8 - @clz(best_diff))) {
            channel = channel_idx;
        }
    }
    return lhs_data[channel] < rhs_data[channel];
}
