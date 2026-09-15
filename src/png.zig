const std = @import("std");
const colours = @import("colours.zig");
const c_png = @cImport({
    @cInclude("png.h");
});

pub const Data = struct {
    data: []colours.Pixel,
    size_x: u32,
    size_y: u32,
};

pub fn readPng(
    allocator: std.mem.Allocator,
    io: std.Io,
    path: []const u8,
) !Data {
    const png_file = try std.Io.Dir.cwd().openFile(io, path, .{});
    const png_file_fd = c_png.fdopen(png_file.handle, "rb");
    if (png_file_fd == null) {
        return error.fdError;
    }
    defer _ = c_png.fclose(png_file_fd);

    const png_ptr = c_png.png_create_read_struct(
        c_png.PNG_LIBPNG_VER_STRING,
        null,
        dieOnPngError,
        null,
    );
    if (png_ptr == null) {
        return error.LibPngError;
    }

    const info_ptr = c_png.png_create_info_struct(png_ptr);
    if (info_ptr == null) {
        c_png.png_destroy_read_struct(
            @constCast(&png_ptr),
            null,
            null,
        );
    }
    defer c_png.png_destroy_read_struct(
        @constCast(&png_ptr),
        @constCast(&info_ptr),
        null,
    );

    c_png.png_init_io(png_ptr, png_file_fd);
    c_png.png_read_info(png_ptr, info_ptr);

    const color_type = c_png.png_get_color_type(png_ptr, info_ptr);
    const bit_depth = c_png.png_get_bit_depth(png_ptr, info_ptr);
    const size_x = c_png.png_get_image_width(png_ptr, info_ptr);
    const size_y = c_png.png_get_image_height(png_ptr, info_ptr);
    const image_size = size_x * size_y;
    if (color_type == c_png.PNG_COLOR_TYPE_GRAY and bit_depth < 8) {
        c_png.png_set_expand_gray_1_2_4_to_8(png_ptr);
    }
    if (color_type == c_png.PNG_COLOR_TYPE_GRAY or color_type == c_png.PNG_COLOR_TYPE_GRAY_ALPHA) {
        c_png.png_set_gray_to_rgb(png_ptr);
    }
    if (color_type == c_png.PNG_COLOR_TYPE_PALETTE) {
        c_png.png_set_palette_to_rgb(png_ptr);
    }
    if ((color_type & c_png.PNG_COLOR_MASK_ALPHA) != 0) {
        c_png.png_set_strip_alpha(png_ptr);
    }
    if (bit_depth == 16) {
        c_png.png_set_strip_16(png_ptr);
    }
    c_png.png_read_update_info(png_ptr, info_ptr);

    var data = try allocator.alloc(colours.Pixel, image_size);
    errdefer allocator.free(data);
    const pixel_size = 3;
    const buffer = try allocator.alloc(u8, image_size * pixel_size);
    defer allocator.free(buffer);
    for (0..size_y) |y| {
        const row = buffer[pixel_size * size_x * y .. (pixel_size * size_x + 1) * y];
        c_png.png_read_row(png_ptr, row.ptr, null);
    }

    for (0..image_size) |idx| {
        const red = buffer[idx * pixel_size];
        const green = buffer[idx * pixel_size + 1];
        const blue = buffer[idx * pixel_size + 2];
        data[idx] = colours.Pixel{
            .red = red,
            .green = green,
            .blue = blue,
            .alpha = 0xFF,
        };
    }
    return .{
        .data = data,
        .size_x = size_x,
        .size_y = size_y,
    };
}

pub fn writePng(
    allocator: std.mem.Allocator,
    io: std.Io,
    data: []colours.Pixel,
    path: []const u8,
    size_x: u32,
    size_y: u32,
) !void {
    const png_ptr = c_png.png_create_write_struct(
        c_png.PNG_LIBPNG_VER_STRING,
        null,
        dieOnPngError,
        null,
    );
    if (png_ptr == null) {
        return error.LibPngError;
    }
    const info_ptr = c_png.png_create_info_struct(png_ptr);
    if (info_ptr == null) {
        c_png.png_destroy_write_struct(@constCast(&png_ptr), null);
        return error.LibPngError;
    }
    defer c_png.png_destroy_write_struct(@constCast(&png_ptr), @constCast(&info_ptr));

    const depth = 8;
    c_png.png_set_IHDR(
        png_ptr,
        info_ptr,
        size_x,
        size_y,
        depth,
        c_png.PNG_COLOR_TYPE_RGB,
        c_png.PNG_INTERLACE_NONE,
        c_png.PNG_COMPRESSION_TYPE_DEFAULT,
        c_png.PNG_FILTER_TYPE_DEFAULT,
    );

    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();
    const arena_alloc = arena.allocator();
    const rows = try arena_alloc.alloc(c_png.png_bytep, size_y);
    const pixel_size = 3;
    var pixel_idx: usize = 0;
    for (0..size_y) |y| {
        const row = try arena_alloc.alloc(c_png.png_byte, size_x * pixel_size);
        for (0..size_x) |x| {
            const pixel = data[pixel_idx];
            const row_data_idx = x * pixel_size;
            row[row_data_idx] = pixel.red;
            row[row_data_idx + 1] = pixel.green;
            row[row_data_idx + 2] = pixel.blue;
            pixel_idx += 1;
        }
        rows[y] = @ptrCast(@constCast(row.ptr));
    }

    const png_file = try std.Io.Dir.cwd().createFile(io, path, .{});
    const png_file_fd = c_png.fdopen(png_file.handle, "wb");
    if (png_file_fd == null) {
        return error.fdError;
    }
    defer _ = c_png.fclose(png_file_fd);

    c_png.png_init_io(png_ptr, png_file_fd);
    c_png.png_write_info(png_ptr, info_ptr);
    c_png.png_write_image(png_ptr, rows.ptr);
    c_png.png_write_end(png_ptr, null);
}

// callback to abort if libpng has an error
fn dieOnPngError(
    png_ptr: c_png.png_structp,
    msg: [*c]const u8,
) callconv(.c) noreturn {
    _ = png_ptr;
    const err_msg = if (msg != null) std.mem.span(msg) else "Unknown png error";
    std.debug.print("libpng fatal error: {s}\n", .{err_msg});

    // libpng requires this function to never return.
    // Abort to prevent longjmp from firing across Zig stack frames.
    std.process.abort();
}
