const std = @import("std");
const colours = @import("colours.zig");
const c_png = @cImport({
    @cInclude("png.h");
});

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
