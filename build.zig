const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimization = b.standardOptimizeOption(.{});
    const main = b.addExecutable(.{
        .name = "main",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = b.graph.host,
            .optimize = optimization,
            .link_libc = true,
        }),
    });
    main.root_module.linkSystemLibrary("libpng", .{});
    const main_check = b.addExecutable(.{
        .name = "main",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = b.graph.host,
            .optimize = optimization,
            .link_libc = true,
        }),
    });
    main_check.root_module.linkSystemLibrary("libpng", .{});

    const check = b.step("check", "Check compilation");
    check.dependOn(&main_check.step);

    b.installArtifact(main);
}
