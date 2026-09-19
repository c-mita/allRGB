const std = @import("std");

pub fn ArgParser(
    comptime flags: anytype,
) type {
    const type_info = @typeInfo(@TypeOf(flags));
    if (type_info != .@"struct" and type_info.@"struct".is_tuple) {
        @compileError("Expected a tuple literal of flags");
    }
    const flag_values = type_info.@"struct".fields;
    var flag_types: [flag_values.len]type = undefined;
    var flag_names: [flag_values.len][]const u8 = undefined;
    var flag_attrs: [flag_values.len]std.builtin.Type.StructField.Attributes = undefined;

    inline for (flags, 0..) |flag, idx| {
        flag_types[idx] = flag.flag_type;
        flag_names[idx] = std.fmt.comptimePrint("{s}", .{flag.flag_name});
        const def = struct {
            comptime {
                _ = flag.flag_default;
            }
            const value: flag.flag_type = @as(flag.flag_type, flag.flag_default);
        };
        flag_attrs[idx] = .{
            .default_value_ptr = @ptrCast(&def.value),
        };
    }

    const parameters_struct: type = @Struct(
        .auto,
        null,
        &flag_names,
        &flag_types,
        &flag_attrs,
    );

    return struct {
        pub fn parse(args: std.process.Args) !parameters_struct {
            var params: parameters_struct = .{};
            var it = args.iterate();
            // Skip executable name
            _ = it.next();
            while (it.next()) |arg| {
                inline for (flags) |flag| {
                    if (flag.matches(arg)) {
                        @field(params, flag.flag_name) = try flag.parse(arg, &it);
                        break;
                    }
                } else {
                    return error.FlagParseError;
                }
            }
            return params;
        }
    };
}

pub fn ValueFlag(
    comptime T: type,
    comptime name: []const u8,
    comptime short_name: ?[]const u8,
    comptime default: T,
) type {
    const full_flag = "--" ++ name;
    const short_flag = if (short_name != null) "-" ++ short_name.? else null;
    const type_info = @typeInfo(T);

    return struct {
        pub const flag_name: []const u8 = name;
        pub const flag_short_name: ?[]const u8 = short_name;
        pub const flag_default: T = default;
        pub const flag_type: type = T;

        pub fn matches(arg: [:0]const u8) bool {
            return (std.mem.eql(
                u8,
                short_flag,
                arg,
            ) or std.mem.eql(
                u8,
                full_flag,
                arg,
            ));
        }

        pub fn parse(
            _: [:0]const u8,
            iterator: *std.process.Args.Iterator,
        ) !T {
            const value = iterator.next() orelse return error.FlagParseError;
            const t_info = if (type_info == .optional) @typeInfo(type_info.optional.child) else type_info;
            return switch (t_info) {
                .array => {
                    if (t_info.array.child == u8) {
                        return value;
                    } else {
                        @compileError("Unsupported type " ++ @typeName(T));
                    }
                },
                .pointer => {
                    if (t_info.pointer.child == u8) {
                        return value;
                    } else {
                        @compileError("Unsupported type " ++ @typeName(T));
                    }
                },
                .int => std.fmt.parseInt(T, value, 10),
                .float => std.fmt.parseFloat(T, value, 10),
                else => @compileError("Unsupported type " ++ @typeName(T)),
            };
        }
    };
}

pub fn BooleanFlag(
    comptime name: []const u8,
    comptime default: bool,
) type {
    const positive_flag = "--" ++ name;
    const negative_flag = "--no" ++ name;

    return struct {
        pub const flag_name: []const u8 = name;
        pub const flag_default: bool = default;
        pub const flag_type: type = bool;

        pub fn matches(arg: [:0]const u8) bool {
            return std.mem.eql(
                u8,
                positive_flag,
                arg,
            ) or std.mem.eql(
                u8,
                negative_flag,
                arg,
            );
        }

        pub fn parse(
            arg: [:0]const u8,
            _: *std.process.Args.Iterator,
        ) !bool {
            if (std.mem.eql(u8, positive_flag, arg)) {
                return true;
            }
            return false;
        }
    };
}

pub fn EnumFlag(
    comptime values: type,
    comptime name: []const u8,
    comptime default: values,
) type {
    const enum_fields = std.enums.values(values);

    return struct {
        pub const flag_name: []const u8 = name;
        pub const flag_default: values = default;
        pub const flag_type: type = values;

        pub fn matches(arg: [:0]const u8) bool {
            inline for (enum_fields) |val| {
                const tag = "--" ++ @tagName(val);
                if (std.mem.eql(u8, tag, arg)) {
                    return true;
                }
            }
            return false;
        }

        pub fn parse(
            arg: [:0]const u8,
            _: *std.process.Args.Iterator,
        ) !values {
            inline for (enum_fields) |val| {
                const tag = "--" ++ @tagName(val);
                if (std.mem.eql(u8, tag, arg)) {
                    return val;
                }
            } else {
                return error.FlagParseError;
            }
        }
    };
}
