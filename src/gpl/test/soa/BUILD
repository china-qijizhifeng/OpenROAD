# SPDX-License-Identifier: BSD-3-Clause
load("@rules_cc//cc:cc_test.bzl", "cc_test")

cc_test(
    name = "storage_test",
    size = "small",
    srcs = ["storage_test.cpp"],
    deps = ["//src/gpl:soa_storage"],
)

cc_test(
    name = "density_test",
    size = "small",
    srcs = ["density_test.cpp"],
    deps = ["//src/gpl:density_soa"],
)

cc_test(
    name = "point_test",
    size = "small",
    srcs = ["point_test.cpp"],
    deps = ["//src/gpl:point_columns"],
)

cc_test(
    name = "metadata_test",
    size = "small",
    srcs = ["metadata_test.cpp"],
    deps = ["//src/gpl:placement_columns"],
)
