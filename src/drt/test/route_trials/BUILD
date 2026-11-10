# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

load("@rules_python//python:defs.bzl", "py_library", "py_test")

py_library(
    name = "route_trials",
    srcs = [
        "feedback.py",
        "joint_adapter.py",
        "select_feedback.py",
    ],
    imports = ["."],
)

[
    py_test(
        name = name,
        srcs = [name + ".py"],
        deps = [":route_trials"],
    )
    for name in [
        "test_feedback",
        "test_joint_adapter",
        "test_select_feedback",
    ]
]

exports_files(["import_joint.tcl"])
