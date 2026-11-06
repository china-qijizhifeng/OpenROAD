# OpenROAD - AI agent context

This file provides project-specific guidance for AI agent sessions working on the OpenROAD codebase.

## Quick Reference

Detailed guides are in `docs/agents/` subdirectory:

- **Build Pitfalls**: See `docs/agents/build.md`
- **Testing Guide**: See `docs/agents/testing.md`
- **Git & CI**: See `docs/agents/ci.md`
- **Coding Patterns**: See `docs/agents/coding.md`

## Critical Rules

1. **Ask before modifying `src/sta/` files** -- OpenSTA is managed upstream (`Sdc.cc`, `Power.cc`, `Sdc.tcl`, etc.). Prefer fixes in OpenROAD code (e.g., src/dbSta/, src/rsz/) when possible. 
2. Run `clang-format -i <files>` for C++ files before commit. **NEVER** for `src/sta/*` and `*.i` files.
3. **Always use `git commit -s`** for DCO compliance.  Do not review DCO compliance; CI handles it.
4. When amending submodule commits, parent repo submodule reference must also be updated via `git submodule update --init --recursive`. It is needed after any merge/pull.
5. **Trace bugs upstream** -- when a bug appears in output (e.g., Verilog), find the data creation point (e.g., `buffer_ports`, `remove_buffers`), not the serialization point (e.g., `VerilogWriter`).

## Build and Test

Bazel is the only build system, on Linux and macOS.
Register new source files, tools and tests in Bazel `BUILD` files.

- Build: `bazelisk build //:openroad` (the binary is `bazel-bin/openroad`)
- Test one module: `bazelisk test //src/<module>/test/...`
- Test one Tcl test: `bazelisk test --test_output=errors //src/<module>/test:<name>-tcl_test`

## FlashEDA profiling environment

- Delegate build, test and profiling execution to subagents so the main agent can continue code and result analysis.
- For parallel profiling, prefer E2B through `hk-linlizhi-jumpserver` using the `SWEEnvMaker 0820` environment. Read the `e2b-code-validation` skill and verify the actual project, template, SDK and resources on HK before execution; do not infer a template ID from this name.
- The verified HK launcher is `/root/SWEEnvMaker/bins/e2b/0820_run_nac_swe_q1_q45_pass4_c40_chenjun_template_id_canary.sh`, calling `/root/code-harbor/scripts/_harbor_entrypoint.py run`. Reuse its environment contract; do not launch its existing SWE benchmark batch for profiling.
- On 2026-11-06 the 0820 template API exposed the ready FlashEDA template `flasheda-r5-cpu8-t16-20261030` (`vshaqeiaoi96u3cp66t3`, 8 vCPU / 16 GiB). Recheck availability and actual sandbox resources before use. Credentials belong only in client options and must not be logged or uploaded.
- Split complete designs across independent sandboxes. Run the native/candidate pair serially inside each sandbox, keep profiling separate from uninstrumented timing, and record source, input, binary, environment and sandbox identities. Do not pool absolute timings from different hosts.

## AI Agent Skills

Skills are located in `.agents/skills/` (with `.claude/skills` symlink for Claude Code).

| Skill | Purpose | Invocation |
|-------|---------|------------|
| `triage-issue` | Reproduce bug and minimize test case with whittle.py | `/triage-issue <issue#>` |
| `fix-bug` | Trace root cause, implement fix, create tests, prepare commit | `/fix-bug <issue#-or-error-code>` |
| `add-test` | Add integration/unit tests registered in Bazel | `/add-test <module> [description]` |
| `review-pr` | Draft local PR review notes (correctness > QoR > testing); human posts | `/review-pr <pr#-or-url>` |
