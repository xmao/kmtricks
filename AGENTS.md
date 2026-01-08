# Repository Guidelines

## Project Structure & Module Organization
- `src/`: core C++ implementation of kmtricks.
- `include/`: public and internal headers.
- `tests/`: GoogleTest suites and CTest integration.
- `plugins/`: optional plugin sources and build targets.
- `cmake/`: CMake helper modules; top-level build config in `CMakeLists.txt`.
- `scripts/`, `doc/`, `docker/`, `conda/`: release tooling, docs, and packaging helpers.
- `thirdparty/`: vendored dependencies (do not edit unless required).

## Build, Test, and Development Commands
- `./install.sh -r Release -t 2 -j 8`: recommended build script; configures, builds, and runs tests (see `./install.sh -h` for flags).
- `mkdir build && cd build && cmake .. -DCOMPILE_TESTS=ON && make -j8`: manual build from source.
- `ctest --verbose`: run tests after a build that enabled `COMPILE_TESTS`.
- `make plugins`: build plugin targets when plugin support is enabled (use `-p` in `install.sh`).

## Coding Style & Naming Conventions
- C++17 is required (`CMAKE_CXX_STANDARD 17`).
- Follow existing formatting: 2-space indentation, braces on their own line, and compact include grouping.
- Naming follows mixed conventions (e.g., `KmDir`, `kmtricksCli`); mirror nearby code instead of introducing new patterns.
- No auto-formatter is enforced; keep diffs minimal and localized.

## Testing Guidelines
- Tests use GoogleTest; files live under `tests/` and end with `_test.cpp`.
- Build test binaries via `-DCOMPILE_TESTS=ON`; run with `ctest --verbose` or directly from `tests/`.
- Tests create temporary data under `tests_tmp/`; clean up if you add new fixtures.

## Commit & Pull Request Guidelines
- Commit history suggests short, direct subjects (e.g., "prepare 1.5.1", "fix merge for hash:bf:bin"). Keep messages concise and descriptive.
- PRs should include: purpose/impact summary, build/test commands run, and any new flags or data requirements.
- Link relevant issues or wiki references when changing behavior or CLI options.

## Notes & Operational Tips
- kmtricks can be disk-intensive; document expected space usage for new workflows.
- If you update CLI behavior, consider adding or adjusting tests in `tests/` to cover regressions.
