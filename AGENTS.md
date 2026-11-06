# AGENTS.md

Guidance for coding agents working in this repository.

## What this is

`t8code_GPU_experiments` is a small experimental C++ project that links against
[t8code](https://github.com/dlr-amr/t8code) (adaptive meshes / forests of space-trees, MPI-parallel)
and adds **OpenMP GPU offloading** (AMD `gfx*` or NVIDIA `sm_*`). It was generated from the
[t8code_project](https://github.com/dlr-amr/t8code_project) template. 

So far it holds template scaffolding, a copy of t8code tutorial **step5** (element data plus ghost exchange),
and a copy of t8code's **advect example** (`src/advect/`). In the advect example, the element loop of the solver
runs on the GPU with OpenMP `target` and unified shared memory (see "GPU offloading in the advect example" below).
It is the only GPU kernel so far.

## Layout

| Path | Purpose |
|---|---|
| `CMakeLists.txt` | Top-level build. Pulls in thirdparty and OpenMP, builds library (`src/`), three executables, tests. |
| `cmake/thirdparty.json` | Dependency list read by `thirdparty.cmake`. Currently only t8code, pinned to a commit hash. |
| `cmake/thirdparty.cmake` | Parses the JSON and runs `FetchContent_Declare/MakeAvailable` for each `git` entry. `depends_on_cmake_option` skips an entry when that option is ON. |
| `cmake/OpenMPOffload.cmake` | GPU detection and the helper `t8_GPU_enable_openmp(<target>)`. |
| `src/CMakeLists.txt` | Library `t8code_GPU_experiments`, VTK handling, install/package config. |
| `src/t8code_GPU_experiments_helloworld.{hxx,cxx}` | Template hello world in namespace `t8code_GPU_experiments`. |
| `src/t8code_GPU_experiments_data.{hxx,cxx}` | t8code step5 example adapted for this project (`t8code_GPU_experiments_main`). |
| `src/t8code_GPU_experiments_data_main.cxx` | `main()` wrapper for the data example. |
| `src/advect/t8code_GPU_experiments_advect.cxx` | Based on t8code `example/advect/t8_advection.cxx` (finite-volume advection solver with adaptive refinement, contains `main`). Element data as struct of arrays, element loop offloaded to the GPU. |
| `src/advect/t8code_GPU_experiments_example_common*.{hxx,cxx}` | Copy of t8code `example/common/` (flow fields, level-set functions, adapt callbacks). |
| `t8code_GPU_experiments_main.cxx` | `main()` for the hello world executable (MPI/sc/t8 init, then finalize). |
| `test/` | GoogleTest tests plus shared `gtest_main` that initializes MPI/sc/t8code. |

### Build targets

- `t8code_GPU_experiments`: library (helloworld + data sources), links `T8CODE::T8` privately.
- `t8code_GPU_experiments_helloworld`: built from `t8code_GPU_experiments_main.cxx`.
- `t8code_GPU_experiments_data`: built from `src/t8code_GPU_experiments_data_main.cxx`. When run, it writes
  `t8code_GPU_experiments_forest*.vtu` and `t8code_GPU_experiments_forest_with_volume_data*.vtu`
  into the working directory.
- `t8code_GPU_experiments_advect`: built from the files in `src/advect/`. It has the same command line as t8code's
  `t8_advection` (`-h` for help). Example: `mpirun -np 2 ./t8code_GPU_experiments_advect -u 1 -e 4 -l 3 -r 1 -T 0.1 -o`.
- Tests: `t8code_GPU_experiments_gtest_example`, `t8code_GPU_experiments_gtest_vtk_linkage`.

## Building

```bash
mkdir build && cd build
cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang   # Clang-based compiler needed for offloading (e.g. ROCm amdclang++)
make -j
ctest            # each test runs as `mpirun -np 4 <exe>`
```

- On the first configure, t8code is downloaded and built through FetchContent. This takes a while and needs network access.
- **`OFFLOAD_ARCH`** cache variable (see `cmake/OpenMPOffload.cmake`):
  - empty (default): autodetect with `amdgpu-arch` / `nvptx-arch`. These are searched in `$ROCM_PATH`, `/opt/rocm`, and `PATH`. `gfx000` (the ROCm CPU agent) is filtered out.
    Target features such as `:xnack-` are stripped, so the image runs with and without `HSA_XNACK=1`.
  - `gfx942`, `gfx90a`, `sm_80`, ... (comma or `;` list): set it explicitly. You need this on login nodes that have no GPU.
  - `NONE`: build host-only OpenMP and silence the "no GPU" warning.
  - If the compiler is not Clang, offloading is turned off with a warning.
- MPI is not found explicitly (the `find_package(MPI)` lines are commented out). MPI support comes through t8code/sc.
  `ctest` still calls `mpirun`, so `mpirun` must be on `PATH`.
- VTK: if t8code was built with `T8CODE_ENABLE_VTK`, this project turns on `T8CODE_GPU_EXPERIMENTS_ENABLE_VTK`
  automatically and checks that both use the same VTK version.

## Conventions

- **Every new compiled target must call `t8_GPU_enable_openmp(<target>)`.** It links `OpenMP::OpenMP_CXX` and adds
  `--offload-arch` to both the compile and link options. Without it, offloaded code silently runs on the host or fails to link.
- New library sources go into `add_library(...)` in `src/CMakeLists.txt`. Public headers must also be added to
  the `install(FILES ...)` list there.
- New tests: add `add_t8code_GPU_experiments_cpp_test(NAME <name> SOURCES <file>.cxx)` in `test/CMakeLists.txt`.
  The shared `t8code_GPU_experiments_gtest_main.cxx` is added automatically, so tests must not define `main`.
- Naming: files and symbols use the prefix `t8code_GPU_experiments_`. C++ code goes in `namespace t8code_GPU_experiments`.
  Macros and CMake variables use `T8CODE_GPU_EXPERIMENTS_`. CMake helpers use `t8_GPU_`.
- Style follows the [t8code coding guidelines](https://github.com/DLR-AMR/t8code/wiki/Coding-Guideline):
  GNU-ish formatting with a space before `(` in calls, `/* */` comments, and doxygen comments on declarations.
  Every program follows the order `sc_MPI_Init` → `sc_init` → `t8_init` … `sc_finalize` → `sc_MPI_Finalize`.
- To bump t8code, change `ref` in `cmake/thirdparty.json`.
- `src/advect/` was copied from t8code at commit `cda4653` (noted in each file). Everything the copy defines was renamed
  from `t8_*` to `t8code_GPU_experiments_*` (prefixes `advect_`, `flow_`, `common_`, `scalar3d_`, `scalar_function_`,
  `example_`, `levelset_`, `incomp_cube_`). To diff against upstream, map these prefixes back first.
  Calls into the t8code API (`t8_forest_*`, `t8_cmesh_*`, `t8_vec*`, ...) keep their names.
- PR titles start with `Documentation:`, `Bugfix:`, `Feature:`, `Improvement:` or `Other:` (see `pull_request_template.md`).
- Commits, issues and PRs created by an AI agent must contain `[created with claude]`.

## GPU offloading in the advect example

- The element loop in `t8code_GPU_experiments_advect_solve` has `#pragma omp target teams distribute parallel for`.
  The file has `#pragma omp requires unified_shared_memory`: the kernel reads host memory directly, with no `map` clauses.
- Element data (`t8code_GPU_experiments_advect_element_data_t`) is a struct of arrays. Per-face fields are
  per-element blocks, e.g. `data->fluxes[ielem][iface][ineigh]`, with at most `MAX_FACE_NEIGHBORS` (4) neighbors per face.
  `phi_values` (`sc_array_t`, one double per local element and ghost) is used through a raw `double *phi`.
- **Code inside the loop must not call t8code or sc functions, use `problem->`, or call function pointers.**
  - Anything that needs t8code (neighbors, face center/normal/area, including child faces of hanging faces) is precomputed on the host in
    `compute_mesh_data` after each adapt or partition (and in `init_elements`).
  - Functions called in the loop are between `#pragma omp declare target` and `#pragma omp end declare target`.
- **Gather, not scatter**: each element writes only its own entries, so each inner face flux is computed from both sides.
  Writing into neighbor data from the loop would cause race conditions.
- The flow field is selected by number (`problem->flow`, `flow_eval` with a `switch`), because `problem->u` is a host
  function pointer. A new flow function must be added to `choose_flow` and to `flow_eval`, and it must be inside the
  `declare target` region in `t8code_GPU_experiments_example_common{.hxx,_functions.cxx}`.
- Running on a GPU needs unified memory support: on AMD gfx90a/gfx942 set `HSA_XNACK=1` (possibly `OFFLOAD_ARCH=gfx90a:xnack+`);
  MI300A supports it natively. On NVIDIA it needs HMM/ATS. Use `OMP_TARGET_OFFLOAD=MANDATORY` to fail instead of falling back
  to the host, and `LIBOMPTARGET_INFO=4` to see the kernel launches. Without offloading (e.g. GCC or `OFFLOAD_ARCH=NONE`) the loop runs on the host.
- Verification: the results must match the original t8code `t8_advection` (`_deps/t8code-src/example/advect/`)
  for quad, hex, tet, prism, hybrid and 1D meshes with several MPI ranks (start/end volume, l_infty, L_2,
  element counts per step, `ghost_sent`, `partition_procs_sent`). The original crashes with `-s`; ours does not.

## Known issues and gotchas

- `src/config.cmake.in` broke during template substitution. `@PACKAGE_INIT@` is missing,
  the file includes `-targets.cmake`, and it calls `check_required_components( )` with no argument. Also, no
  `install(EXPORT ...)` exists. As a result, `find_package(T8CODE_GPU_EXPERIMENTS)` from an installed copy will not work.
  `make install` copies the files, but the package config cannot be used.
- The tests link `GTest::gtest_main`, but `thirdparty.json` does not fetch googletest (even though README.md says it does).
  The target has to come from t8code's own build. If it is missing, add googletest to `thirdparty.json`.
- The data example builds a **uniform** forest (level 3) and does not adapt it. The step5 adapt callback was not carried over.
- License headers are inconsistent: GPLv3 in template files, GPLv2+ in files copied from t8code. Contributors need a
  BSD statement in `doc/author_NAME.txt` (see CONTRIBUTING.md). The `doc/` folder does not exist yet.
