# t8code_GPU_experiments

A project using t8code build from the t8code_project template (https://github.com/dlr-amr/t8code_project).


## Build

After project setup, build with:

```bash
mkdir build && cd build
cmake ..
make
```

Dependencies (t8code and googletest) are automatically fetched from GitHub using CMake's FetchContent.
No manual installation is required.

We also recommend to use:

```bash
-DCMAKE_INSTALL_PREFIX=/path/to/install
-DCMAKE_INSTALL_LIBDIR=/path/to/install
```

when you want to locally install this package.

To install, simply call:

```bash
make install
```

## Folders and Files

### CMakeLists.txt

The basic CMake setting file.
Note that src/ and test/ folders have their own additional CMake files.

### t8code_GPU_experiments_main.cxx

A main program to show how to use and link against t8code_GPU_experiments and how to set up t8code.

### src/

Contains the source and header files of the library.

src/t8code_GPU_experiments_helloworld.cxx and .hxx define a very simple hello world function.

### src/advect/

A copy of the advect example of t8code (example/advect and example/common), built as the executable
`t8code_GPU_experiments_advect`. It solves the advection equation with adaptive mesh refinement and
is the starting point for the GPU development. Run it with `-h` for a list of options.

### test/

Contains tests.

### cmake/

Contains CMake configuration files for fetching thirdparty dependencies (t8code, googletest) using FetchContent.
