# IMPF

A C library for numerical optimization and approximations.



[![Build & Test](https://github.com/zhuanglinsheng/impf/actions/workflows/cmake.yml/badge.svg)](https://github.com/zhuanglinsheng/impf/actions/workflows/cmake.yml) 

[Documentations](https://zhuanglinsheng.github.io/impf/) 



## Features

- The whole package is writen in ISO C90 (ANSI C89). 
- The interfaces follow the classical Fortran subroutine style. 
- No dynamically memory allocation is used. 
- Currently, this package only deal with `double` precision float numbers. 



## Usage

### Dependencies

Building IMPF requires the following to be installed:

- [CMake](https://cmake.org/) 
- A C compiler, e.g. [GCC](https://gcc.gnu.org/) or [Clang](https://clang.llvm.org/) 
- [BLAS](https://netlib.org/blas/) and [LAPACK](https://netlib.org/lapack/), eg [OpenBLAS](https://www.openblas.net/) 

### Installation

Download the source code from Github repository

```shell
git clone https://github.com/zhuanglinsheng/impf.git
cd impf
```

Compile the CMake project and install the shared library 

```shell
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install
```

For Linux, the default installation locations of the dynamic sharing library and headers are `/usr/local/lib` and `/usr/local/include`. Users can choose preferred locations in their C/C++ compiler's searching path by editing the `CMakeLists.txt` file.

When users call methods from IMPF, computer will load the dynamic linked library `impf` from the place where the library  is installed (such thas `/usr/local/include`). Make sure your operation system recognizes this place. One easy way is to changing the environment variable `LD_LIBRARY_PATH` (for Linux) in the `~/.bashrc` file by 

```shell
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"
```

### Using IMPF routines 

The file [example/dts/eg_dts_impf_df_norms.c](example/dts/eg_dts_impf_df_norms.c) calculate the norms of an array with IMPF. To compile the code, we execute the command below 

```shell
cc eg_dts_impf_df_norms.c -limpf -lm -lblas -o eg_dts_impf_df_norms
```

Generally, users need to link to the shared libraries `impf`, `m` and `blas`. 



## Contact

Please just contact me if there were any bugs. Many thanks for any feedbacks. 

Pull requests are also welcomed, and please discuss with me first for the pull requests. 

