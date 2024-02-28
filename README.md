# IMPF

Numerical algorithm demos in pure C

- Numerical algorithms implemented in ISO C90 with complete testings
- APIs are designed in the explicit way (no `typedef`)
- No third-party reliance

This repository is just for fun. We don't do any unnecessary optimations to speed up the routines, leaving all optimizations to compiler. To test this package, you need to compile and install the library on your machine just like other C libraries. Make sure that you have a copy of [cmake](https://cmake.org) (>=3.12), then

```shell
mkdir build && cd build
cmake ..
make
make install
```

After installation, you can test the file [test_diff_mf1.c](./test/test_diff_mf1.c):

```shell
cc test_diff_mf1.c -limpf -o test_diff_mf1
./test_diff_mf1
```

By the way, the name of this library `impf` means *implicit functions*, as it was originally planned to specifically solve low-dimensional implicit functions (that's why you see a lot of routines for `2f2`) and later becomes a hodgepodge.

**Note:**

1. The suffix `_mfn` of all methods in this package means either "m functions of n variants" for problems of dimension <= 3, or "n-variant function that returns m-dimensionay array" for general problems, depending on specific propotype.
2. Generally, we follow the [coding style](https://www.kernel.org/doc/html/v4.10/process/coding-style.html). The only exception is when function's prototype is too long. In this case, we use 120 columns instead of 80 columns.

### Differentiation ([diff.h](./include/impf/diff.h))

- [impf_diff_1f1](./src/diff/1f1.c), [impf_diff_1f2](./src/diff/1f2.c), [impf_diff_1f3](./src/diff/1f3.c), [impf_diff_1fn](./src/diff/1fn.c)
- [impf_diff_2f1](./src/diff/2f1.c), [impf_diff_2f2](./src/diff/2f2.c), [impf_diff_2f3](./src/diff/2f3.c), [impf_diff_2fn](./src/diff/2fn.c)
- [impf_diff_3f1](./src/diff/3f1.c), [impf_diff_3f2](./src/diff/3f2.c), [impf_diff_3f3](./src/diff/3f3.c), [impf_diff_3fn](./src/diff/3fn.c), [impf_diff_mfn](./src/diff/mfn.c)

Example:

- [test_diff_mf1](./test/test_diff_mf1.c)(m = 1, 2, 3), [test_diff_mf2](./test/test_diff_mf2.c)(m = 1, 2, 3)
- [test_diff_mf3](./test/test_diff_mf3.c)(m = 1, 2, 3), [test_diff_mfn](./test/test_diff_mfn.c)(m = 1, 2, 3, m)

### Root Finding ([root.h](./include/impf/root.h))

- [impf_root_1f1_bisection](./src/root/1f1_bisection.c), [impf_root_1f1_secant](./src/root/1f1_secant.c), [impf_root_1f1_newton](./src/root/1f1_newton.c)
- [impf_root_2f2_newton](./src/root/2f2_newton.c), [impf_root_3f3_newton](./src/root/3f3_newton.c), [impf_root_nfn_newton](./src/root/nfn_newton.c)

Example:
- [test_root_1f1.c](test/test_root_1f1.c), [test_root_2f2.c](test/test_root_2f2.c), [test_root_3f3.c](test/test_root_3f3.c), [test_root_nfn.c](test/test_root_nfn.c)

### Linear Programming ([fmin.h](./include/impf/fmin.h) & [lp.h](./include/impf/fmin_lp.h))

- [impf_lp_readmps](./src/lp/readmps.c): import large scale LP model from MPS file
- [impf_lp_simplex](./src/lp/simplex_gen.c): simplex algorithm

Example:
- [test_lp_simplex_1.c](test/test_lp_simplex_1.c), [test_lp_simplex_2.c](test/test_lp_simplex_2.c), [test_lp_simplex_3.c](test/test_lp_simplex_3.c)
- [test_lp_simplex_4.c](test/test_lp_simplex_4.c), [test_lp_simplex_5.c](test/test_lp_simplex_5.c), [test_lp_simplex_6.c](test/test_lp_simplex_6.c)

Some [Netlib LP Benchmarks](https://www.netlib.org/lp/data/index.html) are tested in [test_lp_simplex_netlib.c](test/test_lp_simplex_netlib.c). Here are the testing results:

| testset  |  m   |  n   | dantzig | bland  | [pan97](https://doi.org/10.1016/S0898-1221(98)00127-8) |
| -------- | :--: | :--: | :-----: | :----: | :----------------------------------------------------: |
| 25fv47   | 821  | 1571 |    F    |   F    |                                                        |
| 80bau3b  | 2262 | 9799 |    F    |   F    |                                                        |
| adlittle |  56  |  97  |  5 ms   | 11 ms  |                                                        |
| afiro    |  27  |  32  | < 0 ms  | < 0 ms |                                                        |
| agg      | 488  | 163  | 297 ms  | 487 ms |                                                        |
| agg2     | 516  | 302  | 318 ms  | 291 ms |                                                        |
| agg3     | 516  | 302  |    F    | 473 ms |                                                        |

