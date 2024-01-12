/*
* Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
* License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
*/
#ifndef __IMPF__MAGICS_H__
#define __IMPF__MAGICS_H__

#include "utils.h"
#include "fmin.h"

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/*
 * The general formulation of quadratic programming
 * 
 * min  x'Qx / 2 + c'x
 * s.t. Ax <= b
 */

/* 
 * The standard formulation of quadratic programming
 *
 * min  x'Qx / 2 + c'x
 * s.t. Ax = b
 *
 * whose Lagrangian is 
 *
 * L = x'Qx / 2 + c'x + lambda'(Ax - b)
 *
 * where `lambda` is a vector
 */



#ifdef __cpluscplus
}
#endif /* __cpluscplus */

