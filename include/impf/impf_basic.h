/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**************************************************************************//**
 * \file impf_basic.h
 * \brief Definitions of commonly used types and methods in IMPF.
 *****************************************************************************/

#ifndef __IMPF_BASIC_H__
#define __IMPF_BASIC_H__

#include <impf/impf_basic_blas.h>
#include <impf/impf_basic_rfs.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_global
 *
 * @{
 *****************************************************************************/

/** \brief The current IMPF version. */
#define __impf_VERSION__ "0.0.1"     

/** \brief The name of this package. */
#define __impf_LIBNAME__ "IMPF"      

/** \brief The current licence of this package. */
#define __impf_LICENCE__ "LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>"

/**
 * \brief (Possibly) the subroutine returns.
 *
 * \details Some subroutine returns a state as the execution condition. This is
 * mostly used as the return of iteration methods.
 */
typedef enum
{
    IMPF_EXIT_SUCCESS = 0,                    /**< 0 */
    IMPF_EXIT_FAILURE = 1,                    /**< 1 */
    IMPF_EXIT_ILLNESS = 2,                    /**< 2 */
} impf_t_subrtstate;

/** @} */  /* end of module_global */

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_BASIC_H__ */
