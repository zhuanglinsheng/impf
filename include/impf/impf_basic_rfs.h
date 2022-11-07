/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**************************************************************************//**
 * \file impf_basic_rfs.h
 * \brief Definitions of commonly used functional types.
 *****************************************************************************/

#ifndef __IMPF_BASIC_RFS_H__
#define __IMPF_BASIC_RFS_H__

#include <impf/impf_basic_blas.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_ftypes_realfunc
 * \details 
 * \anchor realfunc 
 * **Real functions** refers to those user-defined functions with real inputs
 * and return that IMPF recognizes.
 *
 * Note that IMPF is a "real" numerical package which focus on `double` float 
 * values. \ref realfunc "Real functions" have `double` inputs and `double` 
 * outputs. The input could be a double value (scalar), several double values, 
 * a double array (pointer), a double vector (pointer), or several double 
 * arrays/vectors. The return could be a double scalar or a double 
 * array/vector.
 * 
 * The \ref realfunc "real functions" of 1-variant, 2-variant, 3-variant and 
 * n-variants that IMPF recognizes are defined here.
 *****************************************************************************/
/**@{*/

/**
 * \brief \ref realfunc "Real function" of 1 variant and 1 return.
 */
typedef double (*impf_t_f1)(const double x);

/**
 * \brief \ref realfunc "Real function" of 2 variants and 1 return.
 */
typedef double (*impf_t_f2)(const double x, const double y);

/**
 * \brief \ref realfunc "Real function" of 3 variants and 1 return.
 */
typedef double (*impf_t_f3)(const double x, const double y, const double z);

/**
 * \brief \ref realfunc "Real function" of n variable and 1 return.
 */
typedef double (*impf_t_fn)(const double * x, const unsigned int n);

/**
 * \brief \ref realfunc "Real function" of n variants and m returns,
 *        \f$y=f(x)\f$ where \f$f:\mathbb{R}^n\to\mathbb{R}^m\f$. 
 * \param [in] m the length of `y`.
 * \param [in] n the length of `x`.
 * \param [in] x the evaluation point.
 * \param [out] y the evaluation.
 */
typedef void (*impf_t_mfn_array)(
    const unsigned int m, const unsigned int n, const double * x, double * y);

/**
 * \brief \ref realfunc "Real function" of n variants and m returns,
 *        \f$y=f(x)\f$ where \f$f:\mathbb{R}^n\to\mathbb{R}^m\f$. 
 * \param [in] x the evaluation point.
 * \param [out] y the evaluation.
 */
typedef void (*impf_t_mfn_vector)(
    const impf_t_vector * x, impf_t_vector * y);
/**@}*/


/**************************************************************************//**
 * \addtogroup module_ftypes_wrapper
 * \details 
 * \anchor wrapper
 * **Wrapper** refers to those user-defined subroutines that takes a
 * \ref realfunc "real function" as an (usually the first) input.
 *
 * \ref wrapper "Wrappers" are often used to represent the outer function that 
 * defines the implicit function. Mathematically, wrappers are functions 
 * \f$F\f$ of the form 
 * \f$F(y(x), x): \mathbb{R}^m\times\mathbb{R}^n\to\mathbb{R}^m\f$, where 
 * \f$y(x): \mathbb{R}^n\to\mathbb{R}^m\f$ is implicitly defined by the
 * functional equation \f$0 = F(y(x), x)\f$. 
 * 
 * The \ref wrapper "wrappers" of 1-variant, 2-variants, 3-variants
 * and n-variants \ref realfunc "real functions" that IMPF recognizes are 
 * defined here. 
 *****************************************************************************/
/**@{*/

/**
 * \brief Wrapper of a \ref realfunc "real function" of 1 variant.
 * \param [in] f the implicit function.
 * \param [in] x the evaluation point.
 * \param [out] re the evaluation of F, a duble scalar.
 * 
 * Example file `example/def/eg_def_wrp_1f1.c`
 * \include example/def/eg_def_wrp_1f1.c
 */
typedef void (*impf_t_1f1_wrapper)(
    const impf_t_f1 f, const double x, double * re);

/**
 * \brief Wrapper of a \ref realfunc "real function" of 2 variants.
 * \param [in] f the implicit function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [out] re the evaluation of F, a duble scalar.
 */
typedef void (*impf_t_1f2_wrapper)(
    const impf_t_f2 f, const double x, const double y, double * re);

/**
 * \brief Wrapper of a \ref realfunc "real function" of 3 variants.
 * \param [in] f the implicit function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [in] z the 3rd scalar of the evaluation point.
 * \param [out] re the evaluation of F, a duble scalar.
 */
typedef void (*impf_t_1f3_wrapper)(
    const impf_t_f3 f, const double x, const double y, const double z, double * re);

/**
 * \brief Wrapper of 2 \ref realfunc "real functions" of 1 variant.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] x the evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 2.
 */
typedef void (*impf_t_2f1_wrapper)(
    const impf_t_f1 f1, const impf_t_f1 f2, const double x,  double * re);

/**
 * \brief Wrapper of 2 \ref realfunc "real functions" of 2 variants.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 2.
 */
typedef void (*impf_t_2f2_wrapper)(
    const impf_t_f2 f1, const impf_t_f2 f2,
    const double x, const double y, double * re);

/**
 * \brief Wrapper of 2 \ref realfunc "real functions" of 3 variants.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [in] z the 3rd scalar of the evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 2.
 */
typedef void (*impf_t_2f3_wrapper)(
    const impf_t_f3 f1, const impf_t_f3 f2,
    const double x, const double y, const double z, double * re);

/**
 * \brief Wrapper of 3 \ref realfunc "real functions" of 1 variant.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] f3 the 3rd implicit function. 
 * \param [in] x the evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 3.
 */
typedef void (*impf_t_3f1_wrapper)(
    const impf_t_f1 f1, const impf_t_f1 f2, const impf_t_f1 f3,
    const double x, double * re);

/**
 * \brief Wrapper of 3 \ref realfunc "real functions" of 2 variants.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] f3 the 3rd implicit function. 
 * \param [in] x the 1st scalar of evaluation point.
 * \param [in] y the 2nd scalar of evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 3.
 */
typedef void (*impf_t_3f2_wrapper)(
    const impf_t_f2 f1, const impf_t_f2 f2, const impf_t_f2 f3,
    const double x, const double y, double * re);

/**
 * \brief Wrapper of 3 \ref realfunc "real functions" of 3 variants.
 * \param [in] f1 the 1st implicit function.
 * \param [in] f2 the 2nd implicit function.
 * \param [in] f3 the 3rd implicit function. 
 * \param [in] x the 1st scalar of evaluation point.
 * \param [in] y the 2nd scalar of evaluation point.
 * \param [in] z the 3rd scalar of evaluation point.
 * \param [out] re the evaluation of F, a duble array of length = 3.
 */
typedef void (*impf_t_3f3_wrapper)(
    const impf_t_f3 f1, const impf_t_f3 f2, const impf_t_f3 f3,
    const double x, const double y, const double z, double * re);

/**
 * \brief Wrapper of m \ref realfunc "real functions" of n variants. 
 * \param [in] mfn the implicit function of n inputs and m returns.
 * \param [in] m the number of returns. 
 * \param [in] n the number of inputs.
 * \param [in] x the evaluation point.
 * \param buf
 * \param [out] re the evaluation of F, a duble array of length = m.
 */
typedef void (*impf_t_mfn_array_wrapper)(
    const impf_t_mfn_array mfn, const unsigned int m, const unsigned int n,
    const double * x, double * buf, double * re);

/**
 * \brief Wrapper of m \ref realfunc "real functions" of n variants. 
 * \param [in] mfn the implicit function of n inputs and m returns.
 * \param [in] m the number of returns. 
 * \param [in] x the evaluation point.
 * \param buf the calculation buffer.
 * \param [out] re the evaluation of F, a duble array of length = m.
 */
typedef void (*impf_t_mfn_vector_wrapper)(
    const impf_t_mfn_vector mfn, const unsigned int m,
    const impf_t_vector * x, double * buf, double * re);
/**@}*/


#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_BASIC_RFS_H__ */
