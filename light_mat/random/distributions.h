/**
 * @file distributions.h
 *
 * Overall header for PRNGs for different distributions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DISTRIBUTIONS_H_
#define LIGHTMAT_DISTRIBUTIONS_H_

// discrete

#include <light_mat/random/uniform_int_distr.h>
#include <light_mat/random/bernoulli_distr.h>
#include <light_mat/random/binomial_distr.h>
#include <light_mat/random/geometric_distr.h>
#include <light_mat/random/discrete_distr.h>

// continuous

#include <light_mat/random/uniform_real_distr.h>
#include <light_mat/random/exponential_distr.h>

#endif 
