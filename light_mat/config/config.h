/**
 * @file config.h
 *
 * @brief The configuration file for Basic Computation Supporting Library
 *
 * @author dhlin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_CONFIG_H_
#define LIGHTMAT_CONFIG_H_

#include <light_mat/config/user_config.h>
#include <light_mat/config/platform_config.h>

#if (LMAT_DIAGNOSIS_LEVEL >= 3)
#define LMAT_ENABLE_INDEX_CHECKING
#endif

#if (LMAT_DIAGNOSIS_LEVEL >= 1)
#define LMAT_ENABLE_DIM_CHECKING
#endif


#endif

