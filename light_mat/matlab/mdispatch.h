/*
 * @file mdispatch.h
 *
 * @brief MATLAB runner dispatch (by type)
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MDISPATCH_H_
#define LIGHTMAT_MDISPATCH_H_

#include <light_mat/matlab/marray.h>

namespace lmat { namespace matlab {


	/********************************************
	 *
	 *  type dispatch tools
	 *
	 ********************************************/

	template<class Program>
	inline void dispatch_for_fptypes(Program& program)
	{
		switch (program.class_id())
		{
		case mxDOUBLE_CLASS:
			program.template run<double>();
			break;

		case mxSINGLE_CLASS:
			program.template run<float>();
			break;

		default:
			mexErrMsgIdAndTxt("light_mat:matlab:type_error",
					"Unexpected class for dispatch_for_fptypes");
		}
	}

	template<class Program>
	inline void dispatch_for_numtypes(Program& program)
	{
		switch (program.class_id())
		{
		case mxDOUBLE_CLASS:
			program.template run<double>();
			break;

		case mxSINGLE_CLASS:
			program.template run<float>();
			break;

		case mxINT32_CLASS:
			program.template run<int32_t>();
			break;

		case mxUINT32_CLASS:
			program.template run<uint32_t>();
			break;

		case mxINT16_CLASS:
			program.template run<int16_t>();
			break;

		case mxUINT16_CLASS:
			program.template run<uint16_t>();
			break;

		case mxINT8_CLASS:
			program.template run<int8_t>();
			break;

		case mxUINT8_CLASS:
			program.template run<uint8_t>();
			break;

		case mxINT64_CLASS:
			program.template run<int64_t>();
			break;

		case mxUINT64_CLASS:
			program.template run<uint64_t>();
			break;

		default:
			mexErrMsgIdAndTxt("light_mat:matlab:type_error",
					"Unexpected class for dispatch_for_numtypes");
		}
	}

} }

#endif /* MDISPATCH_H_ */
