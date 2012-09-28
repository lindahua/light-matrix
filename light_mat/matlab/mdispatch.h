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

#include <marray.h>

namespace lmat { namespace matlab {


	/********************************************
	 *
	 *  type dispatch tools
	 *
	 ********************************************/

	template<class Args, template<typename T> class Runner>
	inline void dispatch_for_fptypes(mxClassID cid, Args& args,
			const char* type_err_id, const char* type_err_msg)
	{
		switch (cid)
		{
		case mxDOUBLE_CLASS:
			Runner<double>::run(args);
			break;

		case mxSINGLE_CLASS:
			Runner<float>::run(args);
			break;

		default:
			mexErrMsgIdAndTxt(type_err_id, type_err_msg);
		}
	}

	template<class Args, template<typename T> class Runner>
	inline void dispatch_for_numtypes(mxClassID cid, Args& args,
			const char* type_err_id, const char* type_err_msg)
	{
		switch (cid)
		{
		case mxDOUBLE_CLASS:
			Runner<double>::run(args);
			break;

		case mxSINGLE_CLASS:
			Runner<float>::run(args);
			break;

		case mxINT32_CLASS:
			Runner<int32_t>::run(args);
			break;

		case mxUINT32_CLASS:
			Runner<uint32_t>::run(args);
			break;

		case mxINT16_CLASS:
			Runner<int16_t>::run(args);
			break;

		case mxUINT16_CLASS:
			Runner<uint16_t>::run(args);
			break;

		case mxINT8_CLASS:
			Runner<int8_t>::run(args);
			break;

		case mxUINT8_CLASS:
			Runner<uint8_t>::run(args);
			break;

		case mxINT64_CLASS:
			Runner<int64_t>::run(args);
			break;

		case mxUINT64_CLASS:
			Runner<uint64_t>::run(args);
			break;

		default:
			mexErrMsgIdAndTxt(type_err_id, type_err_msg);
		}
	}

} }

#endif /* MDISPATCH_H_ */
