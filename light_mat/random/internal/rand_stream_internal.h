/**
 * @file rand_stream_internal.h
 *
 * @brief Internal helper for random stream implementation
 *
 * @author Dahua Lin
 */

#ifndef RAND_STREAM_INTERNAL_H_
#define RAND_STREAM_INTERNAL_H_

#include <light_mat/random/rand_stream.h>
#include <cstring>

namespace lmat { namespace random { namespace internal {

	template<class State, typename T, typename BDTag>
	void gen_rand_seq(State& s, stream_tracker<T>& trk, BDTag bd_tag, void *buf, size_t nbytes)
	{
		// pre-condition: nbytes > 0

		const char *ps = reinterpret_cast<const char*>(s.ptr_base());
		char *pd = reinterpret_cast<char*>(buf);

		const size_t sbytes = trk.length() * sizeof(T);

		trk.to_boundary(bd_tag);

		if (!trk.is_end()) // copy from remaining stream
		{
			const size_t rbytes = trk.remain() * sizeof(T);

			if (nbytes <= rbytes)
			{
				std::memcpy(pd, ps, nbytes);
				nbytes = 0;
				trk.forward_bytes(nbytes);
			}
			else
			{
				std::memcpy(pd, ps, rbytes);
				nbytes -= rbytes;
				pd += rbytes;
				trk.set_end();
			}
		}

		// if nbytes > 0 at this point, then tracker.is_end()

		if (nbytes >= sbytes)  // whole-state copying (without tracking)
		{
			size_t m = nbytes / sbytes;
			for (size_t i = 0; i < m; ++i)
			{
				s.next();
				std::memcpy(pd, ps, sbytes);
			}

			size_t l = m * sbytes;
			nbytes -= l;
			pd += l;
		}

		if (nbytes)  // process remaining
		{
			s.next();
			s.rewind();

			std::memcpy(pd, ps, nbytes);
			trk.forward_bytes(nbytes);
		}
	}

} } }

#endif /* RAND_STREAM_INTERNAL_H_ */
