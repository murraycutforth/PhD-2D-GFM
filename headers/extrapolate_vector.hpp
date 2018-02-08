#ifndef EXTRAPOLATE_VECTOR_H
#define EXTRAPOLATE_VECTOR_H

#include "sim_info.hpp"
#include "GFM_ITM_interface.hpp"
#include "typedefs.hpp"

void extrapolate_vector
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	grideuler2type& grid1, 
	grideuler2type& grid2
);

void extrapolate_vector_mgfm
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	grideuler2type& grid1, 
	grideuler2type& grid2
);

#endif
