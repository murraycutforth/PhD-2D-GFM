#ifndef EXTENSION_VFIELD_H
#define EXTENSION_VFIELD_H

#include "sim_info.hpp"
#include "GFM_ITM_interface.hpp"
#include "typedefs.hpp"

void extrapolate_extension_vfield
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	gridVector2dtype& velocities
);

#endif
