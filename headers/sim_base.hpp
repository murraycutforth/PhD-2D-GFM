#ifndef SIM_BASE_H
#define SIM_BASE_H

#include "GFM_settingsfile.hpp"

class sim_base {

	public:

	virtual void run_sim (GFM_settingsfile SF) =0;
};


#endif
