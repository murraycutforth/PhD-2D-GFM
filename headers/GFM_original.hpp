#ifndef OGFM_H
#define OGFM_H

#include "GFM_base.hpp"

class GFM_original : public GFM_base {

	public:
	
	GFM_original (const sim_info& params)
	:
		GFM_base (params)
	{}
	
	
	virtual void set_ghost_states
	(
		const sim_info& params, 
		const GFM_ITM_interface& ls,
		const binarySGparams& eosparams,
		const double t,
		grideuler2type& grid1, 
		grideuler2type& grid2,
		BBrange& realcells1,
		BBrange& realcells2
	) override;

};


#endif
