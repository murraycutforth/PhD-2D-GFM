#ifndef RIEMANNGFM_H
#define RIEMANNGFM_H

#include "GFM_base.hpp"

class GFM_riemann : public GFM_base {

	public:
	
	GFM_riemann (const sim_info& params)
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
