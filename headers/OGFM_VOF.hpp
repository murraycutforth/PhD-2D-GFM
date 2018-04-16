#ifndef OGFMVOF_H
#define OGFMVOF_H

#include "GFM_base.hpp"

class OGFM_VOF : public GFM_base {

	public:
	
	OGFM_VOF (const sim_info& params)
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
