#ifndef GFM_BASE_H
#define GFM_BASE_H

#include "velocity_field.hpp"
#include "velocity_field_firstorder.hpp"
#include "velocity_field_invdist.hpp"
#include "mixed_RS_base.hpp"
#include "sim_info.hpp"
#include "GFM_ITM_interface.hpp"
#include "typedefs.hpp"
#include "stiffened_gas_eos.hpp"
#include "mixed_RS_exact.hpp"
#include "BBrange.hpp"
#include <memory>

class GFM_base {

	public:
	
	bool use_extension_vfield;
	
	std::shared_ptr<solved_velocity_field_base> vfield;
	
	std::shared_ptr<mixed_RS_base> mixed_RS;
	

	GFM_base (const sim_info& params)
	:
		use_extension_vfield (true),
		vfield	(std::make_shared<solved_velocity_field_firstorder>(params)),
		mixed_RS	(std::make_shared<mixed_RS_exact>())
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
	) =0;
	
	
	void get_interpolated_mixedRPstates
	(
		const sim_info& params, 
		const grideuler2type& prims1, 
		const grideuler2type& prims2,
		const GFM_ITM_interface& ls,
		const int i,
		const int j,
		const double ds,
		vec4type& Lstate,
		vec4type& Rstate
	);
	
	void get_interpolated_mixedRPstates_VOF
	(
		const sim_info& params, 
		const grideuler2type& prims1, 
		const grideuler2type& prims2,
		const Eigen::Vector2d normal,
		const Eigen::Vector2d interfacelocation,
		const double ds,
		vec4type& Lstate,
		vec4type& Rstate
	);

};


#endif
