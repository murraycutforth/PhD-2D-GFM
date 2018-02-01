#ifndef FS_BASE_H
#define FS_BASE_H


#include "pure_RS_base.hpp"
#include "sim_info.hpp"
#include "BBrange.hpp"
#include "typedefs.hpp"
#include <memory>


class FS_base {

	public:
	
	std::shared_ptr<pure_RS_base> RS;

	flow_solver_base (std::shared_ptr<pure_RS_base> RS)
	:
		RS	(RS)
	{}

	virtual void single_fluid_update 
	(
		const sim_info& params,
		const double gamma,
		const double pinf,
		const double dt,
		const BBrange& realcells,
		const gridtype& grid,
		gridtype& future_grid
	) =0;
	
};


#endif
