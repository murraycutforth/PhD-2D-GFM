#ifndef GODUNOV_H
#define GODUNOV_H


#include "FS_base.hpp"
#include <memory>


class FS_godunov : public FS_base {

	public:

	FS_godunov (std::shared_ptr<pure_RS_base> RS)
	:
		FS_base	(RS)
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
	);
	
};

#endif
