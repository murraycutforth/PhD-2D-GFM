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

	virtual void pure_fluid_update 
	(
		const sim_info& params,
		const double gamma,
		const double pinf,
		const double dt,
		const BBrange& realcells,
		grideuler2type& grid,
		grideuler2type& future_grid
	) const override;
	
};

#endif
