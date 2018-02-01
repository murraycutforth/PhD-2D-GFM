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

	FS_base (std::shared_ptr<pure_RS_base> RS)
	:
		RS	(RS)
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
	) const =0;
	
};


#endif
