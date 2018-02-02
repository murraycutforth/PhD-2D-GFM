#ifndef MUSCL_H
#define MUSCL_H


#include "FS_base.hpp"
#include <memory>


class FS_MUSCL : public FS_base {
	
	private:
	
	vec4type flux_x;
	vec4type flux_y;
	vec4type slope_C_x;
	vec4type slope_C_y;
	vec4type slope_L_x;
	vec4type slope_L_y;
	vec4type slope_B_x;
	vec4type slope_B_y;
	vec4type C_BEV_L;
	vec4type C_BEV_R;
	vec4type C_BEV_B;
	vec4type C_BEV_T;
	vec4type B_BEV_L;
	vec4type B_BEV_R;
	vec4type B_BEV_B;
	vec4type B_BEV_T;
	vec4type L_BEV_L;
	vec4type L_BEV_R;
	vec4type L_BEV_B;
	vec4type L_BEV_T;
	vec4type C_BEV_L_evolved;
	vec4type C_BEV_B_evolved;
	vec4type L_BEV_R_evolved;
	vec4type B_BEV_T_evolved;
	
	double limited_slope 
	(
		double beta, 
		double del_L, 
		double del_R
	) const;


	public:

	FS_MUSCL (std::shared_ptr<pure_RS_base> RS)
	:
		FS_base	(RS),
		flux_x	(4),
		flux_y	(4),
		slope_C_x (4),
		slope_C_y (4),
		slope_L_x (4),
		slope_L_y (4),
		slope_B_x (4),
		slope_B_y (4),
		C_BEV_L (4),
		C_BEV_R (4),
		C_BEV_B (4),
		C_BEV_T (4),
		B_BEV_L (4),
		B_BEV_R (4),
		B_BEV_B (4),
		B_BEV_T (4),
		L_BEV_L (4),
		L_BEV_R (4),
		L_BEV_B (4),
		L_BEV_T (4),
		C_BEV_L_evolved (4),
		C_BEV_B_evolved (4),
		L_BEV_R_evolved (4),
		B_BEV_T_evolved (4)
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
	) override;
	
};

#endif
