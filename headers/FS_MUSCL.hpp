#ifndef MUSCL_H
#define MUSCL_H


#include "FS_base.hpp"
#include <memory>


class FS_MUSCL : public FS_base {
	
	private:
	
	vec4type flux_x;
	vec4type flux_y;
	vec4type W_i_j;
	vec4type W_im_j;
	vec4type W_i_jm;
	vec4type W_imm_j;
	vec4type W_i_jmm;
	vec4type W_ip_j;
	vec4type W_i_jp;
	vec4type W_im_jm;
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
	matrix4d A_C;
	matrix4d A_L;
	matrix4d A_B;
	matrix4d B_C;
	matrix4d B_L;
	matrix4d B_B;
	
	double limited_slope 
	(
		double beta, 
		double del_L, 
		double del_R
	) const;
	
	vec4type minmod_limited_slope (const vec4type& diff_L, const vec4type& diff_R) const;


	public:

	FS_MUSCL (std::shared_ptr<pure_RS_base> RS)
	:
		FS_base	(RS),
		flux_x	(),
		flux_y	(),
		W_i_j (),
		W_im_j (),
		W_i_jm (),
		W_imm_j (),
		W_i_jmm (),
		W_ip_j (),
		W_i_jp (),
		W_im_jm (),
		slope_C_x (),
		slope_C_y (),
		slope_L_x (),
		slope_L_y (),
		slope_B_x (),
		slope_B_y (),
		C_BEV_L (),
		C_BEV_R (),
		C_BEV_B (),
		C_BEV_T (),
		B_BEV_L (),
		B_BEV_R (),
		B_BEV_B (),
		B_BEV_T (),
		L_BEV_L (),
		L_BEV_R (),
		L_BEV_B (),
		L_BEV_T (),
		C_BEV_L_evolved (),
		C_BEV_B_evolved (),
		L_BEV_R_evolved (),
		B_BEV_T_evolved (),
		A_C (),
		A_L (),
		A_B (),
		B_C (),
		B_L (),
		B_B ()
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
