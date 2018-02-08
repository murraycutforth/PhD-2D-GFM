#include "GFM_base.hpp"
#include "interpolation.hpp"

void GFM_base :: get_interpolated_mixedRPstates
(
	const sim_info& params, 
	const grideuler2type& prims1, 
	const grideuler2type& prims2,
	const GFM_ITM_interface& ls,
	const int i,
	const int j,
	const double ds,
	vec4type& Lprims,
	vec4type& Rprims
)
{
	Eigen::Vector2d CC = params.cellcentre_coord(i, j);
	double lsval = ls.get_sdf(i, j);
	Eigen::Vector2d normal = ls.get_normal(CC);
	double normmag = normal.norm();
	
	Eigen::Vector2d interfacelocation = CC - lsval * (normal / normmag);
	Eigen::Vector2d sampleposL = interfacelocation + ds * (normal / normmag);
	Eigen::Vector2d sampleposR = interfacelocation - ds * (normal / normmag);
	
	// sampleposL is in direction of positive normal, where phi>0 so fluid 2 is real
	
	Lprims = grid_bilinear_interpolate<vec4type>(params, prims2, sampleposL);
	Rprims = grid_bilinear_interpolate<vec4type>(params, prims1, sampleposR);
}
