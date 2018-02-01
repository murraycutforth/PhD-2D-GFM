#include "FS_godunov.hpp"
#include "euler_misc.hpp"


void FS_godunov :: pure_fluid_update 
(
	const sim_info& params,
	const double gamma,
	const double pinf,
	const double dt,
	const BBrange& realcells,
	grideuler2type& grid,
	grideuler2type& future_grid
) const
{
	const double dtodx = dt/params.dx;
	const double dtody = dt/params.dy;
	vec4type flux_x (4);
	vec4type flux_y (4);
	
	for (int i=realcells.imin; i<realcells.imax; i++)
	{
		for (int j=realcells.jmin; j<realcells.jmax; j++)
		{
			future_grid[i][j] = grid[i][j];
			assert(misc::is_physical_state(gamma, pinf, grid[i][j]));
		}
	}
	
	for (int i=realcells.imin; i<realcells.imax + 1; i++)
	{
		for (int j=realcells.jmin; j<realcells.jmax + 1; j++)
		{			
			RS->solve_RP(true, gamma, pinf, grid[i][j-1], grid[i][j], flux_x);
			RS->solve_RP(false, gamma, pinf, grid[i-1][j], grid[i][j], flux_y);
			
			future_grid[i][j] += dtodx * flux_x + dtody * flux_y;
			future_grid[i][j-1] -= dtodx * flux_x;
			future_grid[i-1][j] -= dtody * flux_y;			
		}
	}
	
	for (int i=realcells.imin; i<realcells.imax; i++)
	{
		for (int j=realcells.jmin; j<realcells.jmax; j++)
		{
			grid[i][j] = future_grid[i][j];
		}
	}
}
