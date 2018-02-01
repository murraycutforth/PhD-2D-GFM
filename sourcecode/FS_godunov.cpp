#include "FS_godunov.hpp"


void FS_godunov :: single_fluid_update 
(
	const sim_info& params,
	const double gamma,
	const double pinf,
	const double dt,
	const BBrange& realcells,
	const gridtype& grid,
	gridtype& future_grid
)
{
	const double dtodx = dt/oldstate.array.dx;
	const double dtody = dt/oldstate.array.dy;
	vectype flux_x (4);
	vectype flux_y (4);
	
	for (int i=realcells.imin; i<realcells.imax; i++)
	{
		for (int j=realcells.ymin; j<realcells.jmax; j++)
		{
			future_grid[i][j] = grid[i][j];
		}
	}
	
	for (int i=realcells.imin; i<realcells.imax + 1; i++)
	{
		for (int j=realcells.ymin; j<realcells.jmax + 1; j++)
		{
			RS->solve_RP(true, gamma, pinf, grid[i][j-1], grid[i][j], flux_x);
			RS->solve_RP(false, gamma, pinf, grid[i-1][j], grid[i][j], flux_y);
			
			future_grid[i][j] += dtodx * flux_x + dtody * flux_y;
			future_grid[i][j-1] -= dtodx * flux_x;
			future_grid[i-1][j] -= dtody * flux_y;
		}
	}
}
