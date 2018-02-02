#include "FS_MUSCL.hpp"
#include "euler_misc.hpp"


void FS_MUSCL :: pure_fluid_update 
(
	const sim_info& params,
	const double gamma,
	const double pinf,
	const double dt,
	const BBrange& realcells,
	grideuler2type& grid,
	grideuler2type& future_grid
)
{
	const double dtodx = dt/params.dx;
	const double dtody = dt/params.dy;
	
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
			// Cell (i,j-1) denoted by L
			// Cell (i,j) denoted by C
			// Cell (i-1,j) denoted by B
			
			double beta = 1.0;
			for (int k=0; k<4; k++)
			{
				slope_C_x(k) = limited_slope(beta, grid[i][j](k) - grid[i][j-1](k), grid[i][j+1](k) - grid[i][j](k));
				slope_C_y(k) = limited_slope(beta, grid[i][j](k) - grid[i-1][j](k), grid[i+1][j](k) - grid[i][j](k));
				slope_L_x(k) = limited_slope(beta, grid[i][j-1](k) - grid[i][j-2](k), grid[i][j](k) - grid[i][j-1](k));
				slope_L_y(k) = limited_slope(beta, grid[i][j-1](k) - grid[i-1][j-1](k), grid[i+1][j-1](k) - grid[i][j-1](k));
				slope_B_x(k) = limited_slope(beta, grid[i-1][j](k) - grid[i-1][j-1](k), grid[i-1][j+1](k) - grid[i-1][j](k));
				slope_B_y(k) = limited_slope(beta, grid[i-1][j](k) - grid[i-2][j](k), grid[i][j](k) - grid[i-1][j](k));
			}
			
			C_BEV_L = grid[i][j] - 0.5*slope_C_x;
			C_BEV_R = grid[i][j] + 0.5*slope_C_x;
			C_BEV_B = grid[i][j] - 0.5*slope_C_y;
			C_BEV_T = grid[i][j] + 0.5*slope_C_y;
			B_BEV_L = grid[i-1][j] - 0.5*slope_B_x;
			B_BEV_R = grid[i-1][j] + 0.5*slope_B_x;
			B_BEV_B = grid[i-1][j] - 0.5*slope_B_y;
			B_BEV_T = grid[i-1][j] + 0.5*slope_B_y;
			L_BEV_L = grid[i][j-1] - 0.5*slope_L_x;
			L_BEV_R = grid[i][j-1] + 0.5*slope_L_x;
			L_BEV_B = grid[i][j-1] - 0.5*slope_L_y;
			L_BEV_T = grid[i][j-1] + 0.5*slope_L_y;
			
			C_BEV_L_evolved = C_BEV_L
					+ 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, C_BEV_L) - misc::flux_x_conserved(gamma, pinf, C_BEV_R))
					+ 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, C_BEV_B) - misc::flux_x_conserved(gamma, pinf, C_BEV_T));
					
			C_BEV_B_evolved = C_BEV_B
					+ 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, C_BEV_L) - misc::flux_x_conserved(gamma, pinf, C_BEV_R))
					+ 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, C_BEV_B) - misc::flux_x_conserved(gamma, pinf, C_BEV_T));
					
			L_BEV_R_evolved = L_BEV_R
					+ 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, L_BEV_L) - misc::flux_x_conserved(gamma, pinf, L_BEV_R))
					+ 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, L_BEV_B) - misc::flux_x_conserved(gamma, pinf, L_BEV_T));
					
			B_BEV_T_evolved = B_BEV_T
					+ 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, B_BEV_L) - misc::flux_x_conserved(gamma, pinf, B_BEV_R))
					+ 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, B_BEV_B) - misc::flux_x_conserved(gamma, pinf, B_BEV_T));
			
			if (! misc::is_physical_state(gamma, pinf, C_BEV_L_evolved) || ! misc::is_physical_state(gamma, pinf, L_BEV_R_evolved))
			{
				C_BEV_L_evolved = grid[i][j];
				L_BEV_R_evolved = grid[i][j-1];
			}
			if (! misc::is_physical_state(gamma, pinf, C_BEV_B_evolved) || ! misc::is_physical_state(gamma, pinf, B_BEV_T_evolved))
			{
				C_BEV_B_evolved = grid[i][j];
				B_BEV_T_evolved = grid[i-1][j];
			}
			
			RS->solve_RP(true, gamma, pinf, L_BEV_R_evolved, C_BEV_L_evolved, flux_x);
			RS->solve_RP(false, gamma, pinf, B_BEV_T_evolved, C_BEV_B_evolved, flux_y);
			
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




double FS_MUSCL :: limited_slope 
(
	double beta, 
	double del_L, 
	double del_R
) const
{
	double result;

	if (del_R > 0.0)
	{
		result = std::max(0.0, std::min(beta*del_L, del_R));
		result = std::max(result, std::min(del_L, beta*del_R));
	}
	else
	{
		result = std::min(0.0, std::max(beta*del_L, del_R));
		result = std::min(result, std::max(del_L, beta*del_R));
	}

	return result;
}
