#include "FS_MUSCL.hpp"
#include "euler_misc.hpp"


//~ void FS_MUSCL :: pure_fluid_update 
//~ (
	//~ const sim_info& params,
	//~ const double gamma,
	//~ const double pinf,
	//~ const double dt,
	//~ const BBrange& realcells,
	//~ grideuler2type& grid,
	//~ grideuler2type& future_grid
//~ )
//~ {
	//~ const double dtodx = dt/params.dx;
	//~ const double dtody = dt/params.dy;
	
	//~ for (int i=realcells.imin; i<realcells.imax; i++)
	//~ {
		//~ for (int j=realcells.jmin; j<realcells.jmax; j++)
		//~ {
			//~ future_grid[i][j] = grid[i][j];
			//~ assert(misc::is_physical_state(gamma, pinf, grid[i][j]));
		//~ }
	//~ }
	
	//~ for (int i=realcells.imin; i<realcells.imax + 1; i++)
	//~ {
		//~ for (int j=realcells.jmin; j<realcells.jmax + 1; j++)
		//~ {
			//~ // Cell (i,j-1) denoted by L
			//~ // Cell (i,j) denoted by C
			//~ // Cell (i-1,j) denoted by B
			
			//~ double beta = 1.0;
			//~ for (int k=0; k<4; k++)
			//~ {
				//~ slope_C_x(k) = limited_slope(beta, grid[i][j](k) - grid[i][j-1](k), grid[i][j+1](k) - grid[i][j](k));
				//~ slope_C_y(k) = limited_slope(beta, grid[i][j](k) - grid[i-1][j](k), grid[i+1][j](k) - grid[i][j](k));
				//~ slope_L_x(k) = limited_slope(beta, grid[i][j-1](k) - grid[i][j-2](k), grid[i][j](k) - grid[i][j-1](k));
				//~ slope_L_y(k) = limited_slope(beta, grid[i][j-1](k) - grid[i-1][j-1](k), grid[i+1][j-1](k) - grid[i][j-1](k));
				//~ slope_B_x(k) = limited_slope(beta, grid[i-1][j](k) - grid[i-1][j-1](k), grid[i-1][j+1](k) - grid[i-1][j](k));
				//~ slope_B_y(k) = limited_slope(beta, grid[i-1][j](k) - grid[i-2][j](k), grid[i][j](k) - grid[i-1][j](k));
			//~ }
			
			//~ C_BEV_L = grid[i][j] - 0.5*slope_C_x;
			//~ C_BEV_R = grid[i][j] + 0.5*slope_C_x;
			//~ C_BEV_B = grid[i][j] - 0.5*slope_C_y;
			//~ C_BEV_T = grid[i][j] + 0.5*slope_C_y;
			//~ B_BEV_L = grid[i-1][j] - 0.5*slope_B_x;
			//~ B_BEV_R = grid[i-1][j] + 0.5*slope_B_x;
			//~ B_BEV_B = grid[i-1][j] - 0.5*slope_B_y;
			//~ B_BEV_T = grid[i-1][j] + 0.5*slope_B_y;
			//~ L_BEV_L = grid[i][j-1] - 0.5*slope_L_x;
			//~ L_BEV_R = grid[i][j-1] + 0.5*slope_L_x;
			//~ L_BEV_B = grid[i][j-1] - 0.5*slope_L_y;
			//~ L_BEV_T = grid[i][j-1] + 0.5*slope_L_y;
			
			//~ C_BEV_L_evolved = C_BEV_L
					//~ + 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, C_BEV_L) - misc::flux_x_conserved(gamma, pinf, C_BEV_R))
					//~ + 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, C_BEV_B) - misc::flux_x_conserved(gamma, pinf, C_BEV_T));
					
			//~ C_BEV_B_evolved = C_BEV_B
					//~ + 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, C_BEV_L) - misc::flux_x_conserved(gamma, pinf, C_BEV_R))
					//~ + 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, C_BEV_B) - misc::flux_x_conserved(gamma, pinf, C_BEV_T));
					
			//~ L_BEV_R_evolved = L_BEV_R
					//~ + 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, L_BEV_L) - misc::flux_x_conserved(gamma, pinf, L_BEV_R))
					//~ + 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, L_BEV_B) - misc::flux_x_conserved(gamma, pinf, L_BEV_T));
					
			//~ B_BEV_T_evolved = B_BEV_T
					//~ + 0.5*dtodx*(misc::flux_x_conserved(gamma, pinf, B_BEV_L) - misc::flux_x_conserved(gamma, pinf, B_BEV_R))
					//~ + 0.5*dtody*(misc::flux_x_conserved(gamma, pinf, B_BEV_B) - misc::flux_x_conserved(gamma, pinf, B_BEV_T));
			
			//~ if (! misc::is_physical_state(gamma, pinf, C_BEV_L_evolved) || ! misc::is_physical_state(gamma, pinf, L_BEV_R_evolved))
			//~ {
				//~ C_BEV_L_evolved = grid[i][j];
				//~ L_BEV_R_evolved = grid[i][j-1];
			//~ }
			//~ if (! misc::is_physical_state(gamma, pinf, C_BEV_B_evolved) || ! misc::is_physical_state(gamma, pinf, B_BEV_T_evolved))
			//~ {
				//~ C_BEV_B_evolved = grid[i][j];
				//~ B_BEV_T_evolved = grid[i-1][j];
			//~ }
			
			//~ RS->solve_RP(true, gamma, pinf, L_BEV_R_evolved, C_BEV_L_evolved, flux_x);
			//~ RS->solve_RP(false, gamma, pinf, B_BEV_T_evolved, C_BEV_B_evolved, flux_y);
			
			//~ future_grid[i][j] += dtodx * flux_x + dtody * flux_y;
			//~ future_grid[i][j-1] -= dtodx * flux_x;
			//~ future_grid[i-1][j] -= dtody * flux_y;
		//~ }
	//~ }
	
	//~ for (int i=realcells.imin; i<realcells.imax; i++)
	//~ {
		//~ for (int j=realcells.jmin; j<realcells.jmax; j++)
		//~ {
			//~ grid[i][j] = future_grid[i][j];
		//~ }
	//~ }
//~ }




//~ double FS_MUSCL :: limited_slope 
//~ (
	//~ double beta, 
	//~ double del_L, 
	//~ double del_R
//~ ) const
//~ {
	//~ double result;

	//~ if (del_R > 0.0)
	//~ {
		//~ result = std::max(0.0, std::min(beta*del_L, del_R));
		//~ result = std::max(result, std::min(del_L, beta*del_R));
	//~ }
	//~ else
	//~ {
		//~ result = std::min(0.0, std::max(beta*del_L, del_R));
		//~ result = std::min(result, std::max(del_L, beta*del_R));
	//~ }

	//~ return result;
//~ }


















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
			
			W_i_j = misc::conserved_to_primitives(gamma, pinf, grid[i][j]);
			W_ip_j = misc::conserved_to_primitives(gamma, pinf, grid[i+1][j]);
			W_im_j = misc::conserved_to_primitives(gamma, pinf, grid[i-1][j]);
			W_imm_j = misc::conserved_to_primitives(gamma, pinf, grid[i-2][j]);
			W_i_jp = misc::conserved_to_primitives(gamma, pinf, grid[i][j+1]);
			W_i_jm = misc::conserved_to_primitives(gamma, pinf, grid[i][j-1]);
			W_i_jmm = misc::conserved_to_primitives(gamma, pinf, grid[i][j-2]);
			W_im_jm = misc::conserved_to_primitives(gamma, pinf, grid[i-1][j-1]);
			
			
			slope_C_x = minmod_limited_slope(W_i_j - W_i_jm, W_i_jp - W_i_j);
			slope_C_y = minmod_limited_slope(W_i_j - W_im_j, W_ip_j - W_i_j);
			
			slope_L_x = minmod_limited_slope(W_i_jm - W_i_jmm, W_i_j - W_i_jm);
			slope_L_y = minmod_limited_slope(W_i_jm - W_im_jm, misc::conserved_to_primitives(gamma, pinf, grid[i+1][j-1]) - W_i_jm);
			
			slope_B_x = minmod_limited_slope(W_im_j - W_im_jm, misc::conserved_to_primitives(gamma, pinf, grid[i-1][j+1]) - W_im_j);
			slope_B_y = minmod_limited_slope(W_im_j - W_imm_j, W_i_j - W_im_j);
			
			
			C_BEV_L = W_i_j - 0.5*slope_C_x;
			C_BEV_R = W_i_j + 0.5*slope_C_x;
			C_BEV_B = W_i_j - 0.5*slope_C_y;
			C_BEV_T = W_i_j + 0.5*slope_C_y;
			B_BEV_L = W_im_j - 0.5*slope_B_x;
			B_BEV_R = W_im_j + 0.5*slope_B_x;
			B_BEV_B = W_im_j - 0.5*slope_B_y;
			B_BEV_T = W_im_j + 0.5*slope_B_y;
			L_BEV_L = W_i_jm - 0.5*slope_L_x;
			L_BEV_R = W_i_jm + 0.5*slope_L_x;
			L_BEV_B = W_i_jm - 0.5*slope_L_y;
			L_BEV_T = W_i_jm + 0.5*slope_L_y;
			
			
			if (!misc::is_physical_primitives(gamma, pinf, C_BEV_L) || !misc::is_physical_primitives(gamma, pinf, C_BEV_R))
			{
				C_BEV_L = W_i_j;
				C_BEV_R = W_i_j;
			}
			
			if (!misc::is_physical_primitives(gamma, pinf, C_BEV_B) || !misc::is_physical_primitives(gamma, pinf, C_BEV_T))
			{
				C_BEV_B = W_i_j;
				C_BEV_T = W_i_j;
			}
			
			if (!misc::is_physical_primitives(gamma, pinf, B_BEV_L) || !misc::is_physical_primitives(gamma, pinf, B_BEV_R))
			{
				B_BEV_L = W_im_j;
				B_BEV_R = W_im_j;
			}
			
			if (!misc::is_physical_primitives(gamma, pinf, B_BEV_B) || !misc::is_physical_primitives(gamma, pinf, B_BEV_T))
			{
				B_BEV_B = W_im_j;
				B_BEV_T = W_im_j;
			}
			
			if (!misc::is_physical_primitives(gamma, pinf, L_BEV_L) || !misc::is_physical_primitives(gamma, pinf, L_BEV_R))
			{
				L_BEV_L = W_i_jm;
				L_BEV_R = W_i_jm;
			}
			
			if (!misc::is_physical_primitives(gamma, pinf, L_BEV_B) || !misc::is_physical_primitives(gamma, pinf, L_BEV_T))
			{
				L_BEV_B = W_i_jm;
				L_BEV_T = W_i_jm;
			}
			
			
			A_C = misc::quasilinear_form_A(gamma, pinf, W_i_j);
			B_C = misc::quasilinear_form_B(gamma, pinf, W_i_j);
			A_L = misc::quasilinear_form_A(gamma, pinf, W_i_jm);
			B_L = misc::quasilinear_form_B(gamma, pinf, W_i_jm);
			A_B = misc::quasilinear_form_A(gamma, pinf, W_im_j);
			B_B = misc::quasilinear_form_B(gamma, pinf, W_im_j);
			
			
			C_BEV_L_evolved = C_BEV_L
					+ 0.5 * dtodx * A_C * (C_BEV_L - C_BEV_R)
					+ 0.5 * dtody * B_C * (C_BEV_B - C_BEV_T);
					
			C_BEV_B_evolved = C_BEV_B 
					+ 0.5 * dtodx * A_C * (C_BEV_L - C_BEV_R)
					+ 0.5 * dtody * B_C * (C_BEV_B - C_BEV_T);
					
			L_BEV_R_evolved = L_BEV_R
					+ 0.5 * dtodx * A_L * (L_BEV_L - L_BEV_R)
					+ 0.5 * dtody * B_L * (L_BEV_B - L_BEV_T);
					
			B_BEV_T_evolved = B_BEV_T
					+ 0.5 * dtodx * A_B * (B_BEV_L - B_BEV_R)
					+ 0.5 * dtody * B_B * (B_BEV_B - B_BEV_T);
			
			
			
			if (! misc::is_physical_primitives(gamma, pinf, C_BEV_L_evolved) || ! misc::is_physical_primitives(gamma, pinf, L_BEV_R_evolved))
			{
				C_BEV_L_evolved = W_i_j;
				L_BEV_R_evolved = W_i_jm;
			}
			if (! misc::is_physical_primitives(gamma, pinf, C_BEV_B_evolved) || ! misc::is_physical_primitives(gamma, pinf, B_BEV_T_evolved))
			{
				C_BEV_B_evolved = W_i_j;
				B_BEV_T_evolved = W_im_j;
			}
			
			RS->solve_RP(true, gamma, pinf, misc::primitives_to_conserved(gamma, pinf, L_BEV_R_evolved), misc::primitives_to_conserved(gamma, pinf, C_BEV_L_evolved), flux_x);
			RS->solve_RP(false, gamma, pinf, misc::primitives_to_conserved(gamma, pinf, B_BEV_T_evolved), misc::primitives_to_conserved(gamma, pinf, C_BEV_B_evolved), flux_y);
			
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
			// assert(misc::is_physical_state(gamma, pinf, grid[i][j]));
		}
	}
}





vec4type FS_MUSCL :: minmod_limited_slope (const vec4type& diff_L, const vec4type& diff_R) const
{
	vec4type slope;
	
	for (int k=0; k<4; k++)
	{
		if (diff_R(k) > 0.0)
		{
			slope(k) = std::max(0.0, std::min(diff_L(k), diff_R(k)));
		}
		else
		{
			slope(k) = std::min(0.0, std::max(diff_L(k), diff_R(k)));
		}
	}
	
	return slope;
}
