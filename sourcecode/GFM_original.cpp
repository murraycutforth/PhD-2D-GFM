#include "GFM_original.hpp"
#include "extrapolate_extension_vfield.hpp"
#include "extrapolate_vector.hpp"
#include "euler_misc.hpp"



void GFM_original :: set_ghost_states
(
	const sim_info& params, 
	const GFM_ITM_interface& ls,
	const binarySGparams& eosparams,
	const double t,
	grideuler2type& grid1, 
	grideuler2type& grid2,
	BBrange& realcells1,
	BBrange& realcells2
)
{
	// Slip BC is not implemented for this GFM
	assert(use_slip_BCs == false);
	
	static gridVector2dtype newvelocities (params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
	static grideuler2type prims1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type prims2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	
	if (int(newvelocities.size()) != params.Ny + 2 * params.numGC)
	{
		newvelocities = gridVector2dtype(params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
		prims1 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
		prims2 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	}
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			assert(misc::is_physical_state(eosparams.gamma1, eosparams.pinf1, grid1[i][j]));
			assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
			
			newvelocities[i][j] = Eigen::Vector2d::Zero();
			prims1[i][j] = misc::conserved_to_primitives(eosparams.gamma1, eosparams.pinf1, grid1[i][j]);
			prims2[i][j] = misc::conserved_to_primitives(eosparams.gamma2, eosparams.pinf2, grid2[i][j]);
		}
	}
	
	realcells1.imin = params.numGC + params.Ny;
	realcells1.imax = params.numGC;
	realcells1.jmin = params.numGC + params.Nx;
	realcells1.jmax = params.numGC;
	
	realcells2.imin = params.numGC + params.Ny;
	realcells2.imax = params.numGC;
	realcells2.jmin = params.numGC + params.Nx;
	realcells2.jmax = params.numGC;
	
	extrapolate_vector(params, ls, 6, prims1, prims2);
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double lsval = ls.get_sdf(i, j);
			
			if (lsval <= 0.0)
			{
				if (ls.is_interfacial_cell(i,j))
				{
					newvelocities[i][j](0) = prims1[i][j](1);
					newvelocities[i][j](1) = prims1[i][j](2);
				}
				
				double u_new = prims1[i][j](1); 
				double v_new = prims1[i][j](2);
				double p_new = prims1[i][j](3);
				double rho_new = eos::isentropic_extrapolation(eosparams.gamma2, eosparams.pinf2, prims2[i][j](0), prims2[i][j](3), p_new);
				
				// assert(p_new >= 0.0);
				assert(rho_new >= 0.0);
				
				prims2[i][j](0) = rho_new;
				prims2[i][j](1) = u_new;
				prims2[i][j](2) = v_new;
				prims2[i][j](3) = p_new;
				
				realcells1.imin = std::min(realcells1.imin, i-2);
				realcells1.imax = std::max(realcells1.imax, i+2);
				realcells1.jmin = std::min(realcells1.jmin, j-2);
				realcells1.jmax = std::max(realcells1.jmax, j+2);
			}
			else
			{
				if (ls.is_interfacial_cell(i,j))
				{
					newvelocities[i][j](0) = prims2[i][j](1);
					newvelocities[i][j](1) = prims2[i][j](2);
				}
				
				double u_new = prims2[i][j](1); 
				double v_new = prims2[i][j](2);
				double p_new = prims2[i][j](3);
				double rho_new = eos::isentropic_extrapolation(eosparams.gamma1, eosparams.pinf1, prims1[i][j](0), prims1[i][j](3), p_new);
				
				// assert(p_new >= 0.0);
				assert(rho_new >= 0.0);
				
				prims1[i][j](0) = rho_new;
				prims1[i][j](1) = u_new;
				prims1[i][j](2) = v_new;
				prims1[i][j](3) = p_new;
				
				realcells2.imin = std::min(realcells2.imin, i-2);
				realcells2.imax = std::max(realcells2.imax, i+2);
				realcells2.jmin = std::min(realcells2.jmin, j-2);
				realcells2.jmax = std::max(realcells2.jmax, j+2);
			}
		}
	}
	
	realcells1.imin = std::max(realcells1.imin, params.numGC);
	realcells1.imax = std::min(realcells1.imax, params.numGC + params.Ny);
	realcells1.jmin = std::max(realcells1.jmin, params.numGC);
	realcells1.jmax = std::min(realcells1.jmax, params.numGC + params.Nx);
	
	realcells2.imin = std::max(realcells2.imin, params.numGC);
	realcells2.imax = std::min(realcells2.imax, params.numGC + params.Ny);
	realcells2.jmin = std::max(realcells2.jmin, params.numGC);
	realcells2.jmax = std::min(realcells2.jmax, params.numGC + params.Nx);
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			grid1[i][j] = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, prims1[i][j]);
			grid2[i][j] = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, prims2[i][j]);
			assert(misc::is_physical_state(eosparams.gamma1, eosparams.pinf1, grid1[i][j]));
			assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
		}
	}
	
	if (use_extension_vfield)
	{ 
		extrapolate_extension_vfield(params, ls, 20, newvelocities);
	}
	else
	{
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				double lsval = ls.get_sdf(i, j);
				
				if (lsval <= 0.0)
				{
					newvelocities[i][j] << prims1[i][j](1), prims1[i][j](2);
				}
				else
				{
					newvelocities[i][j] << prims2[i][j](1), prims2[i][j](2);
				}
			}
		}
	}
		
	vfield->store_velocity_field(newvelocities, t);
}
