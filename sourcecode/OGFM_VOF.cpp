#include "OGFM_VOF.hpp"
#include "euler_misc.hpp"
#include "euler_bc.hpp"
#include "GFM_VOF_interface.hpp"



void OGFM_VOF :: set_ghost_states
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
	
	const GFM_VOF_interface& vof = dynamic_cast<const GFM_VOF_interface&>(ls);
	
	static gridVector2dtype newvelocities (params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
	static grideuler2type prims1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type prims2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type MRPsoln1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type MRPsoln2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static std::vector<std::vector<int>> flag (params.Ny + 2 * params.numGC, std::vector<int>(params.Nx + 2 * params.numGC));
	double ds = 1.0 * std::min(params.dx, params.dy);
	
	int ghostrad = 3;

	if (int(newvelocities.size()) != params.Ny + 2 * params.numGC)
	{
		newvelocities = gridVector2dtype(params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
		prims1 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
		prims2 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
		MRPsoln1 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
		MRPsoln2 = grideuler2type(params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
		flag = std::vector<std::vector<int>>(params.Ny + 2 * params.numGC, std::vector<int>(params.Nx + 2 * params.numGC));
	}
	
	realcells1.imin = params.numGC + params.Ny;
	realcells1.imax = params.numGC;
	realcells1.jmin = params.numGC + params.Nx;
	realcells1.jmax = params.numGC;
	
	realcells2.imin = params.numGC + params.Ny;
	realcells2.imax = params.numGC;
	realcells2.jmin = params.numGC + params.Nx;
	realcells2.jmax = params.numGC;

	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			double z = vof.get_z(i, j);
			newvelocities[i][j] = Eigen::Vector2d::Zero();
			prims1[i][j] = misc::conserved_to_primitives(eosparams.gamma1, eosparams.pinf1, grid1[i][j]);
			prims2[i][j] = misc::conserved_to_primitives(eosparams.gamma2, eosparams.pinf2, grid2[i][j]);
			
			if (z >= 0.5)
			{
				assert(misc::is_physical_state(eosparams.gamma1, eosparams.pinf1, grid1[i][j]));
				
				realcells1.imin = std::min(realcells1.imin, i-2);
				realcells1.imax = std::max(realcells1.imax, i+2);
				realcells1.jmin = std::min(realcells1.jmin, j-2);
				realcells1.jmax = std::max(realcells1.jmax, j+2);
			}
			else
			{
				assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
				
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
			flag[i][j] = -1;
		}
	}
	


	// This loop set flag values and stores real densities/pressures on interface
	
	for (int i=params.numGC-ghostrad; i<params.Ny + params.numGC+ghostrad; i++)
	{
		for (int j=params.numGC-ghostrad; j<params.Nx + params.numGC+ghostrad; j++)
		{
			double z = vof.get_z(i, j);
			
			if (0.0 < z && z < 1.0)
			{
				// Set flag in cells within ghostrad of this mixed cell
				
				for (int xinc=-ghostrad; xinc<=ghostrad; xinc++)
				{
					for (int yinc=-ghostrad; yinc<=ghostrad; yinc++)
					{
						if (ls.get_z(i+yinc, j+xinc) == 0.0)
						{
							flag[i+yinc][j+xinc] = 1;
						}
						else if (ls.get_z(i+yinc, j+xinc) == 1.0)
						{
							flag[i+yinc][j+xinc] = 2;
						}
						else
						{
							flag[i+yinc][j+xinc] = 0;
						}
					}
				}
				
				
				// Save interpolated states for later use in extrapolating entropy

				Eigen::Vector2d normal = vof.get_mixedcell_normal(i,j);
				Eigen::Vector2d interfacelocation = vof.get_interfacelocation(i,j);

				vec4type interfaceprims1;
				vec4type interfaceprims2;
							
				get_interpolated_mixedRPstates_VOF(params, prims1, prims2, normal, interfacelocation, ds, interfaceprims2, interfaceprims1);
				
				MRPsoln1[i][j] = interfaceprims1;
				MRPsoln2[i][j] = interfaceprims2;
				
				
				// Set ghost fluid state in this mixed cell
				
				if (z > 0.5)
				{
					prims2[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma2, eosparams.pinf2, MRPsoln2[i][j](0), MRPsoln2[i][j](3), prims1[i][j](3));
					prims2[i][j](1) = prims1[i][j](1);
					prims2[i][j](2) = prims1[i][j](2);
					prims2[i][j](3) = prims1[i][j](3);
				}
				else
				{
					prims1[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma1, eosparams.pinf1, MRPsoln1[i][j](0), MRPsoln1[i][j](3), prims2[i][j](3));
					prims1[i][j](1) = prims2[i][j](1);
					prims1[i][j](2) = prims2[i][j](2);
					prims1[i][j](3) = prims2[i][j](3);
				}
			}
			
			
			// Store real fluid velocity
			
			if (z >= 0.5)
			{
				newvelocities[i][j] << prims1[i][j](1), prims1[i][j](2);
			}
			else
			{
				newvelocities[i][j] << prims2[i][j](1), prims2[i][j](2);
			}
		}
	}
							
	
	
	// Set ghost fluid state in non-mixed cells around interface			
				
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double z = vof.get_z(i, j);
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);

			if (flag[i][j] == 1)
			{
				// Need to set fluid 1 ghost state

				assert(z == 0.0);
				double suminvdist = 0.0;
				vec4type sumprims = Eigen::Vector4d::Zero();

				for (int xinc=-ghostrad; xinc<=ghostrad; xinc++)
				{
					for (int yinc=-ghostrad; yinc<=ghostrad; yinc++)
					{
						if (flag[i+yinc][j+xinc] == 0)
						{
							double thisz = vof.get_z(i+yinc, j+xinc);
							assert(0.0 < thisz && thisz < 1.0);
							Eigen::Vector2d thisCC = params.cellcentre_coord(i+yinc, j+xinc);

							double thisinvdist = 1.0 / (thisCC - CC).norm();

							suminvdist += thisinvdist;
							sumprims += thisinvdist * MRPsoln1[i+yinc][j+xinc];
						}
					}
				}

				sumprims /= suminvdist;
				
				
				prims1[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma1, eosparams.pinf1, sumprims(0), sumprims(3), prims2[i][j](3));
				prims1[i][j](1) = prims2[i][j](1);
				prims1[i][j](2) = prims2[i][j](2);
				prims1[i][j](3) = prims2[i][j](3);
			}
			else if (flag[i][j] == 2)
			{
				// Need to set fluid 2 ghost state

				assert(z == 1.0);
				double suminvdist = 0.0;
				vec4type sumprims = Eigen::Vector4d::Zero();

				for (int xinc=-ghostrad; xinc<=ghostrad; xinc++)
				{
					for (int yinc=-ghostrad; yinc<=ghostrad; yinc++)
					{
						if (flag[i+yinc][j+xinc] == 0)
						{
							double thisz = vof.get_z(i+yinc, j+xinc);
							assert(0.0 < thisz && thisz < 1.0);
							Eigen::Vector2d thisCC = params.cellcentre_coord(i+yinc, j+xinc);

							double thisinvdist = 1.0 / (thisCC - CC).norm();

							suminvdist += thisinvdist;
							sumprims += thisinvdist * MRPsoln2[i+yinc][j+xinc];
						}
					}
				}

				sumprims /= suminvdist;
				
				prims2[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma2, eosparams.pinf2, sumprims(0), sumprims(3), prims1[i][j](3));
				prims2[i][j](1) = prims1[i][j](1);
				prims2[i][j](2) = prims1[i][j](2);
				prims2[i][j](3) = prims1[i][j](3);
			}
			
			grid1[i][j] = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, prims1[i][j]);
			grid2[i][j] = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, prims2[i][j]);
			assert(misc::is_physical_state(eosparams.gamma1, eosparams.pinf1, grid1[i][j]));
			assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
		}
	}
	
	
	apply_BCs_euler(params, grid1);
	apply_BCs_euler(params, grid2);
	vfield->store_velocity_field(newvelocities, t);
}
