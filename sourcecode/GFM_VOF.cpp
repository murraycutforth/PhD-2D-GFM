#include "GFM_VOF.hpp"
#include "euler_misc.hpp"
#include "mixed_RS_exact.hpp"
#include "euler_bc.hpp"
#include "GFM_VOF_interface.hpp"


void GFM_VOF :: set_ghost_states
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
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			flag[i][j] = -1;
		}
	}

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
			}
			else
			{
				assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
			}
		}
	}


	// This loop sets ghost states in mixed cells and stores MRP solution
	
	for (int i=params.numGC-ghostrad; i<params.Ny + params.numGC+ghostrad; i++)
	{
		for (int j=params.numGC-ghostrad; j<params.Nx + params.numGC+ghostrad; j++)
		{
			double z = vof.get_z(i, j);
			
			if (0.0 < z && z < 1.0)
			{
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
			}

			if (0.0 < z && z < 1.0)
			{
				Eigen::Vector2d normal = vof.get_mixedcell_normal(i,j);
				Eigen::Vector2d interfacelocation = vof.get_interfacelocation(i,j);

				vec4type interfaceprims1;
				vec4type interfaceprims2;
							
				get_interpolated_mixedRPstates_VOF(params, prims1, prims2, normal, interfacelocation, ds, interfaceprims2, interfaceprims1);
								
				
				// Fluid velocities in Cartesian frame
				
				Eigen::Vector2d u1, u2;
				u1 << interfaceprims1(1), interfaceprims1(2);
				u2 << interfaceprims2(1), interfaceprims2(2);
				

				// Fluid velocities normal to interface
				
				Eigen::Vector2d u1_normal = u1.dot(normal) * normal;
				Eigen::Vector2d u2_normal = u2.dot(normal) * normal;
				
	
				// Fluid velocities parallel to interface
				
				Eigen::Vector2d u1_tang = u1 - u1_normal;
				Eigen::Vector2d u2_tang = u2 - u2_normal;
				

								
				// Solve mixed Riemann problem in frame tangential to interface
				
				double p_star, u_star, rho_star_L, rho_star_R;
				vec4type Lstate;
				vec4type Rstate;
				vec4type Lprims;
				vec4type Rprims;
				
				Lprims(0) = interfaceprims1(0);
				Lprims(1) = u1.dot(normal);
				Lprims(2) = 0.0;
				Lprims(3) = std::max(1e-6, interfaceprims1(3));
				
				Rprims(0) = interfaceprims2(0);
				Rprims(1) = u2.dot(normal);
				Rprims(2) = 0.0;
				Rprims(3) = std::max(1e-6, interfaceprims2(3));
				
				Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
				Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);
				
				mixed_RS->solve_mixed_RP(eosparams.gamma1, eosparams.pinf1, eosparams.gamma2, eosparams.pinf2, Lstate, Rstate, p_star, u_star, rho_star_L, rho_star_R);
				
				
				// Compute star-state velocities in Cartesian frame
				
				Eigen::Vector2d u1_star = u_star * normal + u1_tang;
				Eigen::Vector2d u2_star = u_star * normal + u2_tang;
				
				
				//~ // OR - apply no-slip interface BC
				
				//~ Eigen::Vector2d u_tang_avg = 0.5 * (u1_tang + u2_tang);
				//~ Eigen::Vector2d u1_star = u_star * normal + u2_tang;
				//~ Eigen::Vector2d u2_star = u_star * normal + u2_tang;
				


				// Store RP solution in this mixed cell

				MRPsoln1[i][j](0) = rho_star_L;
				MRPsoln1[i][j](1) = u1_star(0);
				MRPsoln1[i][j](2) = u1_star(1);
				MRPsoln1[i][j](3) = p_star;

				MRPsoln2[i][j](0) = rho_star_R;
				MRPsoln2[i][j](1) = u2_star(0);
				MRPsoln2[i][j](2) = u2_star(1);
				MRPsoln2[i][j](3) = p_star;


				// Set real/ghost states in this mixed cell

				if (z > 0.5)
				{
					prims1[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma1, eosparams.pinf1, rho_star_L, p_star, prims1[i][j](3));
					
					prims2[i][j](0) = rho_star_R;
					prims2[i][j](1) = u2_star(0);
					prims2[i][j](2) = u2_star(1);
					prims2[i][j](3) = p_star;
				}
				else
				{
					prims2[i][j](0) = eos::isentropic_extrapolation(eosparams.gamma2, eosparams.pinf2, rho_star_R, p_star, prims2[i][j](3));
					
					prims1[i][j](0) = rho_star_L;
					prims1[i][j](1) = u1_star(0);
					prims1[i][j](2) = u1_star(1);
					prims1[i][j](3) = p_star;
				}
			}
		}
	}
	

	// This loop set ghost states around interface
	
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
				prims1[i][j] = sumprims;
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
				prims2[i][j] = sumprims;
			}
		}
	}


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
	
	
	apply_BCs_euler(params, grid1);
	apply_BCs_euler(params, grid2);
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			if (ls.get_z(i,j) >= 0.5)
			{
				newvelocities[i][j] << prims1[i][j](1), prims1[i][j](2);
			}
			else
			{
				newvelocities[i][j] << prims2[i][j](1), prims2[i][j](2);
			}
		}
	}
	
	
	vfield->store_velocity_field(newvelocities, t);
	
	
}
