#include "GFM_real.hpp"
#include "extrapolate_extension_vfield.hpp"
#include "extrapolate_vector.hpp"
#include "euler_misc.hpp"
#include "mixed_RS_exact.hpp"


void GFM_real :: set_ghost_states
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
			double lsval = ls.get_sdf(i, j);
			newvelocities[i][j] = Eigen::Vector2d::Zero();
			prims1[i][j] = misc::conserved_to_primitives(eosparams.gamma1, eosparams.pinf1, grid1[i][j]);
			prims2[i][j] = misc::conserved_to_primitives(eosparams.gamma2, eosparams.pinf2, grid2[i][j]);
			
			if (lsval <= 0.0)
			{
				assert(misc::is_physical_state(eosparams.gamma1, eosparams.pinf1, grid1[i][j]));
			}
			else
			{
				assert(misc::is_physical_state(eosparams.gamma2, eosparams.pinf2, grid2[i][j]));
			}
		}
	}
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			if (ls.is_interfacial_cell(i,j))
			{
				double lsval = ls.get_sdf(i, j);
				Eigen::Vector2d normal = ls.get_normal(params.cellcentre_coord(i, j));
					
				if (normal.norm() < 1e-12)
				{
					normal(0) = 1.0;
					normal(1) = 0.0;
				}
				else
				{
					normal /= normal.norm();
				}
				
				vec4type interfaceprims1;
				vec4type interfaceprims2;
				
				
				// Search over all neighbouring cells for one which minimises the difference in angles between normals
				
				double mindiff = 1e100;
				int i_opp = -1000;
				int j_opp = -1000;
				
				for (int xinc=-1; xinc<=1; xinc++)
				{
					for (int yinc=-1; yinc<=1; yinc++)
					{
						double otherlsval = ls.get_sdf(i+yinc, j+xinc);
						
						if (lsval * otherlsval <= 0.0)
						{
							Eigen::Vector2d othernormal = ls.get_normal(params.cellcentre_coord(i+yinc, j+xinc));
							
							if (othernormal.norm() < 1e-12)
							{
								othernormal(0) = 1.0;
								othernormal(1) = 0.0;
							}
							else
							{
								othernormal /= othernormal.norm();
							}
							
							double dotprod = normal.dot(othernormal);
							dotprod = std::max(dotprod, -1.0);
							dotprod = std::min(dotprod, 1.0);
							double anglediff = acos(dotprod);
							
							if (anglediff < mindiff)
							{
								i_opp = i + yinc;
								j_opp = j + xinc;
								mindiff = anglediff;
							}
						}	
					}
				}
				
				assert(i_opp != -1000);
				assert(j_opp != -1000);
				
				if (lsval <= 0.0)
				{
					interfaceprims1 = prims1[i][j];
					interfaceprims2 = prims2[i_opp][j_opp];
				}
				else
				{
					interfaceprims2 = prims2[i][j];
					interfaceprims1 = prims1[i_opp][j_opp];
				}
								
				
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
				
				if (p_star != p_star)
				{
					std::cout << "Exact RS returned NaN" << std::endl;
				}
				
				
				// Compute star-state velocities in Cartesian frame
				
				Eigen::Vector2d u1_star = u_star * normal + u1_tang;
				Eigen::Vector2d u2_star = u_star * normal + u2_tang;
				
				
				// Set the real fluid primitives here to the star-state primitives from MRP
				
				if (lsval <= 0.0)
				{
					prims1[i][j](0) = rho_star_L;
					prims1[i][j](1) = u1_star(0);
					prims1[i][j](2) = u1_star(1);
					prims1[i][j](3) = p_star;
				}
				else
				{
					prims2[i][j](0) = rho_star_R;
					prims2[i][j](1) = u2_star(0);
					prims2[i][j](2) = u2_star(1);
					prims2[i][j](3) = p_star;
				}
				
				
				// Set the interface extension vfield in this cell
				
				newvelocities[i][j] = u_star * normal;
			}
		}
	}
	
	extrapolate_vector(params, ls, 6, prims1, prims2);
	
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
}
	
	
