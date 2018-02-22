#include "extrapolate_extension_vfield.hpp"
#include "boundary_conditions.hpp"

void extrapolate_extension_vfield
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	gridVector2dtype& velocities
)
{
	static gridVector2dtype updated_velocities (params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
	Eigen::Vector2d del_v_x;
	Eigen::Vector2d del_v_y;
	double dt = 0.5 * std::min(params.dx, params.dy);
	
	if (int(updated_velocities.size()) != params.Ny + 2 * params.numGC)
	{
		updated_velocities = gridVector2dtype(params.Ny + 2 * params.numGC, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>(params.Nx + 2 * params.numGC));
	}
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			updated_velocities[i][j] = velocities[i][j];
		}
	}
	
	
	for (int num=0; num<N; num++)
	{
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				velocities[i][j] = updated_velocities[i][j];
			}
		}
		
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=params.numGC; j<params.Nx + params.numGC; j++)
			{
				if (!ls.is_interfacial_cell(i,j))
				{
					double lsval = ls.get_sdf(i, j);
					Eigen::Vector2d normal = ls.get_normal(params.cellcentre_coord(i, j));
					

					if (normal.norm() > 1e-12)
					{
						
						normal /= normal.norm();
						
						if (lsval <= 0.0)
						{
							// Fluid 2 is ghost here - information flows in direction of -normal
							
							if (normal(0) > 0.0)
							{
								del_v_x = velocities[i][j+1] - velocities[i][j];
							}
							else
							{
								del_v_x = velocities[i][j] - velocities[i][j-1];
							}
							
							if (normal(1) > 0.0)
							{
								del_v_y = velocities[i+1][j] - velocities[i][j];
							}
							else
							{
								del_v_y = velocities[i][j] - velocities[i-1][j];
							}
							
							updated_velocities[i][j] = velocities[i][j] + (dt / params.dx) * normal(0) * del_v_x + (dt / params.dy) * normal(1) * del_v_y;
						}
						else
						{
							// Fluid 1 is ghost here - information flows in direction of +normal
							
							if (normal(0) > 0.0)
							{
								del_v_x = velocities[i][j] - velocities[i][j-1];
							}
							else
							{
								del_v_x = velocities[i][j+1] - velocities[i][j];
							}
							
							if (normal(1) > 0.0)
							{
								del_v_y = velocities[i][j] - velocities[i-1][j];
							}
							else
							{
								del_v_y = velocities[i+1][j] - velocities[i][j];
							}
							
							updated_velocities[i][j] = velocities[i][j] - (dt / params.dx) * normal(0) * del_v_x - (dt / params.dy) * normal(1) * del_v_y;
						}

					}
				}
			}
		}
		
		apply_gridVector2dtype_BCs(params, updated_velocities);
	}
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			velocities[i][j] = updated_velocities[i][j];
		}
	}
}
