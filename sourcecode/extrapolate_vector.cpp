#include "extrapolate_vector.hpp"
#include "euler_bc.hpp"


void extrapolate_vector
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	grideuler2type& grid1, 
	grideuler2type& grid2
)
{	
	static grideuler2type updated_grid1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type updated_grid2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	vec4type del_vec_x;
	vec4type del_vec_y;
	double dt = 0.5 * std::min(params.dx, params.dy);


	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			updated_grid1[i][j] = grid1[i][j];
			updated_grid2[i][j] = grid2[i][j];
			assert(grid1[i][j](0) >= 0.0);
			assert(grid1[i][j](3) >= 0.0);
			assert(grid2[i][j](0) >= 0.0);
			assert(grid2[i][j](3) >= 0.0);
		}
	}
	
	
	for (int num=0; num<N; num++)
	{
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				grid1[i][j] = updated_grid1[i][j];
				grid2[i][j] = updated_grid2[i][j];
			}
		}
		
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=params.numGC; j<params.Nx + params.numGC; j++)
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
				
				if (lsval <= 0.0)
				{
					// Fluid 2 is ghost here - information flows in direction of -normal
					
					if (normal(0) > 0.0)
					{
						del_vec_x = grid2[i][j+1] - grid2[i][j];
					}
					else
					{
						del_vec_x = grid2[i][j] - grid2[i][j-1];
					}
					
					if (normal(1) > 0.0)
					{
						del_vec_y = grid2[i+1][j] - grid2[i][j];
					}
					else
					{
						del_vec_y = grid2[i][j] - grid2[i-1][j];
					}
					
					updated_grid2[i][j] = grid2[i][j] + (dt / params.dx) * normal(0) * del_vec_x + (dt / params.dy) * normal(1) * del_vec_y;
				}
				else
				{
					// Fluid 1 is ghost here - information flows in direction of +normal
					
					if (normal(0) > 0.0)
					{
						del_vec_x = grid1[i][j] - grid1[i][j-1];
					}
					else
					{
						del_vec_x = grid1[i][j+1] - grid1[i][j];
					}
					
					if (normal(1) > 0.0)
					{
						del_vec_y = grid1[i][j] - grid1[i-1][j];
					}
					else
					{
						del_vec_y = grid1[i+1][j] - grid1[i][j];
					}
					
					updated_grid1[i][j] = grid1[i][j] - (dt / params.dx) * normal(0) * del_vec_x - (dt / params.dy) * normal(1) * del_vec_y;
				}
				
				assert(updated_grid1[i][j](0) >= 0.0);
				assert(updated_grid1[i][j](3) >= 0.0);
				assert(updated_grid2[i][j](0) >= 0.0);
				assert(updated_grid2[i][j](3) >= 0.0);
			}
		}
		
		apply_BCs_euler(params, updated_grid1);
		apply_BCs_euler(params, updated_grid2);
	}
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			grid1[i][j] = updated_grid1[i][j];
			grid2[i][j] = updated_grid2[i][j];
		}
	}
}








void extrapolate_vector_mgfm
(
	const sim_info& params,
	const GFM_ITM_interface& ls,
	const int N,
	grideuler2type& grid1, 
	grideuler2type& grid2
)
{	
	static grideuler2type updated_grid1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	static grideuler2type updated_grid2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
	vec4type del_vec_x;
	vec4type del_vec_y;
	double dt = 0.5 * std::min(params.dx, params.dy);


	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			updated_grid1[i][j] = grid1[i][j];
			updated_grid2[i][j] = grid2[i][j];
			assert(grid1[i][j](0) >= 0.0);
			assert(grid1[i][j](3) >= 0.0);
			assert(grid2[i][j](0) >= 0.0);
			assert(grid2[i][j](3) >= 0.0);
		}
	}
	
	
	for (int num=0; num<N; num++)
	{
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				grid1[i][j] = updated_grid1[i][j];
				grid2[i][j] = updated_grid2[i][j];
			}
		}
		
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=params.numGC; j<params.Nx + params.numGC; j++)
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
				
				if (!ls.is_interfacial_cell(i, j))
				{
				
					if (lsval <= 0.0)
					{
						// Fluid 2 is ghost here - information flows in direction of -normal
						
						if (normal(0) > 0.0)
						{
							del_vec_x = grid2[i][j+1] - grid2[i][j];
						}
						else
						{
							del_vec_x = grid2[i][j] - grid2[i][j-1];
						}
						
						if (normal(1) > 0.0)
						{
							del_vec_y = grid2[i+1][j] - grid2[i][j];
						}
						else
						{
							del_vec_y = grid2[i][j] - grid2[i-1][j];
						}
						
						updated_grid2[i][j] = grid2[i][j] + (dt / params.dx) * normal(0) * del_vec_x + (dt / params.dy) * normal(1) * del_vec_y;
					}
					else
					{
						// Fluid 1 is ghost here - information flows in direction of +normal
						
						if (normal(0) > 0.0)
						{
							del_vec_x = grid1[i][j] - grid1[i][j-1];
						}
						else
						{
							del_vec_x = grid1[i][j+1] - grid1[i][j];
						}
						
						if (normal(1) > 0.0)
						{
							del_vec_y = grid1[i][j] - grid1[i-1][j];
						}
						else
						{
							del_vec_y = grid1[i+1][j] - grid1[i][j];
						}
						
						updated_grid1[i][j] = grid1[i][j] - (dt / params.dx) * normal(0) * del_vec_x - (dt / params.dy) * normal(1) * del_vec_y;
					}
					
					assert(updated_grid1[i][j](0) >= 0.0);
					assert(updated_grid1[i][j](3) >= 0.0);
					assert(updated_grid2[i][j](0) >= 0.0);
					assert(updated_grid2[i][j](3) >= 0.0);
				}
			}
		}
		
		apply_BCs_euler(params, updated_grid1);
		apply_BCs_euler(params, updated_grid2);
	}
	
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			grid1[i][j] = updated_grid1[i][j];
			grid2[i][j] = updated_grid2[i][j];
		}
	}
}
