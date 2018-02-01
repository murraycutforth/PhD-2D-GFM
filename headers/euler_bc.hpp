#ifndef EULER_BC_H
#define EULER_BC_H


void apply_BCs_euler (const sim_info& params, gridtype& grid)
{
	if (	params.BC_T == "periodic"
		&& params.BC_L == "periodic"
		&& params.BC_B == "periodic"
		&& params.BC_R == "periodic" )
	{
	
		for (int i=0; i<params.Nx + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Ny + 2 * params.numGC; j++)
			{
				grid[i][j] = grid[(i-params.numGC+params.Ny) % params.Ny + params.numGC][(j-params.numGC+params.Nx) % params.Nx + params.numGC];
			}
		}
	}
	else 
	{
		if (params.BC_L == "transmissive" || params.BC_L == "reflective")
		{
			for (int i=params.numGC; i<params.Ny + params.numGC; i++)
			{
				for (int j=0; j<params.numGC; j++)
				{
					grid[i][j] = grid[i][2 * params.numGC - j - 1];
					
					if (params.BC_L == "reflective")
					{
						grid[i][j](1) = - grid[i][j](1);
					}
				}
			}
		}
		else
		{
			assert(!"[boundary_conditions] Invalid BC_L");
		}
		
		if (params.BC_R == "transmissive" || params.BC_R == "reflective")
		{
			for (int i=params.numGC; i<params.Ny + params.numGC; i++)
			{
				for (int j=0; j<params.numGC; j++)
				{
					grid[i][j + params.Nx + params.numGC] = grid[i][params.Nx + params.numGC - 1 - j];
					
					if (params.BC_L == "reflective")
					{
						grid[i][j + params.Nx + params.numGC](1) = - grid[i][j + params.Nx + params.numGC](1);
					}
				}
			}
		}
		else
		{
			assert(!"[boundary_conditions] Invalid BC_R");
		}
		
		if (params.BC_B == "transmissive" || params.BC_B == "reflective")
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				for (int i=0; i<params.numGC; i++)
				{
					grid[i][j] = grid[2 * params.numGC - i - 1][j];
					
					if (params.BC_B == "reflective")
					{
						grid[i][j](2) = - grid[i][j](2);
					}
				}
			}
		}
		else
		{
			assert(!"[interface-tracking-methods] Invalid BC_B");
		}
		
		if (params.BC_T == "transmissive" || params.BC_T == "reflective")
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				for (int i=0; i<params.numGC; i++)
				{
					grid[i + params.Ny + params.numGC][j] = grid[params.Ny + params.numGC - 1 - i][j];
					
					if (params.BC_T == "reflective")
					{
						grid[i + params.Ny + params.numGC][j](2) = - grid[i + params.Ny + params.numGC][j](2);
					}
				}
			}
		}
		else
		{
			assert(!"[boundary_conditions] Invalid BC_T");
		}
	}
}


#endif
