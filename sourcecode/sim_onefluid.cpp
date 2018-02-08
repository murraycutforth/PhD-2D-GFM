#include "sim_onefluid.hpp"
#include "euler_misc.hpp"
#include "BBrange.hpp"
#include "FS_godunov.hpp"
#include "FS_MUSCL.hpp"
#include "pure_RS_exact.hpp"
#include "euler_bc.hpp"
#include <cassert>
#include <fstream>


void sim_onefluid :: run_sim (GFM_settingsfile SF)
{
	sim_info params;
	binarySGparams eosparams;
	set_sim_parameters(SF, params, eosparams);
	
	
	std::shared_ptr<pure_RS_base> RS;
	std::shared_ptr<FS_base> FS;
	set_sim_methods(SF, RS, FS);
	
	
	grideuler2type grid (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC, vec4type(4)));
	grideuler2type future_grid (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC, vec4type(4)));
	BBrange realcells (params.numGC, params.Ny + params.numGC, params.numGC, params.Nx + params.numGC);
	set_sim_ICs(SF, params, eosparams, grid);
	

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;
	
	output(numsteps, t, params, eosparams, grid, SF.basename);

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations with CFL = " << SF.CFL << "." << std::endl;

	double tstart = omp_get_wtime();

	while (t < params.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		dt = compute_dt(CFL, params, grid, eosparams, params.T, t);

		FS->pure_fluid_update(params, eosparams.gamma1, eosparams.pinf1, dt, realcells, grid, future_grid);

		apply_BCs_euler (params, grid);
		
		numsteps++;
		t += dt;
		if (numsteps % params.output_freq == 0) output(numsteps, t, params, eosparams, grid, SF.basename);
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}
	
	double tend = omp_get_wtime();

	
	
	output(numsteps, t, params, eosparams, grid, SF.basename);

	std::cout << "[" << SF.basename << "] Simulation complete." << std::endl;
	std::cout << "[" << SF.basename << "] Simulation took " << tend - tstart << "s." << std::endl;
}



double sim_onefluid :: compute_dt (double CFL, const sim_info& params, const grideuler2type& grid, const binarySGparams& eosparams, double T, double t)
{
	double maxu = 0.0, maxv = 0.0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0);
			double e = misc::specific_ie(grid[i][j]);
			double p = eos::pressure(eosparams.gamma1, eosparams.pinf1, e, rho);
			double soundspeed = eos::soundspeed(eosparams.gamma1, eosparams.pinf1, p, rho);
			
			maxu = std::max(maxu, fabs(grid[i][j](1))/rho + soundspeed);
			maxv = std::max(maxv, fabs(grid[i][j](2))/rho + soundspeed);
		}
	}
			
	double dt = CFL * std::min(params.dx/maxu, params.dy/maxv);

	if (t + dt > T) dt = T - t;

	return dt;
}




void sim_onefluid :: set_sim_parameters (const GFM_settingsfile& SF, sim_info& params, binarySGparams& eosparams)
{
	assert(SF.sim_type == "onefluid");
	params.numGC = 2;
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.output_freq = 1e9;
	eosparams.gamma1 = 1.4;
	eosparams.gamma2 = 0.0;
	eosparams.pinf1 = 0.0;
	eosparams.pinf2 = 0.0;
	params.x0 = 0.0;
	params.y0 = 0.0;
	params.dx = 1.0/params.Nx;
	params.dy = 1.0/params.Ny;
	params.BC_L = "transmissive";
	params.BC_T = "transmissive";
	params.BC_R = "transmissive";
	params.BC_B = "transmissive";
	
	if (SF.test_case == "TTC1_x" || SF.test_case == "TTC1_y")
	{
		params.T = 0.25;
	}
	else if (SF.test_case == "TTC2_x" || SF.test_case == "TTC2_y")
	{
		params.T = 0.15;
	}
	else if (SF.test_case == "TTC3_x" || SF.test_case == "TTC3_y")
	{
		params.T = 0.012;
	}
	else if (SF.test_case == "TTC4_x" || SF.test_case == "TTC4_y")
	{
		params.T = 0.035;
	}
	else if (SF.test_case == "TTC5_x" || SF.test_case == "TTC5_y")
	{
		params.T = 0.035;
	}
	else if (SF.test_case == "circularexplosion")
	{
		params.dx = 2.0/params.Nx;
		params.dy = 2.0/params.Ny;
		params.T = 0.25;
	}
	else
	{
		assert(!"Invalid test_case choice.");
	}
}


void sim_onefluid :: set_sim_methods (const GFM_settingsfile& SF, std::shared_ptr<pure_RS_base>& RS, std::shared_ptr<FS_base>& FS)
{
	RS = std::make_shared<pure_RS_exact>();
	
	if (SF.FS == "Godunov")
	{
		FS = std::make_shared<FS_godunov>(RS);
	}
	else if (SF.FS == "MUSCL")
	{
		FS = std::make_shared<FS_MUSCL>(RS);
	}
	else
	{
		assert(!"Invalid method choice");
	}
}
	

	
void sim_onefluid :: output (int numsteps, double t, const sim_info& params, const binarySGparams& eosparams, const grideuler2type& grid, std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename + "-" + std::to_string(numsteps) + "-prims.dat");
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);
			
			double rho = grid[i][j](0);
			double u = grid[i][j](1) / rho;
			double v = grid[i][j](2) / rho;
			double e = misc::specific_ie(grid[i][j]);
			double p = eos::pressure(eosparams.gamma1, eosparams.pinf1, e, rho);
			
			outfile << CC(0) << " " << CC(1) << " " << rho << " " << u << " " << v << " " << e << " " << p << std::endl;
		}
		outfile << std::endl;
	}
	
	outfile.close();	
}



void sim_onefluid :: set_sim_ICs (const GFM_settingsfile& SF, const sim_info& params, const binarySGparams& eosparams, grideuler2type& grid)
{
	if (SF.test_case == "TTC1_x" || SF.test_case == "TTC1_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1.0;
		
		Rprims(0) = 0.125;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 0.1;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		if (SF.test_case == "TTC1_x")
		{
			set_planar_IC(Lstate, Rstate, 1.0, 0.0, 0.5, grid, params);
		}
		else
		{
			set_planar_IC(Lstate, Rstate, 0.0, 1.0, 0.5, grid, params);
		}
	}
	else if (SF.test_case == "TTC2_x")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = -2.0;
		Lprims(2) = 0.0;
		Lprims(3) = 0.4;
		
		Rprims(0) = 1.0;
		Rprims(1) = 2.0;
		Rprims(2) = 0.0;
		Rprims(3) = 0.4;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		set_planar_IC(Lstate, Rstate, 1.0, 0.0, 0.5, grid, params);
	}
	else if (SF.test_case == "TTC2_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = -2.0;
		Lprims(3) = 0.4;
		
		Rprims(0) = 1.0;
		Rprims(1) = 0.0;
		Rprims(2) = 2.0;
		Rprims(3) = 0.4;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		set_planar_IC(Lstate, Rstate, 0.0, 1.0, 0.5, grid, params);
	}
	else if (SF.test_case == "TTC3_x" || SF.test_case == "TTC3_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1000.0;
		
		Rprims(0) = 1.0;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 0.01;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		if (SF.test_case == "TTC3_x")
		{
			set_planar_IC(Lstate, Rstate, 1.0, 0.0, 0.5, grid, params);
		}
		else
		{
			set_planar_IC(Lstate, Rstate, 0.0, 1.0, 0.5, grid, params);
		}
	}
	else if (SF.test_case == "TTC4_x" || SF.test_case == "TTC4_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 0.1;
		
		Rprims(0) = 1.0;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 100.0;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		if (SF.test_case == "TTC4_x")
		{
			set_planar_IC(Lstate, Rstate, 1.0, 0.0, 0.5, grid, params);
		}
		else
		{
			set_planar_IC(Lstate, Rstate, 0.0, 1.0, 0.5, grid, params);
		}
	}
	else if (SF.test_case == "TTC5_x")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 5.99924;
		Lprims(1) = 19.5975;
		Lprims(2) = 0.0;
		Lprims(3) = 460.894;
		
		Rprims(0) = 5.99242;
		Rprims(1) = -6.19633;
		Rprims(2) = 0.0;
		Rprims(3) = 46.0950;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		set_planar_IC(Lstate, Rstate, 1.0, 0.0, 0.5, grid, params);
	}
	else if (SF.test_case == "TTC5_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 5.99924;
		Lprims(1) = 0.0;
		Lprims(2) = 19.5975;
		Lprims(3) = 460.894;
		
		Rprims(0) = 5.99242;
		Rprims(1) = 0.0;
		Rprims(2) = -6.19633;
		Rprims(3) = 46.0950;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		set_planar_IC(Lstate, Rstate, 0.0, 1.0, 0.5, grid, params);
	}
	else if (SF.test_case == "circularexplosion")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1.0;
		
		Rprims(0) = 0.125;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 0.1;
		
		int N = 10;
		Eigen::Vector2d centre;
		centre << 1.0, 1.0;
		double R = 0.4;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{				
				int totalnumsamples = N*N;
				int numinside = 0;
				double delx = params.dx/N;
				double dely = params.dy/N;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				Eigen::Vector2d samplepos;
				
				for (int a=0; a<N; a++)
				{
					for (int b=0; b<N; b++)
					{
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos -= centre;
						
						if (samplepos.norm() <= R) numinside++;
					}
				}
	
				double frac = double(numinside)/totalnumsamples;
				
				vec4type W = frac * Lprims + (1.0 - frac) * Rprims;
				
				grid[i][j] = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, W);
			}
		}

	}
	else
	{
		assert(!"Invalid test_case choice.");
	}
	
	apply_BCs_euler (params, grid);
}



void sim_onefluid :: set_planar_IC (const vec4type& U_under, const vec4type& U_over, const double a, const double b, const double c, grideuler2type& grid, const sim_info& params)
{
	/*
	 * Set every cell where a*x + b*y <= c to the state U_under,
	 * otherwise U_over. Useful for planar initial conditions.
	 */
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			
			if (cc(0)*a + cc(1)*b <= c)
			{
				grid[i][j] = U_under;
			}
			else
			{
				grid[i][j] = U_over;
			}
		}
	}
}
