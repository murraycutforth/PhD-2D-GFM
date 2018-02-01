#include "sim_onefluid.hpp"
#include "euler_misc.hpp"
#include "BBrange.hpp"
#include "FS_godunov.hpp"
#include "pure_RS_exact.hpp"
#include "euler_bc.hpp"


void sim_onefluid :: run_sim (GFM_settingsfile SF)
{
	sim_info params;
	binarySGparams eosparams;
	set_sim_parameters(SF, params, eosparams);
	
	
	std::shared_ptr<pure_RS_base> RS;
	std::shared_ptr<flow_solver_base> FS;
	set_sim_methods(SF, RS, FS);
	
	
	gridtype grid (params.Ny + 2 * params.numGC, rowtype(params.Nx + 2 * params.numGC, vectype(4)));
	gridtype future_grid (params.Ny + 2 * params.numGC, rowtype(params.Nx + 2 * params.numGC, vectype(4)));
	BBrange realcells (params.numGC, params.Ny + params.numGC, params.numGC, params.Nx + params.numGC);
	set_sim_ICs(SF, params, eosparams, grid);
	

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;
	
	output(numsteps, t, params, eosparams, grid, SF.basename);

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations.." << std::endl;

	double tstart = omp_get_wtime();

	while (t < params.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		dt = compute_dt(CFL, params, grid, eosparams, params.T, t);

		FS->pure_fluid_update(params, eosparams, dt, realcells, grid, future_grid);

		grid.swap(future_grid);
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



double sim_onefluid :: compute_dt (double CFL, const sim_info& params, const gridtype& grid, const binarySGparams& eosparams, double T, double t)
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
	
	if (SF.test_case == "TTC1_x" || SF.test_case == "TTC1_y")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 0.0;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
		params.Nx = SF.Nx;
		params.Ny = SF.Ny;
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.25;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
		params.output_freq = 1e9;
	}
	else
	{
		assert(!"Invalid test_case choice.");
	}
}


void sim_onefluid :: set_sim_methods (const GFM_settingsfile& SF, std::shared_ptr<pure_RS_base>& RS, std::shared_ptr<flow_solver_base>& FS)
{
	RS = std::make_shared<pure_RS_exact>();
	
	if (SF.method == "Godunov")
	{
		FS = std::make_shared<FS_godunov>(RS);
	}
	else if (SF.method == "MUSCL")
	{
		//FS = std::make_shared<>();
	}
	else
	{
		assert(!"Invalid method choice");
	}
}
	
	
void sim_onefluid :: set_sim_ICs (const GFM_settingsfile& SF, const sim_info& params, const binarySGparams& eosparams, gridtype& grid)
{
	// TODO
}

	
void sim_onefluid :: output (int numsteps, double t, const GFM_settingsfile& SF, const sim_info& params, const binarySGparams& eosparams, const gridtype& grid, std::string filename)
{
	// TODO
}
