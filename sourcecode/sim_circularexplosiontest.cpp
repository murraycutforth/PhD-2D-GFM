#include "sim_circularexplosiontest.hpp"
#include "euler_misc.hpp"
#include "pure_RS_exact.hpp"
#include "FS_MUSCL.hpp"
#include "euler_bc.hpp"
#include <fstream> 



void sim_circularexplosiontest :: run_sim (GFM_settingsfile SF)
{
	// Set up single material problem with 4000x4000 resolution
	// Advance to end of time steps
	
	sim_info params1;
	binarySGparams eosparams1;
	
	params1.numGC = 2;
	params1.Nx = 4000;
	params1.Ny = 4000;
	params1.output_freq = 1e9;
	eosparams1.gamma1 = 1.4;
	eosparams1.gamma2 = 0.0;
	eosparams1.pinf1 = 0.0;
	eosparams1.pinf2 = 0.0;
	params1.x0 = 0.0;
	params1.y0 = 0.0;
	params1.BC_L = "transmissive";
	params1.BC_T = "transmissive";
	params1.BC_R = "transmissive";
	params1.BC_B = "transmissive";
	params1.dx = 2.0/params1.Nx;
	params1.dy = 2.0/params1.Ny;
	params1.T = 0.25;
	
	std::shared_ptr<pure_RS_base> RS = std::make_shared<pure_RS_exact>();
	std::shared_ptr<FS_base> FS = std::make_shared<FS_MUSCL>(RS);
	
	grideuler2type grid (params1.Ny + 2 * params1.numGC, roweuler2type(params1.Nx + 2 * params1.numGC, vec4type(4)));
	grideuler2type future_grid (params1.Ny + 2 * params1.numGC, roweuler2type(params1.Nx + 2 * params1.numGC, vec4type(4)));
	BBrange realcells (params1.numGC, params1.Ny + params1.numGC, params1.numGC, params1.Nx + params1.numGC);
	
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
	
	for (int i=0; i<params1.Ny + 2 * params1.numGC; i++)
	{
		for (int j=0; j<params1.Nx + 2 * params1.numGC; j++)
		{				
			int totalnumsamples = N*N;
			int numinside = 0;
			double delx = params1.dx/N;
			double dely = params1.dy/N;
			
			Eigen::Vector2d cc = params1.cellcentre_coord(i, j);
			Eigen::Vector2d BL;
			BL(0) = cc(0) - 0.5 * params1.dx;
			BL(1) = cc(1) - 0.5 * params1.dy;
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
			
			grid[i][j] = misc::primitives_to_conserved(eosparams1.gamma1, eosparams1.pinf1, W);
		}
	}
	
	apply_BCs_euler (params1, grid);
	
	std::cout << "Initialisation of single fluid complete. Beginning time iterations with CFL = " << SF.CFL << "." << std::endl;
	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;

	while (t < params1.T)
	{
		CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
		dt = compute_onefluid_dt(CFL, params1, grid, eosparams1, params1.T, t);

		FS->pure_fluid_update(params1, eosparams1.gamma1, eosparams1.pinf1, dt, realcells, grid, future_grid);

		apply_BCs_euler (params1, grid);
		
		numsteps++;
		t += dt;
		
		std::cout << "[CEtest_onefluid] Time step " << numsteps << " complete. t = " << t << std::endl;
	}
		
	
	
	
	
	// For each GFM
		// For each resolution
			// Setup and run GFM to end
			// Compute L1, L2, and Linf errors compared to reference solution
	
	std::vector<std::string> GFMs = {"OGFM", "MGFM", "realGFM", "riemannGFM"};
	std::vector<int> res = {50, 100, 200, 400};
	
	for (std::string thisGFM : GFMs)
	{
		std::ofstream outfile;
		outfile.open("./output/CEtest-" + thisGFM + ".dat");
		
		SF.test_case = "circularexplosion";
		SF.FS = "MUSCL";
		SF.GFM = thisGFM;
		
		for (int thisres : res)
		{
			SF.Nx = thisres;
			SF.Ny = thisres;
			
			sim_info params;
			binarySGparams eosparams;
			
			params.numGC = 5;
			params.Nx = SF.Nx;
			params.Ny = SF.Ny;
			params.output_freq = 1e9;
			eosparams.gamma1 = 1.4;
			eosparams.gamma2 = 1.4;
			eosparams.pinf1 = 0.0;
			eosparams.pinf2 = 0.0;
			params.x0 = 0.0;
			params.y0 = 0.0;
			params.BC_L = "transmissive";
			params.BC_T = "transmissive";
			params.BC_R = "transmissive";
			params.BC_B = "transmissive";
			params.dx = 2.0/params.Nx;
			params.dy = 2.0/params.Ny;
			params.T = 0.25;
			
			std::shared_ptr<GFM_base> GFM;
			std::shared_ptr<GFM_ITM_interface> ls;
			set_sim_methods(SF, params, FS, GFM, ls);
			
			grideuler2type grid1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
			grideuler2type future_grid1 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
			grideuler2type grid2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
			grideuler2type future_grid2 (params.Ny + 2 * params.numGC, roweuler2type(params.Nx + 2 * params.numGC));
			BBrange realcells1 (params.numGC, params.Ny + params.numGC, params.numGC, params.Nx + params.numGC);
			BBrange realcells2 (params.numGC, params.Ny + params.numGC, params.numGC, params.Nx + params.numGC);
			set_sim_ICs(SF, params, eosparams, grid1, grid2, *ls);
			ls->initialisation();
			
			numsteps = 0;
			t = 0.0;
					
			std::cout << "Initialisation of GFM complete. Beginning time iterations with CFL = " << SF.CFL << "." << std::endl;
				
			while (t < params.T)
			{
				CFL = (numsteps < 5) ? std::min(SF.CFL, 0.2) : SF.CFL;
				dt = compute_dt(CFL, params, grid1, grid2, *ls, eosparams, params.T, t);
		
				GFM->set_ghost_states(params, *ls, eosparams, t, grid1, grid2, realcells1, realcells2);
				FS->pure_fluid_update(params, eosparams.gamma1, eosparams.pinf1, dt, realcells1, grid1, future_grid1);
				FS->pure_fluid_update(params, eosparams.gamma2, eosparams.pinf2, dt, realcells2, grid2, future_grid2);
				ls->update_interface_coupled(*(GFM->vfield), t, dt);
		
				apply_BCs_euler (params, grid1);
				apply_BCs_euler (params, grid2);
				
				numsteps++;
				t += dt;
				
				std::cout << "Time step " << numsteps << " complete. t = " << t << std::endl;
			}
			
			// Call function to get the L1, L2 and Linf errors here
		}
		
		outfile.close();
	}
	
	
}


double sim_circularexplosiontest :: compute_onefluid_dt (double CFL, const sim_info& params, const grideuler2type& grid, const binarySGparams& eosparams, double T, double t)
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
