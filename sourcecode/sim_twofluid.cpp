#include "sim_twofluid.hpp"
#include "levelset.hpp"
#include "GFM_ITM_interface.hpp"
#include "GFM_base.hpp"
#include "GFM_original.hpp"
#include "GFM_real.hpp"
#include "GFM_riemann.hpp"
#include "GFM_modified.hpp"
#include "euler_bc.hpp"
#include "euler_misc.hpp"
#include "pure_RS_exact.hpp"
#include "mixed_RS_exact.hpp"
#include "FS_godunov.hpp"
#include "FS_MUSCL.hpp"




void sim_twofluid :: run_sim (GFM_settingsfile SF)
{
	sim_info params;
	binarySGparams eosparams;
	set_sim_parameters(SF, params, eosparams);
	
	
	std::shared_ptr<FS_base> FS;
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
	

	int numsteps = 0;
	double t = 0.0;
	double CFL, dt;
	
	output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));

	std::cout << "[" << SF.basename << "] Initialisation complete. Beginning time iterations with CFL = " << SF.CFL << "." << std::endl;

	double tstart = omp_get_wtime();

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
		if (numsteps % params.output_freq == 0) output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}
	
	double tend = omp_get_wtime();

	
	
	output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));

	std::cout << "[" << SF.basename << "] Simulation complete." << std::endl;
	std::cout << "[" << SF.basename << "] Simulation took " << tend - tstart << "s." << std::endl;
}




double sim_twofluid :: compute_dt 
(
	double CFL, 
	const sim_info& params, 
	const grideuler2type& grid1, 
	const grideuler2type& grid2, 
	const GFM_ITM_interface& ls, 
	const binarySGparams& eosparams, 
	double T, 
	double t
)
{
	double maxu = 0.0, maxv = 0.0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double phi = ls.get_sdf(i,j);
			
			if (phi < 0.0)
			{
				double rho = grid1[i][j](0);
				double e = misc::specific_ie(grid1[i][j]);
				double p = eos::pressure(eosparams.gamma1, eosparams.pinf1, e, rho);
				double soundspeed = eos::soundspeed(eosparams.gamma1, eosparams.pinf1, p, rho);
				
				maxu = std::max(maxu, fabs(grid1[i][j](1))/rho + soundspeed);
				maxv = std::max(maxv, fabs(grid1[i][j](2))/rho + soundspeed);
			}
			else
			{
				double rho = grid2[i][j](0);
				double e = misc::specific_ie(grid2[i][j]);
				double p = eos::pressure(eosparams.gamma2, eosparams.pinf2, e, rho);
				double soundspeed = eos::soundspeed(eosparams.gamma2, eosparams.pinf2, p, rho);
				
				maxu = std::max(maxu, fabs(grid2[i][j](1))/rho + soundspeed);
				maxv = std::max(maxv, fabs(grid2[i][j](2))/rho + soundspeed);
			}
		}
	}
			
	double dt = CFL * std::min(params.dx/maxu, params.dy/maxv);

	if (t + dt > T) dt = T - t;

	assert(dt > 0.0);
	return dt;
}





void sim_twofluid :: output 
(
	int numsteps, 
	double t, 
	const sim_info& params, 
	const binarySGparams& eosparams,
	const grideuler2type& grid1, 
	const grideuler2type& grid2, 
	const GFM_ITM_interface& ls, 
	std::string filename,
	const velocity_field_base& vfield
)
{
	std::ofstream outfile;
	outfile.open(filename + "-" + std::to_string(numsteps) + "-prims.dat");
	std::ofstream outfile2;
	outfile2.open(filename + "-" + std::to_string(numsteps) + "-ls.dat");
	std::ofstream outfile3;
	outfile3.open(filename + "-" + std::to_string(numsteps) + "-vfield.dat");
	
	for (int i=1; i<params.Ny + 2 * params.numGC - 1; i++)
	{
		for (int j=1; j<params.Nx + 2 * params.numGC - 1; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);
			
			double phi = ls.get_sdf(i,j);
			outfile2 << CC(0) << " " << CC(1) << " " << phi << std::endl;
			outfile3 << CC(0) << " " << CC(1) << " " << vfield.get_velocity(CC, t)(0) << " " << vfield.get_velocity(CC, t)(1) << std::endl;
			
			if (phi < 0.0)
			{
				double rho = grid1[i][j](0);
				double u = grid1[i][j](1) / rho;
				double v = grid1[i][j](2) / rho;
				double e = misc::specific_ie(grid1[i][j]);
				double p = eos::pressure(eosparams.gamma1, eosparams.pinf1, e, rho);
				
				outfile << CC(0) << " " << CC(1) << " " << rho << " " << u << " " << v << " " << e << " " << p << std::endl;
			}
			else
			{
				double rho = grid2[i][j](0);
				double u = grid2[i][j](1) / rho;
				double v = grid2[i][j](2) / rho;
				double e = misc::specific_ie(grid2[i][j]);
				double p = eos::pressure(eosparams.gamma2, eosparams.pinf2, e, rho);
				
				outfile << CC(0) << " " << CC(1) << " " << rho << " " << u << " " << v << " " << e << " " << p << std::endl;
			}
		}
		outfile << std::endl;
		outfile2 << std::endl;
		outfile3 << std::endl;
	}
	
	outfile.close();	
	outfile2.close();
	outfile3.close();
}




void sim_twofluid :: set_sim_parameters 
(
	const GFM_settingsfile& SF, 
	sim_info& params, 
	binarySGparams& eosparams
)
{
	assert(SF.sim_type == "twofluid");
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
	params.dx = 1.0/params.Nx;
	params.dy = 1.0/params.Ny;
	params.BC_L = "transmissive";
	params.BC_T = "transmissive";
	params.BC_R = "transmissive";
	params.BC_B = "transmissive";
	
	if (SF.test_case == "circularexplosion")
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




void sim_twofluid :: set_sim_methods 
(
	const GFM_settingsfile& SF, 
	const sim_info& params,
	std::shared_ptr<FS_base>& FS,
	std::shared_ptr<GFM_base>& GFM,
	std::shared_ptr<GFM_ITM_interface>& ls
)
{
	std::shared_ptr<pure_RS_base> RS = std::make_shared<pure_RS_exact>();
	
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
		assert(!"Invalid FS choice");
	}
	
	
	std::shared_ptr<mixed_RS_base> mixed_RS = std::make_shared<mixed_RS_exact>();
	
	if (SF.GFM == "OGFM")
	{
		GFM = std::make_shared<GFM_original>(params);
	}
	else if (SF.GFM == "riemannGFM")
	{
		GFM = std::make_shared<GFM_riemann>(params);
	}
	else if (SF.GFM == "realGFM")
	{
		GFM = std::make_shared<GFM_real>(params);
	}
	else if (SF.GFM == "MGFM")
	{
		GFM = std::make_shared<GFM_modified>(params);
	}
	else
	{
		assert(!"Invalid GFM choice");
	}
	
	
	if (SF.ITM == "LS1")
	{
		ls = std::make_shared<Levelset>(params, 1);
	}
	else if (SF.ITM == "LS3")
	{
		ls = std::make_shared<Levelset>(params, 3);
	}
	else if (SF.ITM == "LS5")
	{
		ls = std::make_shared<Levelset>(params, 5);
	}
	else if (SF.ITM == "LS9")
	{
		ls = std::make_shared<Levelset>(params, 9);
	}
	else
	{
		assert(!"Invalid ITM choice");
	}
}




void sim_twofluid :: set_sim_ICs 
(
	const GFM_settingsfile& SF, 
	const sim_info& params, 
	const binarySGparams& eosparams, 
	grideuler2type& grid1, 
	grideuler2type& grid2, 
	GFM_ITM_interface& ls
)
{
	if (SF.test_case == "circularexplosion")
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
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);
		
		Eigen::Vector2d centre;
		centre << 1.0, 1.0;
		double R = 0.4;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				double lsval = (CC - centre).norm() - R;
				
				ls.set_sdf(i, j, lsval);
				
				grid1[i][j] = Lstate;
				grid2[i][j] = Rstate;
			}
		}
	}
	else
	{
		assert(!"Invalid test_case choice.");
	}
	
	apply_BCs_euler (params, grid1);
	apply_BCs_euler (params, grid2);
}

