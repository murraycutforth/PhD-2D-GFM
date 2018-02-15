#include "sim_twofluid.hpp"
#include "levelset.hpp"
#include "CLSVOF.hpp"
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
#include <cmath>



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
	ls->initialisation();
	

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
	std::ofstream outfile4;
	if (numsteps == 0) outfile4.open(filename + "-masschange.dat");
	else outfile4.open(filename + "-masschange.dat", std::ofstream::app);

	static double initialmass1;
	static double initialmass2;
	
	if (numsteps == 0)
	{
		initialmass1 = fluid_mass(-1.0, params, grid1, ls);
		initialmass2 = fluid_mass(1.0, params, grid2, ls);
	}

	double currentmass1 = fluid_mass(-1.0, params, grid1, ls);
	double currentmass2 = fluid_mass(1.0, params, grid2, ls);
	double masserror1 = (currentmass1 - initialmass1) / initialmass1;
	double masserror2 = (currentmass2 - initialmass2) / initialmass2;

	outfile4 << t << " " << masserror1 << " " << masserror2 << std::endl;
 
	
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
	else if (SF.test_case == "ST1_x")
	{
		eosparams.gamma2 = 1.2;
		params.T = 0.0007;
	}
	else if (SF.test_case == "ST2_y")
	{
		eosparams.gamma1 = 4.4;
		eosparams.pinf1 = 600000000.0;

		params.x0 = -2.0;
		params.y0 = -2.0;
		params.dx = 4.0 / params.Nx;
		params.dy = 4.0 / params.Ny;
		params.T = 0.0009;
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		eosparams.gamma2 = 1.667;

		params.dx = 325.0 / params.Nx;
		params.dy = 89.0 / params.Ny;
		params.T = 280.0;
		params.BC_T = "reflective";
		params.BC_B = "reflective";
		params.output_freq = params.Nx;
	}
	else if (SF.test_case == "shocked_SF6")
	{
		eosparams.gamma2 = 1.076;

		params.dx = 0.45 / params.Nx;
		params.dy = 0.2 / params.Ny;
		params.T = 0.25;
		params.BC_T = "reflective";
		params.BC_R = "reflective";
		params.BC_B = "reflective";
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		eosparams.gamma1 = 7.15;
		eosparams.pinf1 = 3309.0;

		params.dx = 12.0 / params.Nx;
		params.dy = 12.0 / params.Ny;
		params.T = 0.05;
		params.output_freq = params.Nx / 2;
	}
	else if (SF.test_case == "underwater_explosion")
	{
		eosparams.gamma2 = 7.15;
		eosparams.pinf2 = 3.309e8;

		params.x0 = -5.0;
		params.y0 = -5.0;
		params.dx = 10.0 / params.Nx;
		params.dy = 10.0 / params.Ny;
		params.T = 0.005;
		params.output_freq = params.Nx / 2;
	}
	else if (SF.test_case == "tin_air_implosion")
	{
		eosparams.gamma1 = 3.27;
		eosparams.pinf1 = 149500.0;
		
		params.x0 = 0.0;
		params.y0 = -25.0;
		params.dx = 25.0/params.Nx;
		params.dy = 50.0/params.Ny;
		params.T = 0.08;
		params.BC_L = "reflective";
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
	else if (SF.ITM == "CLSVOF")
	{
		ls = std::make_shared<CLSVOF>(params);
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
				
				int N = 10;
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
				
				double zval = frac;
				
				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);
				
				grid1[i][j] = Lstate;
				grid2[i][j] = Rstate;
			}
		}
	}
	else if (SF.test_case == "ST1_x")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 100000.0;
		
		Rprims(0) = 0.125;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 10000.0;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				double lsval = CC(0) - 0.5;
				double zval = std::min(1.0, std::max(0.0, (0.5 - CC(0) + 0.5 * params.dx) / params.dx));

				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);

				grid1[i][j] = Lstate;
				grid2[i][j] = Rstate;
			}
		}
	}
	else if (SF.test_case == "ST2_y")
	{
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		Lprims(0) = 1000.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1000000000.0;
		
		Rprims(0) = 50.0;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 100000.0;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				double lsval = CC(1) - 0.7;
				double zval = std::min(1.0, std::max(0.0, (0.7 - CC(1) + 0.5 * params.dy) / params.dy));

				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);

				grid1[i][j] = Lstate;
				grid2[i][j] = Rstate;
			}
		}
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		// Fluid 2 is helium 
		
		vec4type heprims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		heprims(0) = 0.138;
		heprims(1) = 0.0;
		heprims(2) = 0.0;
		heprims(3) = 1.0;
		
		Lprims(0) = 1.3764;
		Lprims(1) = -0.394;
		Lprims(2) = 0.0;
		Lprims(3) = 1.5698;
		
		Rprims(0) = 1.0;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 1.0;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		vec4type hestate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, heprims);
		
		Eigen::Vector2d centre;
		centre << 175.0, 44.5;
		double R = 25.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				
				
				// Set air state
				
				if (CC(0) < 225.0) grid1[i][j] = Rstate;
				else grid1[i][j] = Lstate;
				
				
				// Set he state
				
				grid2[i][j] = hestate;
				
				
				// Set ls
				
				int N = 10;
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
						
						if (samplepos.norm() > R) numinside++;
					}
				}
	
				double frac = double(numinside)/totalnumsamples;
				
				double zval = frac;
				
				double lsval = R - (CC - centre).norm();
				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);
			}
		}
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		// Fluid 1 is water
		
		vec4type airprims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		airprims(0) = 0.0012;
		airprims(1) = 0.0;
		airprims(2) = 0.0;
		airprims(3) = 1.0;
		
		Lprims(0) = 1.31;
		Lprims(1) = 67.32;
		Lprims(2) = 0.0;
		Lprims(3) = 19000.0;
		
		Rprims(0) = 1.0;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 1.0;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		vec4type airstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, airprims);
		
		Eigen::Vector2d centre;
		centre << 6.0, 6.0;
		double R = 3.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				
				
				// Set air state
				
				if (CC(0) < 2.4) grid1[i][j] = Lstate;
				else grid1[i][j] = Rstate;
				
				
				// Set he state
				
				grid2[i][j] = airstate;
				
				
				// Set ls
				
				int N = 10;
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
				
				double zval = 1.0 - frac;
				
				double lsval = R - (CC - centre).norm();
				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);
			}
		}
	}
	else if (SF.test_case == "underwater_explosion")
	{
		// Fluid 1 is air
		
		vec4type waterprims (4);
		vec4type bubbleprims (4);
		vec4type ambientprims (4);
		
		waterprims(0) = 1000.0;
		waterprims(1) = 0.0;
		waterprims(2) = 0.0;
		waterprims(3) = 1.0e5;
		
		bubbleprims(0) = 1270.0;
		bubbleprims(1) = 0.0;
		bubbleprims(2) = 0.0;
		bubbleprims(3) = 8.29e8;
		
		ambientprims(0) = 1.0;
		ambientprims(1) = 0.0;
		ambientprims(2) = 0.0;
		ambientprims(3) = 1.0e5;
		
		vec4type waterstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, waterprims);
		vec4type bubblestate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, bubbleprims);
		vec4type ambientstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, ambientprims);
		
		Eigen::Vector2d centre;
		centre << 0.0, 0.0;
		double R = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				
				
				// Set water state
				
				grid2[i][j] = waterstate;
				
				
				// Set air state
				
				if (CC(1) < 2.0)
				{
					grid1[i][j] = bubblestate;
				}
				else
				{
					grid1[i][j] = ambientstate;
				}
				
				
				// Set level set value
				
				double lsbubbleval = (CC - centre).norm() - R;
				double lssurfaceval = 2.5 - CC(1);
				double lsval = std::min(lsbubbleval, lssurfaceval);
				ls.set_sdf(i, j, lsval);
				
				
				// Set volume fraction
				
				if (CC(1) < 2.0)
				{
					int N = 10;
					int totalnumsamples = N*N;
					int numinside = 0;
					double delx = params.dx/N;
					double dely = params.dy/N;
					
					Eigen::Vector2d BL;
					BL(0) = CC(0) - 0.5 * params.dx;
					BL(1) = CC(1) - 0.5 * params.dy;
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
	
					double zval = double(numinside)/totalnumsamples;
					ls.set_z(i, j, zval);
				}
				else
				{
					double zval = 1.0 - std::min(1.0, std::max(0.0, (2.5 - CC(1) + 0.5 * params.dy) / params.dy));
					ls.set_z(i, j, zval);
				}
			}
		}
	}
	else if (SF.test_case == "tin_air_implosion")
	{
		// Fluid 2 is air
		
		vec4type airprims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);

		airprims(0) = 0.001;
		airprims(1) = 0.0;
		airprims(2) = 0.0;
		airprims(3) = 1.0;
		
		Lprims(0) = 11.84;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1000000.0;
		
		Rprims(0) = 7.28;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 1.0;
		
		vec4type airstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, airprims);
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{				
				
				// Set tin state as fraction of area inside circle of radius 24
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						if (samplepos.norm() <= 24.0) numinside++;
					}
				}
				
				double insideratio = double(numinside) / totalnumsamples;
				
				vec4type tinprims = insideratio * Rprims + (1.0 - insideratio) * Lprims;
				grid1[i][j] = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, tinprims);
				
				
				
				// Set air state
				
				grid2[i][j] = airstate;
				
				
				
				// Set z and phi
				
				numinside = 0;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						double theta = atan2(samplepos(1), samplepos(0));
						theta -= atan(1) * 2;
						double r_interface = 20.0 + 0.4 * cos(22 * theta) + 0.4 * cos(17 * theta) + 0.3 * cos(29 * theta);
						
						if (samplepos.norm() <= r_interface) numinside++;
					}
				}

				double z = 1.0 - double(numinside)/totalnumsamples;
				ls.set_z(i, j, z);
				
				double mindist = 1e100;
				
				for (double theta=0.0; theta<M_PI; theta += 0.005)
				{
					double angle = theta - atan(1) * 2;
					double r_interface = 20.0 + 0.4 * cos(22 * angle) + 0.4 * cos(17 * angle) + 0.3 * cos(29 * angle);
					
					Eigen::Vector2d interfacepos;
					interfacepos << r_interface * cos(angle), r_interface * sin(angle);
					
					double sign = std::copysign(1.0, r_interface - cc.norm());
					double dist = (cc - interfacepos).norm() * sign;
					
					mindist = fabs(dist) < fabs(mindist) ? dist : mindist;
				}
				
				ls.set_sdf(i, j, mindist);
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



double sim_twofluid :: fluid_mass
(
	const double interiorsign,
	const sim_info& params,
	const grideuler2type& grid,
	const GFM_ITM_interface& ls
)
{
	double totalmass = 0.0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double lsval = ls.get_sdf(i,j);

			if (interiorsign * lsval > 0.0)
			{
				totalmass += grid[i][j](0) * params.dx * params.dy;
			}
		}
	}

	return totalmass;
}


