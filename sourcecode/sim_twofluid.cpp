#include "sim_twofluid.hpp"
#include "levelset.hpp"
#include "VOF.hpp"
#include "MOF.hpp"
#include "MOF_advection.hpp"
#include "EMOF2_reconstruction.hpp"
#include "CLSVOF.hpp"
#include "GFM_ITM_interface.hpp"
#include "GFM_base.hpp"
#include "GFM_original.hpp"
#include "GFM_real.hpp"
#include "GFM_VOF.hpp"
#include "GFM_riemann.hpp"
#include "GFM_modified.hpp"
#include "OGFM_VOF.hpp"
#include "euler_bc.hpp"
#include "velocity_field_ODE_solver.hpp"
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
	double outputinterval = params.T / params.output_freq;
	double lastoutputtime = 0.0;
	
	output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));
	output_mass_error(numsteps, t, params, grid1, grid2, *ls, SF.basename);
	if (SF.test_case == "RMI_pert") output_perturbation_amplitude(numsteps, t, 0.0, *(GFM->vfield), params, grid1, grid2, *ls, SF.basename);

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
		if (t - lastoutputtime > outputinterval) 
		{
			output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));
			lastoutputtime = t;
		}
		output_mass_error(numsteps, t, params, grid1, grid2, *ls, SF.basename);
		if (SF.test_case == "RMI_pert") output_perturbation_amplitude(numsteps, t, dt, *(GFM->vfield), params, grid1, grid2, *ls, SF.basename);
		
		std::cout << "[" << SF.basename << "] Time step " << numsteps << " complete. t = " << t << std::endl;
	}
	
	double tend = omp_get_wtime();

	
	
	output(numsteps, t, params, eosparams, grid1, grid2, *ls, SF.basename, *(GFM->vfield));

	std::cout << "[" << SF.basename << "] Simulation complete." << std::endl;
	std::cout << "[" << SF.basename << "] Simulation took " << tend - tstart << "s." << std::endl;
	
	std::ofstream outfile;
	outfile.open(SF.basename + "-runtime.dat");
	outfile << tend - tstart << std::endl;
	outfile.close();
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
	outfile4.open(filename + "-" + std::to_string(numsteps) + "-z.dat");
	std::ofstream outfile5;
	outfile5.open(filename + "-" + std::to_string(numsteps) + "-schlieren.dat");
	std::ofstream outfile6;
	outfile6.open(filename + "-GDAerr.dat");
	
	
	for (int i=2; i<params.Ny + 2 * params.numGC - 2; i++)
	{
		for (int j=2; j<params.Nx + 2 * params.numGC - 2; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);
			
			double phi = ls.get_sdf(i,j);
			outfile2 << CC(0) << " " << CC(1) << " " << phi << std::endl;
			outfile3 << CC(0) << " " << CC(1) << " " << vfield.get_velocity(CC, t)(0) << " " << vfield.get_velocity(CC, t)(1) << std::endl;
			outfile4 << CC(0) << " " << CC(1) << " " << ls.get_z(i,j) << std::endl;
			
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
		outfile4 << std::endl;
	}
	
	outfile.close();	
	outfile2.close();
	outfile3.close();
	outfile4.close();
	
	
	std::vector<double> allgrads;
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double phi_L = ls.get_sdf(i,j-1);
			double phi_R = ls.get_sdf(i,j+1);
			double phi_T = ls.get_sdf(i+1,j);
			double phi_B = ls.get_sdf(i-1,j);
			
			double rho_L = (phi_L <= 0.0) ? grid1[i][j-1](0) : grid2[i][j-1](0);
			double rho_R = (phi_R <= 0.0) ? grid1[i][j+1](0) : grid2[i][j+1](0);
			double rho_T = (phi_T <= 0.0) ? grid1[i+1][j](0) : grid2[i+1][j](0);
			double rho_B = (phi_B <= 0.0) ? grid1[i-1][j](0) : grid2[i-1][j](0);
			
			double gradx = (rho_R - rho_L) / (2.0 * params.dx);
			double grady = (rho_T - rho_B) / (2.0 * params.dy);

			allgrads.push_back(sqrt(gradx*gradx + grady*grady));
		}
	}

	double maxgrad = *std::max_element(allgrads.begin(), allgrads.end());
			
	int counter = 0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);

			outfile5 << CC(0) << " " << CC(1) << " " << allgrads[counter] / maxgrad << std::endl;
			counter++;
		}
		outfile5 << std::endl;
	}
	
	outfile5.close();
	
	
	double L1err = 0.0;
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			Eigen::Vector2d mu;
			mu << 0.5, 0.5;
			
			double rho = grid1[i][j](0);
			double rho_exact = 1.0 + 1000.0 * exp(- ((cc - mu).squaredNorm()) / (2.0 * 0.1 * 0.1));
			
			L1err += fabs(rho - rho_exact) * params.dx * params.dy;
		}
	}
	
	outfile6 << params.Nx << " " << L1err << std::endl;
	outfile6.close();
}





void sim_twofluid :: output_mass_error
(
	int numsteps,
	double t,
	const sim_info& params,
	const grideuler2type& grid1, 
	const grideuler2type& grid2, 
	const GFM_ITM_interface& ls, 
	std::string filename
)
{
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
}




void sim_twofluid :: set_sim_parameters 
(
	const GFM_settingsfile& SF, 
	sim_info& params, 
	binarySGparams& eosparams
)
{
	assert(SF.sim_type == "twofluid");
	params.numGC = 8;
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.output_freq = 4;
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
	else if (SF.test_case == "GDA")
	{
		
		params.BC_L = "periodic";
		params.BC_T = "periodic";
		params.BC_R = "periodic";
		params.BC_B = "periodic";
		params.T = 10.0;
		
		params.output_freq = 0;
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		eosparams.gamma1 = 1.667;

		params.dx = 325.0 / params.Nx;
		params.dy = 89.0 / params.Ny;
		params.T = 280.0;
		params.BC_T = "reflective";
		params.BC_B = "reflective";
		
		params.output_freq = 2;
	}
	else if (SF.test_case == "shocked_R22_bubble")
	{
		eosparams.gamma1 = 1.249;

		params.dx = 0.445 / params.Nx;
		params.dy = 0.089 / params.Ny;
		params.T = 0.0014;
		params.BC_T = "reflective";
		params.BC_B = "reflective";
		
		params.output_freq = 20;
	}
	else if (SF.test_case == "shocked_SF6")
	{
		eosparams.gamma2 = 1.076;

		params.dx = 0.45 / params.Nx;
		params.dy = 0.2 / params.Ny;
		params.T = 0.2;
		params.BC_T = "reflective";
		params.BC_R = "reflective";
		params.BC_B = "reflective";
		params.output_freq = 10;
	}
	else if (SF.test_case == "RMI" || SF.test_case == "RMI_pert" || SF.test_case == "RMI_unpert")
	{
		eosparams.gamma2 = 1.093;

		params.dx = 4.0 / params.Nx;
		params.dy = 1.0 / params.Ny;
		params.T = 10.0;
		params.BC_T = "reflective";
		params.BC_B = "reflective";
		params.output_freq = 10;
		
		if (SF.test_case == "RMI_pert") params.T = 4.0;
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		eosparams.gamma1 = 7.15;
		eosparams.pinf1 = 3309.0;

		params.dx = 12.0 / params.Nx;
		params.dy = 12.0 / params.Ny;
		params.T = 0.05;
		params.output_freq = 8;
	}
	else if (SF.test_case == "underwater_explosion")
	{
		eosparams.gamma2 = 7.15;
		eosparams.pinf2 = 3.309e8;
		params.output_freq = 10;

		params.x0 = -5.0;
		params.y0 = -5.0;
		params.dx = 10.0 / params.Nx;
		params.dy = 10.0 / params.Ny;
		params.T = 0.005;
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
	else if (SF.test_case == "TSTM")
	{
		eosparams.gamma1 = 1.5;
		
		params.dx = 7.0/params.Nx;
		params.dy = 3.0/params.Ny;
		params.T = 8.0;
		params.output_freq = 10;
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
	else if (SF.GFM == "VOFGFM")
	{
		GFM = std::make_shared<GFM_VOF>(params);
	}
	else if (SF.GFM == "OGFMVOF")
	{
		GFM = std::make_shared<OGFM_VOF>(params);
	}
	else
	{
		assert(!"Invalid GFM choice");
	}
	
	GFM->use_extension_vfield = false;
	GFM->use_slip_BCs = false;
	
	
	if (SF.ITM == "LS1")
	{
		ls = std::make_shared<Levelset>(params, 1, 1);
	}
	else if (SF.ITM == "LS3")
	{
		ls = std::make_shared<Levelset>(params, 3);
	}
	else if (SF.ITM == "LS5")
	{
		ls = std::make_shared<Levelset>(params, 5, 1);
	}
	else if (SF.ITM == "LS9")
	{
		ls = std::make_shared<Levelset>(params, 9, 1);
	}
	else if (SF.ITM == "CLSVOF")
	{
		ls = std::make_shared<CLSVOF>(params);
	}
	else if (SF.ITM == "PYVOF")
	{
		ls = std::make_shared<VOF>(params);
		GFM->use_extension_vfield = false;
	}
	else if (SF.ITM == "MOF")
	{
		ls = std::make_shared<MOF>(params);
		GFM->use_extension_vfield = false;
	}
	else if (SF.ITM == "EMOF2")
	{
		ls = std::make_shared<MOF>(params, std::make_shared<EMOF2_reconstruction>(6400,1));
		GFM->use_extension_vfield = false;
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
	else if (SF.test_case == "GDA")
	{
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				
				Eigen::Vector2d mu;
				mu << 0.5, 0.5;
				
				double sigma = 0.1;
				
				vec4type prims;
				
				prims(0) = 1.0 + 1000.0 * exp(- ((cc - mu).squaredNorm()) / (2.0 * sigma * sigma));
				prims(1) = 1.0;
				prims(2) = 1.0;
				prims(3) = 0.0001;
				
				vec4type state = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, prims);
				
				grid1[i][j] = state;
				grid2[i][j] = state;
				
				ls.set_sdf(i, j, -1.0);
				ls.set_z(i, j, 1.0);
			}
		}
			
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		// Fluid 1 is helium 
		
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
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);
		vec4type hestate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, heprims);
		
		Eigen::Vector2d centre;
		centre << 175.0, 44.5;
		double R = 25.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				
				
				// Set air state
				
				if (CC(0) < 225.0) grid2[i][j] = Rstate;
				else grid2[i][j] = Lstate;
				
				
				// Set he state
				
				grid1[i][j] = hestate;
				
				
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
				
				double zval = frac;
				
				double lsval = (CC - centre).norm() - R;
				ls.set_sdf(i, j, lsval);
				ls.set_z(i, j, zval);
			}
		}
	}
	else if (SF.test_case == "shocked_R22_bubble")
	{
		// Fluid 1 is R22 
		
		vec4type r22prims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);
		
		r22prims(0) = 3.863;
		r22prims(1) = 0.0;
		r22prims(2) = 0.0;
		r22prims(3) = 1.01325e5;
		
		Lprims(0) = 1.686;
		Lprims(1) = -113.5;
		Lprims(2) = 0.0;
		Lprims(3) = 1.59e5;
		
		Rprims(0) = 1.225;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 1.01325e5;
		
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, Rprims);
		vec4type r22state = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, r22prims);
		
		Eigen::Vector2d centre;
		centre << 0.225, 0.0445;
		double R = 0.025;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d CC = params.cellcentre_coord(i, j);
				
				
				// Set air state
				
				if (CC(0) < 0.275) grid2[i][j] = Rstate;
				else grid2[i][j] = Lstate;
				
				
				// Set he state
				
				grid1[i][j] = r22state;
				
				
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
				
				double zval = frac;
				
				double lsval = (CC - centre).norm() - R;
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
	else if (SF.test_case == "RMI")
	{
		// Air is fluid 1
		
		vec4type sf6prims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);

		sf6prims(0) = 5.04;
		sf6prims(1) = 0.0;
		sf6prims(2) = 0.0;
		sf6prims(3) = 1.0;
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1.0;
		
		Rprims(0) = 1.411;
		Rprims(1) = -0.39;
		Rprims(2) = 0.0;
		Rprims(3) = 1.628;
		
		vec4type sf6state = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, sf6prims);
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (cc(0) < 3.2)
				{
					grid1[i][j] = Lstate;
				}
				else
				{
					grid1[i][j] = Rstate;
				}

				grid2[i][j] = sf6state;
				
				
				// Set z as fraction of area inside sinusoidal interface

				double x1 = 2.9;
				double eps = 0.2;
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double pi = atan(1.0) * 4.0;
				double dely = params.dy/numsamples;
				
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
						
						double interfacex = x1 - eps * sin(2.0 * pi * (samplepos(1) + 0.25));
						
						
						if (interfacex > samplepos(0)) 
						{
							numinside++;
						}
					}
				}

				double zval = 1.0 - double(numinside)/totalnumsamples;

				ls.set_z(i, j, zval);


				// Set level set by iterating along sinusoidal interface

				double mindist = 1e100;

				for (double ypos = 0.0; ypos <=0.5; ypos += params.dy * 0.1)
				{
					double xpos = x1 - eps * sin(2.0 * M_PI * (ypos + 0.25));

					Eigen::Vector2d interfacepos;
					interfacepos << xpos, ypos;

					mindist = std::min(mindist, (interfacepos - cc).norm());
				}

				double lsval = std::copysign(mindist, 0.5 - zval);

				ls.set_sdf(i, j, lsval);
			}
		}
	}
	else if (SF.test_case == "RMI_pert")
	{
		// Air is fluid 1
		
		vec4type sf6prims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);

		sf6prims(0) = 5.04;
		sf6prims(1) = 0.0;
		sf6prims(2) = 0.0;
		sf6prims(3) = 1.0;
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1.0;
		
		Rprims(0) = 1.411;
		Rprims(1) = -0.39;
		Rprims(2) = 0.0;
		Rprims(3) = 1.628;
		
		vec4type sf6state = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, sf6prims);
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (cc(0) < 3.2)
				{
					grid1[i][j] = Lstate;
				}
				else
				{
					grid1[i][j] = Rstate;
				}

				grid2[i][j] = sf6state;
				
				
				// Set z as fraction of area inside sinusoidal interface

				double x1 = 2.9;
				double eps = 0.02;
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double pi = atan(1.0) * 4.0;
				double dely = params.dy/numsamples;
				
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
						
						double interfacex = x1 - eps * sin(2.0 * pi * (samplepos(1) + 0.25));
						
						
						if (interfacex > samplepos(0)) 
						{
							numinside++;
						}
					}
				}

				double zval = 1.0 - double(numinside)/totalnumsamples;

				ls.set_z(i, j, zval);


				// Set level set by iterating along sinusoidal interface

				double mindist = 1e100;

				for (double ypos = 0.0; ypos <=0.5; ypos += params.dy * 0.1)
				{
					double xpos = x1 - eps * sin(2.0 * M_PI * (ypos + 0.25));

					Eigen::Vector2d interfacepos;
					interfacepos << xpos, ypos;

					mindist = std::min(mindist, (interfacepos - cc).norm());
				}

				double lsval = std::copysign(mindist, 0.5 - zval);

				ls.set_sdf(i, j, lsval);
			}
		}
	}
	else if (SF.test_case == "RMI_unpert")
	{
		// Air is fluid 1
		
		vec4type sf6prims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);

		sf6prims(0) = 5.04;
		sf6prims(1) = 0.0;
		sf6prims(2) = 0.0;
		sf6prims(3) = 1.0;
		
		Lprims(0) = 1.0;
		Lprims(1) = 0.0;
		Lprims(2) = 0.0;
		Lprims(3) = 1.0;
		
		Rprims(0) = 1.411;
		Rprims(1) = -0.39;
		Rprims(2) = 0.0;
		Rprims(3) = 1.628;
		
		vec4type sf6state = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, sf6prims);
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);

		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (cc(0) < 3.2)
				{
					grid1[i][j] = Lstate;
				}
				else
				{
					grid1[i][j] = Rstate;
				}

				grid2[i][j] = sf6state;
				
				
				// Set z as fraction of area inside sinusoidal interface

				double x1 = 2.9;
				double eps = 0.0;
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double pi = atan(1.0) * 4.0;
				double dely = params.dy/numsamples;
				
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
						
						double interfacex = x1 - eps * sin(2.0 * pi * (samplepos(1) + 0.25));
						
						
						if (interfacex > samplepos(0)) 
						{
							numinside++;
						}
					}
				}

				double zval = 1.0 - double(numinside)/totalnumsamples;

				ls.set_z(i, j, zval);


				// Set level set by iterating along sinusoidal interface

				double mindist = 1e100;

				for (double ypos = 0.0; ypos <=0.5; ypos += params.dy * 0.1)
				{
					double xpos = x1 - eps * sin(2.0 * M_PI * (ypos + 0.25));

					Eigen::Vector2d interfacepos;
					interfacepos << xpos, ypos;

					mindist = std::min(mindist, (interfacepos - cc).norm());
				}

				double lsval = std::copysign(mindist, 0.5 - zval);

				ls.set_sdf(i, j, lsval);
			}
		}
	}
	else if (SF.test_case == "shocked_SF6")
	{
		// Air is fluid 1
		
		vec4type sf6prims (4);
		vec4type Lprims (4);
		vec4type Rprims (4);

		sf6prims(0) = 5.805;
		sf6prims(1) = 0.0;
		sf6prims(2) = 0.0;
		sf6prims(3) = 9.6856;
		
		Lprims(0) = 1.6672;
		Lprims(1) = 1.33273;
		Lprims(2) = 0.0;
		Lprims(3) = 16.3256;
		
		Rprims(0) = 1.153;
		Rprims(1) = 0.0;
		Rprims(2) = 0.0;
		Rprims(3) = 9.6856;
		
		vec4type sf6state = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, sf6prims);
		vec4type Lstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims);
		vec4type Rstate = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims);
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				grid2[i][j] = sf6state;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
								
				if (cc(0) < 0.05)
				{
					grid1[i][j] = Lstate;
				}
				else
				{
					grid1[i][j] = Rstate;
				}
				
				
				// Set z as fraction of area outside rectangular region [0.1, 0.25] x [0.0, 0.1]
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
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
						
						if (samplepos(0) >= 0.1 && samplepos(0) <= 0.25
							&& samplepos(1) >= 0.0 && samplepos(1) <= 0.1) 
						{
							numinside++;
						}
					}
				}

				double z = 1.0 - double(numinside)/totalnumsamples;
				ls.set_z(i, j, z);
				
				
				// Set level set
				
				double y1 = 0.1 - cc(1);
				double xL = cc(0) - 0.1;
				double xR = 0.25 - cc(0);
				Eigen::Vector2d TL;
				TL << 0.1, 0.1;
				Eigen::Vector2d TR;
				TR << 0.25, 0.1;
				
				if (0.1 <= cc(0) && cc(0) <= 0.25 && cc(1) <= 0.1)
				{
					ls.set_sdf(i, j, std::min({y1, xL, xR}));
				}
				else if (cc(0) <= 0.1)
				{
					if (cc(1) <= 0.1)
					{
						ls.set_sdf(i, j, xL);
					}
					else
					{
						ls.set_sdf(i, j, -(cc - TL).norm());
					}
				}
				else if (cc(0) >= 0.25)
				{
					if (cc(1) <= 0.1)
					{
						ls.set_sdf(i, j, xR);
					}
					else
					{
						ls.set_sdf(i, j, -(cc - TR).norm());
					}
				}
				else
				{
					ls.set_sdf(i, j, y1);
				}
			}
		}
	}
	else if (SF.test_case == "TSTM")
	{
		vec4type prims2 (4);
		vec4type Lprims1 (4);
		vec4type Rprims1 (4);

		prims2(0) = 1.0;
		prims2(1) = 0.0;
		prims2(2) = 0.0;
		prims2(3) = 0.1;
		
		Lprims1(0) = 1.0;
		Lprims1(1) = 0.0;
		Lprims1(2) = 0.0;
		Lprims1(3) = 1.0;
		
		Rprims1(0) = 0.125;
		Rprims1(1) = 0.0;
		Rprims1(2) = 0.0;
		Rprims1(3) = 0.1;
		
		vec4type state2 = misc::primitives_to_conserved(eosparams.gamma2, eosparams.pinf2, prims2);
		vec4type Lstate1 = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Lprims1);
		vec4type Rstate1 = misc::primitives_to_conserved(eosparams.gamma1, eosparams.pinf1, Rprims1);
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				grid2[i][j] = state2;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
								
				if (cc(0) < 1.0)
				{
					grid1[i][j] = Lstate1;
				}
				else
				{
					grid1[i][j] = Rstate1;
				}
				
				
				// Set z as fraction of area outside rectangular region [1, 7] x [0.0, 1.5]
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
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
						
						if (samplepos(0) >= 1.0 && samplepos(0) <= 8.0
							&& samplepos(1) >= -1.0 && samplepos(1) <= 1.5) 
						{
							numinside++;
						}
					}
				}

				double z = 1.0 - double(numinside)/totalnumsamples;
				ls.set_z(i, j, z);
				
				
				// Set level set
				
				double y1 = 1.5 - cc(1);
				double xL = cc(0) - 1.0;
				double xR = 8.0 - cc(0);
				Eigen::Vector2d TL;
				TL << 1.0, 1.5;
				Eigen::Vector2d TR;
				TR << 8.0, 1.5;
				
				if (1.0 <= cc(0) && cc(0) <= 8.0 && cc(1) <= 1.5)
				{
					ls.set_sdf(i, j, std::min({y1, xL, xR}));
				}
				else if (cc(0) <= 1.0)
				{
					if (cc(1) <= 1.5)
					{
						ls.set_sdf(i, j, xL);
					}
					else
					{
						ls.set_sdf(i, j, -(cc - TL).norm());
					}
				}
				else if (cc(0) >= 8.0)
				{
					if (cc(1) <= 1.5)
					{
						ls.set_sdf(i, j, xR);
					}
					else
					{
						ls.set_sdf(i, j, -(cc - TR).norm());
					}
				}
				else
				{
					ls.set_sdf(i, j, y1);
				}
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



void sim_twofluid :: output_perturbation_amplitude
(
	int numsteps,
	double t,
	double dt,
	const solved_velocity_field_base& vfield,
	const sim_info& params,
	const grideuler2type& grid1,
	const grideuler2type& grid2,
	const GFM_ITM_interface& ls,
	std::string filename
)
{
	static Eigen::Vector2d x0;
	static Eigen::Vector2d x1;
	
	if (numsteps == 0)
	{
		std::ofstream outfile1;
		outfile1.open(filename + "-pertamp.dat");
		
		std::ofstream outfile2;
		outfile2.open(filename + "-pertampderiv.dat");
		
		x0 << 2.88, 1.0;
		x1 << 2.92, 0.5;
	}
	else
	{
		std::ofstream outfile1;
		outfile1.open(filename + "-pertamp.dat", std::ofstream::app);
		
		std::ofstream outfile2;
		outfile2.open(filename + "-pertampderiv.dat", std::ofstream::app);
		
		
		// Output current amplitude
		
		outfile1 << t << " " << 0.5 * (x1(0) - x0(0)) << " " << x1(0) << " " << x0(0) << std::endl;
		
		
		// Output rate of change of amplitude
		
		outfile2 << t << " " << 0.5 * (vfield.get_velocity(x1, t)(0) - vfield.get_velocity(x0, t)(0)) << " " << vfield.get_velocity(x1, t)(0) << " " << vfield.get_velocity(x0, t)(0) << std::endl;
		
		
		// Advect points in vfield
		
		time_step_RK4(vfield, t, dt, x0, x0);
		time_step_RK4(vfield, t, dt, x1, x1);
		
		// Keep points in correct y-plane
		
		x0(1) = 1.0;
		x1(1) = 0.5;
	}	
}

