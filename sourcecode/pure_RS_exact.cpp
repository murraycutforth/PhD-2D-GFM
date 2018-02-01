#include "pure_RS_exact.hpp"
#include "exact_RS_stiffenedgas.hpp"
#include <cassert>
#include <blitz/array.h>


void pure_RS_exact :: solve_RP (	
	const bool paralleltox
	const double gamma,
	const double pinf,
	const vectype& Lstate, 
	const vectype& Rstate, 
	vectype& flux,			
) const
{
	assert(misc::is_physical_state(gamma, pinf, Lstate));
	assert(misc::is_physical_state(gamma, pinf, Rstate));
	
	vectype Lprims = misc::conserved_to_primitives(gamma, pinf, Lstate);
	vectype Rprims = misc::conserved_to_primitives(gamma, pinf, Rstate);
	vectype soln (4);
	
	blitz::Array<double,1> Lprims_blitz (4);
	blitz::Array<double,1> Rprims_blitz (4);
	blitz::Array<double,1> soln_blitz (4);
	
	for (int i=0; i<4; i++)
	{
		Lprims_blitz(i) = Lprims(i);
		Rprims_blitz(i) = Rprims(i);
	}
	
	exact_rs_stiffenedgas RS (gamma, gamma, pinf, pinf);
	
	if (paralleltox)
	{
		RS.celledge_primitives_2D(Lprims_blitz, Rprims_blitz, soln_blitz);
		
		for (int i=0; i<4; i++)
		{
			soln(i) = soln_blitz(i);
		}
		
		flux = misc::flux_primitive_var(gamma, pinf, soln);
	}
	else
	{
		double temp = Lprims_blitz(1);
		Lprims_blitz(1) = Lprims_blitz(2);
		Lprims_blitz(2) = temp;
		temp = Rprims_blitz(1);
		Rprims_blitz(1) = Rprims_blitz(2);
		Rprims_blitz(2) = temp;
		
		RS.celledge_primitives_2D(Lprims_blitz, Rprims_blitz, soln_blitz);
		
		for (int i=0; i<4; i++)
		{
			soln(i) = soln_blitz(i);
		}
		
		flux = misc::flux_primitive_var(gamma, pinf, soln);
		
		double temp = flux(1);
		flux(1) = flux(2);
		flux(2) = temp;
	}
}
