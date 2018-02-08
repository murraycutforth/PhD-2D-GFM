#include "pure_RS_exact.hpp"
#include "exact_RS_stiffenedgas.hpp"
//~ #include "exact_RS_idealgas.hpp"
#include "euler_misc.hpp"
#include <cassert>
#include <blitz/array.h>


void pure_RS_exact :: solve_RP (	
	const bool paralleltox,
	const double gamma,
	const double pinf,
	const vec4type& Lstate, 
	const vec4type& Rstate, 
	vec4type& flux
) const
{
	assert(misc::is_physical_state(gamma, pinf, Lstate));
	assert(misc::is_physical_state(gamma, pinf, Rstate));
	
	vec4type Lprims = misc::conserved_to_primitives(gamma, pinf, Lstate);
	vec4type Rprims = misc::conserved_to_primitives(gamma, pinf, Rstate);
	vec4type soln (4);
	
	blitz::Array<double,1> Lprims_blitz (3);
	blitz::Array<double,1> Rprims_blitz (3);
	blitz::Array<double,1> soln_blitz (3);
	
	exact_rs_stiffenedgas RS (gamma, gamma, pinf, pinf);
		
	if (paralleltox)
	{
		Lprims_blitz(0) = Lprims(0);
		Lprims_blitz(1) = Lprims(1);
		Lprims_blitz(2) = Lprims(3);
		Rprims_blitz(0) = Rprims(0);
		Rprims_blitz(1) = Rprims(1);
		Rprims_blitz(2) = Rprims(3);
				
		RS.solve_RP(Lprims_blitz, Rprims_blitz);
		soln_blitz = RS.sample_solution(Lprims_blitz, Rprims_blitz, 0.0);
		
		soln(0) = soln_blitz(0);
		soln(1) = soln_blitz(1);
		soln(2) = RS.S_STAR <= 0.0 ? Rprims(2) : Lprims(2);
		soln(3) = soln_blitz(2);
		
		flux = misc::flux_primitive_var(gamma, pinf, soln);
	}
	else
	{
		Lprims_blitz(0) = Lprims(0);
		Lprims_blitz(1) = Lprims(2);
		Lprims_blitz(2) = Lprims(3);
		Rprims_blitz(0) = Rprims(0);
		Rprims_blitz(1) = Rprims(2);
		Rprims_blitz(2) = Rprims(3);
		
		RS.solve_RP(Lprims_blitz, Rprims_blitz);
		soln_blitz = RS.sample_solution(Lprims_blitz, Rprims_blitz, 0.0);
		
		soln(0) = soln_blitz(0);
		soln(1) = soln_blitz(1);
		soln(2) = RS.S_STAR <= 0.0 ? Rprims(1) : Lprims(1);
		soln(3) = soln_blitz(2);
		
		flux = misc::flux_primitive_var(gamma, pinf, soln);
		
		double temp = flux(1);
		flux(1) = flux(2);
		flux(2) = temp;
	}
}
