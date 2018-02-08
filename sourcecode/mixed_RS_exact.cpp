#include "mixed_RS_exact.hpp"
#include "exact_RS_stiffenedgas.hpp"
#include "euler_misc.hpp"
#include <cassert>
#include <blitz/array.h>


void mixed_RS_exact :: solve_mixed_RP 
(
	const double gammaL,
	const double pinfL,
	const double gammaR,
	const double pinfR,
	const vec4type& Lstate, 
	const vec4type& Rstate, 
	double& p_star,
	double& u_star,
	double& rho_star_L,
	double& rho_star_R
) const
{
	assert(misc::is_physical_state(gammaL, pinfL, Lstate));
	assert(misc::is_physical_state(gammaR, pinfR, Rstate));
	
	vec4type Lprims = misc::conserved_to_primitives(gammaL, pinfL, Lstate);
	vec4type Rprims = misc::conserved_to_primitives(gammaR, pinfR, Rstate);
	vec4type soln (4);
	
	blitz::Array<double,1> Lprims_blitz (3);
	blitz::Array<double,1> Rprims_blitz (3);
	blitz::Array<double,1> soln_blitz (3);
	
	exact_rs_stiffenedgas RS (gammaL, gammaR, pinfL, pinfR);
	
	Lprims_blitz(0) = Lprims(0);
	Lprims_blitz(1) = Lprims(1);
	Lprims_blitz(2) = Lprims(3);
	Rprims_blitz(0) = Rprims(0);
	Rprims_blitz(1) = Rprims(1);
	Rprims_blitz(2) = Rprims(3);
				
	RS.solve_RP(Lprims_blitz, Rprims_blitz);
	
	p_star = RS.P_STAR;
	u_star = RS.S_STAR;
	rho_star_L = RS.rho_star_L;
	rho_star_R = RS.rho_star_R;
}


