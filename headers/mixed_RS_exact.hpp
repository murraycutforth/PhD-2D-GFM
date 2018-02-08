#ifndef MIXED_RS_EXACT_H
#define MIXED_RS_EXACT_H


#include "typedefs.hpp"
#include "mixed_RS_base.hpp"


class mixed_RS_exact : public mixed_RS_base {

	public:

	virtual void solve_mixed_RP 
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
	) const override;
};

#endif
