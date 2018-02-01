#ifndef EXACT_SG_RS_H
#define EXACT_SG_RS_H

#include "pure_RS_base.hpp"

class pure_RS_exact : public pure_RS_base {
	
	public:

	virtual void solve_RP (	
		const bool paralleltox
		const double gamma,
		const double pinf,
		const vectype& Lstate, 
		const vectype& Rstate, 
		vectype& flux,			
	) const;
	
};

#endif
