#ifndef PURE_RS_BASE_H
#define PURE_RS_BASE_H


#include "typedefs.hpp"
#include "stiffened_gas_eos.hpp"


class pure_RS_base {

	public:

	virtual void solve_RP (	
		const bool paralleltox,
		const double gamma,
		const double pinf,
		const vec4type& Lstate, 
		const vec4type& Rstate, 
		vec4type& flux		
	) const =0;
};

#endif
