#ifndef PURE_RS_BASE_H
#define PURE_RS_BASE_H


#include "typedefs.hpp"
#include "stiffened_gas_eos.hpp"


class pure_RS_base {

	public:

	virtual void solve_RP (	
		const bool paralleltox
		const binarySGparams& eosparams,
		const vectype& Lstate, 
		const vectype& Rstate, 
		vectype& flux,			
	) const =0;
};

#endif
