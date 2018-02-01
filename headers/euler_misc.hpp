/*
 *	DESCRIPTION:	Miscellaneous useful functions for the 2D Euler equations
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef MISC_EULER_H
#define MISC_EULER_H

#include "typedefs.hpp"
#include "stiffened_gas_eos.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cassert>


namespace misc {
	

	inline double specific_ie (const vectype& U)
	{
		double rho = U(0);
		double u = U(1) / rho;
		double v = U(2) / rho;
		return (U(3) / rho) - 0.5 * (u*u + v*v);
	}
	
	
	inline bool is_physical_state (double gamma, double pinf, const vectype& U)
	{
		double rho = U(0);
		double u = U(1) / rho;
		double v = U(2) / rho;
		double e = (U(3) / rho) - 0.5 * (u*u + v*v);
		double p = eos::pressure(gamma, pinf, e, rho);
			
		return rho >= 0.0 && e >= 0.0 && p >= 0.0;
	}
	
	
	inline vectype conserved_to_primitives (double gamma, double pinf, const vectype& U)
	{
		vectype prims (4);
		double rho = U(0);
		double u = U(1) / rho;
		double v = U(2) / rho;
		double e = (U(3) / rho) - 0.5 * (u*u + v*v);
		double p = eos::pressure(gamma, pinf, e, rho);
		
		prims(0) = rho;
		prims(1) = u;
		prims(2) = v;
		prims(3) = p;
		
		return prims;
	}
	
	
	inline vectype flux_primitive_var (double gamma, double pinf, const vectype& W)
	{
		vectype flux (4);
		double rho = W(0);
		double u = W(1);
		double v = W(2);
		double p = W(3);
		double e = eos::specific_ie(gamma, pinf, p, rho);
		double E = rho * e + 0.5 * rho * (u*u + v*v);
		
		flux(0) = rho * u;
		flux(1) = rho * u * u + p;
		flux(2) = rho * u * v;
		flux(3) = u * (E + p);
		
		return flux;
	}

}


// TODO: functions below this!





inline vectype flux_conserved_var (const binarySGparams& eosparams, const vectype& U)
{
	/*
	 * The flux of conserved variables in the x-direction
	 */
	 
	vectype flux (6);
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double p = allairemodel::mixture_pressure(eosparams, rho, e, U(5));
	
	flux(0) = U(0) * u;
	flux(1) = U(1) * u;
	flux(2) = U(2) * u + p;
	flux(3) = U(3) * u;
	flux(4) = (U(4) + p) * u;
	flux(5) = 0.0;
	
	return flux;
}





inline vectype primitives_to_conserved (const binarySGparams& eosparams, const vectype& W)
{
	vectype conserved (6);
	double rho = W(0) + W(1);
	double e = allairemodel::mixture_specific_ie(eosparams, rho, W(4), W(5));
	
	conserved(0) = W(0);
	conserved(1) = W(1);
	conserved(2) = rho * W(2);
	conserved(3) = rho * W(3);
	conserved(4) = rho * e + 0.5 * rho * (W(2) * W(2) + W(3) * W(3));
	conserved(5) = W(5);
	
	return conserved;
}

inline vectype primitives_to_conserved (const binarySGparams& eosparams, double rhoz1, double rhoz2, double u, double v, double p, double z)
{
	vectype W (6);
	W << rhoz1, rhoz2, u, v, p, z;
	return primitives_to_conserved(eosparams, W);
}

inline void A_primitive_vars (const binarySGparams& eosparams, const vectype& W, Matrix6d& A)
{
	/*
	 * Use vector of primitive variables (W) to set the value of the
	 * matrix A which is the Jacobian of the system in primitive
	 * variable form.
	 */
	 
	double rho = W(0) + W(1);
	double c = allairemodel::mixture_soundspeed(eosparams, rho, W(4), W(5));
	
	A = W(2) * Eigen::Matrix<double, 6, 6>::Identity();
	A(0,2) = W(0);
	A(1,2) = W(1);
	A(2,4) = 1.0 / rho;
	A(4,2) = rho * c * c;
}




#endif
