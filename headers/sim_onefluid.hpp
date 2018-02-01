#ifndef ONEFLUID_SIM_H
#define ONEFLUID_SIM_H


#include "sim_info.hpp"
#include "typedefs.hpp"
#include "sim_base.hpp"
#include "GFM_settingsfile.hpp"
#include "stiffened_gas_eos.hpp"
#include "pure_RS_base.hpp"
#include "FS_base.hpp"
#include <memory>



class sim_onefluid : public sim_base {
	
	private:
	
	void set_sim_parameters (const GFM_settingsfile& SF, sim_info& params, binarySGparams& eosparams);
	
	void set_sim_methods (const GFM_settingsfile& SF, std::shared_ptr<pure_RS_base>& RS, std::shared_ptr<FS_base>& FS);
	
	void set_sim_ICs (const GFM_settingsfile& SF, const sim_info& params, const binarySGparams& eosparams, grideuler2type& grid);
	
	void output (int numsteps, double t, const sim_info& params, const binarySGparams& eosparams, const grideuler2type& grid, std::string filename);
	
	void set_planar_IC (const vec4type& U_under, const vec4type& U_over, const double a, const double b, const double c, grideuler2type& grid, const sim_info& params);
	

	public:

	void run_sim (GFM_settingsfile SF);

	double compute_dt (double CFL, const sim_info& params, const grideuler2type& grid, const binarySGparams& eos, double T, double t);
};

#endif
