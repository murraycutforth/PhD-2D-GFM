#ifndef ONEFLUID_SIM_H
#define ONEFLUID_SIM_H


#include "sim_info.hpp"
#include "typedefs.hpp"
#include "sim_base.hpp"
#include "GFM_settingsfile.hpp"
#include "stiffened_gas_eos.hpp"


class sim_onefluid : public sim_base {
	
	private:
	
	void set_sim_parameters (const GFM_settingsfile& SF, sim_info& params, binarySGparams& eosparams);
	
	void set_sim_methods (const GFM_settingsfile& SF, std::shared_ptr<pure_RS_base>& RS, std::shared_ptr<flow_solver_base>& FS);
	
	void set_sim_ICs (const GFM_settingsfile& SF, const sim_info& params, const binarySGparams& eosparams, gridtype& grid);
	
	void output ();
	

	public:

	void run_sim (GFM_settingsfile SF);

	double compute_dt (double CFL, const sim_info& params, const gridtype& grid, const binarySGparams& eos, double T, double t);
};

#endif
