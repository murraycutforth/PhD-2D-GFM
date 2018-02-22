#ifndef TWOFLUID_SIM_H
#define TWOFLUID_SIM_H


#include "sim_info.hpp"
#include "typedefs.hpp"
#include "sim_base.hpp"
#include "GFM_settingsfile.hpp"
#include "GFM_base.hpp"
#include "GFM_ITM_interface.hpp"
#include "stiffened_gas_eos.hpp"
#include "pure_RS_base.hpp"
#include "FS_base.hpp"
#include "levelset.hpp"
#include <memory>


class sim_twofluid : public sim_base {
	
	protected:
	
	void set_sim_parameters 
	(
		const GFM_settingsfile& SF, 
		sim_info& params, 
		binarySGparams& eosparams
	);
	
	void set_sim_methods 
	(
		const GFM_settingsfile& SF, 
		const sim_info& params,
		std::shared_ptr<FS_base>& FS,
		std::shared_ptr<GFM_base>& GFM,
		std::shared_ptr<GFM_ITM_interface>& ls
	);
	
	void set_sim_ICs 
	(
		const GFM_settingsfile& SF, 
		const sim_info& params, 
		const binarySGparams& eosparams, 
		grideuler2type& grid1, 
		grideuler2type& grid2, 
		GFM_ITM_interface& ls
	);
	
	void output 
	(
		int numsteps, 
		double t, 
		const sim_info& params, 
		const binarySGparams& eosparams,
		const grideuler2type& grid1, 
		const grideuler2type& grid2, 
		const GFM_ITM_interface& ls, 
		std::string filename,
		const velocity_field_base& vfield
	);
	
	void output_mass_error
	(
		int numsteps,
		double t,
		const sim_info& params,
		const grideuler2type& grid1, 
		const grideuler2type& grid2, 
		const GFM_ITM_interface& ls, 
		std::string filename
	);


	double fluid_mass
	(
		const double interiorsign,
		const sim_info& params,
		const grideuler2type& grid,
		const GFM_ITM_interface& ls
	);
		

	public:

	virtual void run_sim (GFM_settingsfile SF) override;

	double compute_dt (double CFL, const sim_info& params, const grideuler2type& grid1, const grideuler2type& grid2, const GFM_ITM_interface& ls, const binarySGparams& eosparams, double T, double t);
};



#endif
