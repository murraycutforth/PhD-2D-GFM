#include "sim_twofluid.hpp"

class sim_circularexplosiontest : public sim_twofluid {
	
	private:
	
	double compute_onefluid_dt (double CFL, const sim_info& params, const grideuler2type& grid, const binarySGparams& eosparams, double T, double t);
	
	
	public:
	
	virtual void run_sim (GFM_settingsfile SF) override;

};
