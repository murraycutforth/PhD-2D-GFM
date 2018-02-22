#include "sim_twofluid.hpp"

class sim_circularexplosiontest : public sim_twofluid {
	
	private:
	
	double compute_onefluid_dt (double CFL, const sim_info& params, const grideuler2type& grid, const binarySGparams& eosparams, double T, double t);

	void compute_errors (double& l1, double& l2, double& linf, const sim_info& params1, const grideuler2type& gridref, const sim_info& params, const grideuler2type& grid1, const grideuler2type& grid2, const GFM_ITM_interface& ls);
	
	
	public:
	
	virtual void run_sim (GFM_settingsfile SF) override;

};
