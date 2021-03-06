#include "GFM_settingsfile.hpp"
#include "sim_base.hpp"
#include "sim_onefluid.hpp"
#include "sim_twofluid.hpp"
#include "sim_circularexplosiontest.hpp"
#include <iostream>

int main (int argc, char* argv[])
{
	std::cout << "Beginning 2D GFM code.." << std::endl;

	GFM_settingsfile SF;
	SF.read_settings_file(argv[1]);

	std::cout << "Settings file loaded." << std::endl;
	std::cout << "Simulation name is " << SF.basename << std::endl;
	

	std::shared_ptr<sim_base> sim;

	if (SF.sim_type == "onefluid")
	{
		sim = std::make_shared<sim_onefluid>();
	}
	else if (SF.sim_type == "twofluid")
	{
		sim = std::make_shared<sim_twofluid>();
	}
	else if (SF.sim_type == "CEtest")
	{
		sim = std::make_shared<sim_circularexplosiontest>();
	}
	else
	{
		assert(!"Invalid sim type");
	}

	sim->run_sim(SF);

	std::cout << "Simulation complete!" << std::endl;
	
	return 0;
}
