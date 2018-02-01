#ifndef GFM_SF_H
#define GFM_SF_H

#include <string>
#include <fstream>
#include <sstream>
#include <cassert>

class GFM_settingsfile {
	
	public:
	
	int Nx;
	int Ny;
	double CFL;
	std::string sim_type;
	std::string method;
	std::string test_case;
	std::string outputpath;
	std::string basename;
	
	void read_settings_file (std::string filename)
	{
		std::ifstream infile(filename);
		std::string line;
	
		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			std::string inputname;
	
			iss >> inputname;
			
			if (inputname == "Nx") iss >> Nx;
			
			else if (inputname == "Ny") iss >> Ny;
			
			else if (inputname == "CFL") iss >> CFL;
			
			else if (inputname == "sim_type") iss >> sim_type;
			
			else if (inputname == "method") iss >> method;
			
			else if (inputname == "test_case") iss >> test_case;
						
			else if (inputname == "outputpath") iss >> outputpath;
						
		}
			
		basename = outputpath + "-" + sim_type + "-" + test_case + "-" + method + "-" + std::to_string(Nx) + "-" + std::to_string(Ny);
		infile.close();
	}
};

#endif
