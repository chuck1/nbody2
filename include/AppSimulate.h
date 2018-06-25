#ifndef APPSIMULATE_H
#define APPSIMULATE_H

#include <vector>

#include <boost/program_options.hpp>

#include "kernel.h"
#include "App.h"
#include "decl.h"
#include "History.h"

class AppSimulate : public App
{
public:
	AppSimulate();
	virtual ~AppSimulate();
	virtual void		run(std::vector<std::string> vec);
	int					simulate();
	void				generate_binary();
	void				generate_cluster();
	void				generate_disc();
	void				generate_read();
	void				generate_pairs();
	void				gen_func();

	void				print_simulation_step(int i);

	Header				header;
	std::vector<Body0>	bodies0;
	std::vector<Body1>	bodies1;
	std::vector<Pair>	pairs;

	History				hist;

	std::string			string_gen_func;

	boost::program_options::variables_map	vm;

	int					m;

	OCLApp *			ocl_app;

	struct
	{
		int n;
		std::string file;
		double velocity_factor;
	} generator_options;

	// profiling
	TimeMeasurement t_kernel, t_readbuffer;
};

#endif