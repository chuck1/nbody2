#include "stdafx.h"

#include "OCL.h"
#include "AppSimulate.h"
#include "OCLApp.h"
#include "Misc.h"
#include "mymath.h"

void clear_bodies(std::vector<Body0> & bodies0)
{
	for (unsigned int i = 0; i < bodies0.size(); ++i)
	{
		bodies0[i].q.v[0] = 1;
		bodies0[i].q.v[1] = 0;
		bodies0[i].q.v[2] = 0;
		bodies0[i].q.v[3] = 0;
	}
}

AppSimulate::AppSimulate()
{
	generator_options.velocity_factor = 1.0;

	ocl_app = 0;

	m = 0;
}
AppSimulate::~AppSimulate()
{
	if (ocl_app) ocl_app->shutdown();
}
void				AppSimulate::run(std::vector<std::string> vec)
{
	boost::program_options::options_description desc("simulation options");

	desc.add_options()
		("help", "produce help message")
		("m", boost::program_options::value<int>(), "simulation steps")
		("n", boost::program_options::value<int>(), "number of bodies")
		("gen", boost::program_options::value<std::string>(), "generator function")
		("file", boost::program_options::value<std::string>(), "filename of simulation to read")
		("vf", boost::program_options::value<double>(), "velocity factor for generator function");

	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(vec).options(desc).run(), vm);
	}
	catch (std::exception & e)
	{
		std::cout << e.what() << std::endl;
		std::cout << desc << std::endl;
		return;
	}

	boost::program_options::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return;
	}

	if (vm.count("m"))
	{
		m = vm["m"].as<int>();
	}

	if (vm.count("n"))
	{
		generator_options.n = vm["n"].as<int>();
	}

	if (vm.count("gen"))
	{
		string_gen_func = vm["gen"].as<std::string>();
	}

	if (vm.count("vf"))
	{
		generator_options.velocity_factor = vm["vf"].as<double>();
	}

	simulate();
}
int					AppSimulate::simulate()
{


	srand(time(NULL));

	hist.folder = "C:\\test";

	header.t = 0;
	header.dt = 0;

	gen_func();

	generate_pairs();

	int n = bodies0.size();
	int p = n*(n - 1) / 2;

	hist.resize(n);

	unsigned int file_size = m * (n * sizeof(Body0)+p * sizeof(Pair));

	printf("file size = %u MB\n", file_size / 1024 / 1024);

	double * dt_partial = new double[p];

	ocl_app = new OCLApp();
	ocl_app->init(bodies0, bodies1, pairs, header);


	std::cout << "kernel start" << std::endl;

	// if this is a new simulation, save the zero frame
	if (hist.frame_times.empty())
	{
		hist.push(ocl_app->memobj_header, ocl_app->memobj_bodies0, ocl_app->memobj_pairs, bodies0.size());
	}
	else
	{
		// think this is done above already
		//header.t = hist.frame_times.back();
	}

	int dt_len = n + p;

	double * dt_input = new double[n + p];

	for (int i = 0; i < m; ++i)
	{
		double t0 = header.t;

		//printf("kernel execution i=%i\n", i);

		t_kernel.start();

		ocl_app->kernel_calc_acc->enqueue_ND_range_kernel(1024, 1024);

		ocl_app->kernel_dt_calc->enqueue_ND_range_kernel(next_power_of_two(n), 1);

		
		ocl_app->routine_min_dt->run();


		ocl_app->kernel_dt_store->enqueue_ND_range_kernel(1, 1);

		ocl_app->kernel_step->enqueue_ND_range_kernel(next_power_of_two(n), 1);

		t_kernel.stop();

		if (should_print_step(i,m) || should_save_frame(i))
		{
			//
		}

		// doing this at the end of step kernel
		//header.t = t0 + header.dt;

		

		// only save every tenth frame
		if (should_save_frame(i)) {

			t_readbuffer.start();

			//om.memobj_header->EnqueueRead(&header, sizeof(Header));
			hist.push(ocl_app->memobj_header, ocl_app->memobj_bodies0, ocl_app->memobj_pairs, n);

			t_readbuffer.stop();

			//printf("%8i %16.2e %16.2e\n", i, t_kernel.dif, t_readbuffer.dif);

		}

		if (should_print_step(i, m)) {
			ocl_app->memobj_header->EnqueueRead(&header, sizeof(Header));
			print_simulation_step(i);
		}

#if 0
		unsigned int ret;
		ocl_app->memobj_ret->EnqueueRead(&ret, sizeof(unsigned int));
		printf("ret = %8u p = %8u (n+p-1) = %8u\n", ret, p, n+p-1);


		ocl_app->memobj_dt_input->EnqueueRead(dt_input, sizeof(double)* (n + p));
		for (int i = 0; i < (n + p); ++i){
			if (dt_input[i] < 1.0E0)
				printf("  %8i %16.2e\n", i, dt_input[i]);
		}
#endif
		// at n=64, there is a bug, dt comes out 0 for all time steps

		if (signaled == 1)
		{
			init_signal();
			signaled = 0;
			break;
		}
	}


	hist.bodies1.resize(n);
	ocl_app->memobj_bodies1->EnqueueRead(&hist.bodies1[0], n * sizeof(Body1));
	hist.write();


	std::cout << "done" << std::endl;

	/*char message[256];
	sprintf_s(message, "simulation complete");
	client("192.168.56.2","4001",message);*/

	return 0;
}
void				AppSimulate::generate_binary()
{
	bodies0.resize(2);
	bodies1.resize(2);

	clear_bodies(bodies0);

	Body0 & b0_0 = bodies0[0];
	Body1 & b1_0 = bodies1[0];

	Body0 & b0_1 = bodies0[1];
	Body1 & b1_1 = bodies1[1];

	double & x0 = b0_0.pos.v[0];
	double & x1 = b0_1.pos.v[0];

	double & m0 = b1_0.mass;
	double & m1 = b1_1.mass;

	b1_0.density = 5000.0;
	b1_1.density = 5000.0;

	m0 = 1.0E11;
	m1 = 1.0E11;

	b0_0.radius = pow(3.0 * b1_0.mass / 4.0 / M_PI / b1_0.density, 1.0 / 3.0);
	b0_1.radius = pow(3.0 * b1_1.mass / 4.0 / M_PI / b1_1.density, 1.0 / 3.0);

	printf("b0.mass   = %16.2e\n", b1_0.mass);
	printf("b0.radius = %16.2e\n", b0_0.radius);
	printf("b1.mass   = %16.2e\n", b1_1.mass);
	printf("b1.radius = %16.2e\n", b0_1.radius);

	double d = 500.0;

	x0 = d * m1 / (m0 + m1);
	x1 = x0 - d;

	if (1)
	{
		b1_0.vel.v[1] = generator_options.velocity_factor * sqrt(M_G * m1 * x0 / d / d);
		b1_1.vel.v[1] = -generator_options.velocity_factor * sqrt(M_G * m1 * -x1 / d / d);
	}


	//b0.w.v[2] = 0.001;

	//b0.q.v[3] = 0.3;

	/*printf("position\n");
	print(b0_0.pos);
	print(b0_1.pos);
	printf("velocity\n");
	print(b1_0.vel);
	print(b1_1.vel);*/

	double F = M_G * m0 * m1 / pow(b0_0.radius + b0_1.radius, 2.0);
	printf("force when just touching = %e\n", F);

	double V = mymath::volume_sphere_sphere_intersection(0.9 * (b0_0.radius + b0_1.radius), b0_0.radius, b0_1.radius);
	printf("intersection volume at 10%% pen = %e\n", V);
	double m = V * (b1_0.density + b1_1.density);
	printf("total mass f that region = %e\n", m);
}
void				AppSimulate::generate_cluster()
{
	int n = 3;

	bodies0.resize(n);
	bodies1.resize(n);

	double mass = 1.0E11;
	double density = 5000.0;
	double radius = mymath::sphere_radius(mass / density);

	clear_bodies(bodies0);



	for (int i = 0; i < n; ++i)
	{
		Body0 & b0 = bodies0[i];
		Body1 & b1 = bodies1[i];

		b0.radius = radius;
		b1.mass = mass;
		b1.density = density;

		double & x = b0.pos.v[0];

		x = ((double)i - (double)(n - 1) / 2.0) * radius * 1.0;
	}
}
void				AppSimulate::generate_disc()
{
	int n = generator_options.n;

	bodies0.resize(n);
	bodies1.resize(n);

	clear_bodies(bodies0);

	double x, y, r;

	double body_mass = 1E8;
	double body_density = 5000.0;
	double body_radius = mymath::sphere_radius(body_mass / body_density);

	double radius = sqrt(10.0 * (double)n) * body_radius;

	double density = body_mass * (double)n / (M_PI * pow(radius, 2.0));

	for (int i = 0; i < n; ++i)
	{
		while (true)
		{
			x = (float)rand() / (float)RAND_MAX * 2.0 - 1.0;
			y = (float)rand() / (float)RAND_MAX * 2.0 - 1.0;
			x = x * radius;
			y = y * radius;
			r = sqrt(x*x + y*y);

			if (r < radius) break;
		}

		Body0 & b0 = bodies0[i];
		Body1 & b1 = bodies1[i];

		b1.mass = body_mass;
		b0.radius = body_radius;
		b1.density = body_density;

		b0.pos.v[0] = x;
		b0.pos.v[1] = y;

		//printf("surface gravity = %16.2e", body_surface_gravity(b0, b1));
	}

	double X = 0;
	double Y = 0;
	double Z = 0;
	double M = 0;
	for (int i = 0; i < n; ++i)
	{
		Body0 & b0 = bodies0[i];
		Body1 & b1 = bodies1[i];
		X += b1.mass * b0.pos.v[0];
		Y += b1.mass * b0.pos.v[1];
		Z += b1.mass * b0.pos.v[2];
		M += b1.mass;
	}
	X /= M;
	Y /= M;
	Z /= M;

	for (int i = 0; i < n; ++i)
	{
		Body0 & b0 = bodies0[i];
		Body1 & b1 = bodies1[i];

		double x = b0.pos.v[0] - X;
		double y = b0.pos.v[1] - Y;

		double r = sqrt(x*x + y*y);

		if (r == 0)
		{
			b1.vel.v[1] = 0;
			b1.vel.v[0] = 0;
		}
		else
		{
			double acc2 = mymath::gravity_uniform_disc(density, radius, r);

			//printf("a = %16e\n", acc2);

			double s = sqrt(acc2 * r);

			s *= generator_options.velocity_factor;

			b1.vel.v[1] = s * x / r;
			b1.vel.v[0] = -s * y / r;
		}
	}

}
void				AppSimulate::generate_read()
{
	if (vm.count("file"))
	{
		generator_options.file = vm["file"].as<std::string>();
	}

	hist.folder = generator_options.file;

	hist.load();

	unsigned int i = hist.frame_times.size() - 1;

	auto f = hist.get_frame(i);

	header = f->header;
	bodies0 = f->bodies0;
	bodies1 = hist.bodies1;
}
void				AppSimulate::generate_pairs()
{
	int n = bodies0.size();

	pairs.resize(n*(n - 1) / 2);

	int k = 0;
	for (int i = 0; i < (n - 1); ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			pairs[k].i = i;
			pairs[k].j = j;
			++k;
		}
	}
}
void				AppSimulate::gen_func()
{
	std::map < std::string, std::function<void()> > func_map;

	func_map["binary"] = std::bind(&AppSimulate::generate_binary, this);
	func_map["cluster"] = std::bind(&AppSimulate::generate_cluster, this);
	func_map["disc"] = std::bind(&AppSimulate::generate_disc, this);
	func_map["read"] = std::bind(&AppSimulate::generate_read, this);

	auto it = func_map.find(string_gen_func);

	if (it == func_map.end())
	{
		throw std::exception("invalid generator function");
	}

	// run generator
	it->second();


	generate_pairs();

	header.bodies_size = bodies0.size();
}


void				AppSimulate::print_simulation_step(int i)
{
	printf("%12i t=%10.2e dt=%10.2e count_pen=%8i read speed=%10.2e\n", i, header.t, header.dt, header.count_pen, (double)(bodies0.size()*sizeof(Body0))/t_readbuffer.dif);
}

