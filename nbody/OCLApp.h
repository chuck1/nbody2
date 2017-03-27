#ifndef OCLAPP_H
#define OCLAPP_H

#include <vector>
#include <memory>

#include "decl.h"

class OCLApp
{
public:
	void	init(
		std::vector<Body0> & bodies0,
		std::vector<Body1> & bodies1,
		std::vector<Pair> & pairs,
		Header & header)
	{
		int n = bodies0.size();
		int p = pairs.size();

		ocl = std::make_shared<OCL::Manager>();
		ocl->init();

		char * files[2] = { "kernel.cl", "kernel.h" };

		std::shared_ptr<OCL::Program> program;

		while (true)
		{
			try
			{
				program = ocl->create_program(files, 1, "");
				break;
			}
			catch (...)
			{
				printf("try again\n");
				getchar();
			}
		}

		kernel_calc_acc = program->create_kernel("calc_acc");
		kernel_step = program->create_kernel("step_pos");
		kernel_dt_calc = program->create_kernel("dt_calc");
		kernel_dt_store = program->create_kernel("store_dt");

		std::cout << "bodies0" << std::endl;
		memobj_bodies0 = ocl->create_buffer(CL_MEM_READ_WRITE, bodies0.size() * sizeof(Body0));
		memobj_bodies0->EnqueueWrite(&bodies0[0], bodies0.size() * sizeof(Body0));

		std::cout << "bodies1" << std::endl;
		memobj_bodies1 = ocl->create_buffer(CL_MEM_READ_WRITE, bodies1.size() * sizeof(Body1));
		memobj_bodies1->EnqueueWrite(&bodies1[0], bodies1.size() * sizeof(Body1));

		std::cout << "pairs" << std::endl;
		memobj_pairs = ocl->create_buffer(CL_MEM_READ_WRITE, pairs.size() * sizeof(Pair));
		memobj_pairs->EnqueueWrite(&pairs[0], pairs.size() * sizeof(Pair));

		memobj_header = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(Header));
		memobj_header->EnqueueWrite(&header, sizeof(Header));

		/*unsigned int counter = 0;
		auto memobj_counter = ocl->create_buffer(CL_MEM_READ_WRITE, sizeof(unsigned int));
		memobj_counter->EnqueueWrite(&counter, sizeof(unsigned int));*/

		auto memobj_dt_input = ocl->create_buffer(CL_MEM_READ_WRITE, (n + p) * sizeof(double));


		// set kernel args

		kernel_calc_acc->set_arg(memobj_header, 0);
		kernel_calc_acc->set_arg(memobj_bodies0, 1);
		kernel_calc_acc->set_arg(memobj_bodies1, 2);
		kernel_calc_acc->set_arg(memobj_pairs, 3);
		kernel_calc_acc->set_arg(memobj_dt_input, 4);
		//kernel_calc_acc->set_arg(memobj_counter, 5);

		kernel_step->set_arg(memobj_header, 0);
		kernel_step->set_arg(memobj_bodies0, 1);
		kernel_step->set_arg(memobj_bodies1, 2);
		kernel_step->set_arg(memobj_pairs, 3);

		kernel_dt_calc->set_arg(memobj_header, 0);
		kernel_dt_calc->set_arg(memobj_bodies0, 1);
		kernel_dt_calc->set_arg(memobj_bodies1, 2);
		kernel_dt_calc->set_arg(memobj_dt_input, 3);
		

		routine_min_dt.reset(new OCL::RoutineArrayReduce<double>);
		routine_min_dt->_M_mgr = ocl;
		routine_min_dt->_M_prg = program;
		routine_min_dt->_M_memobj_in = memobj_dt_input;
		routine_min_dt->init("k_min_double", n + p, ocl->device_max_work_group_size);


		kernel_dt_store->set_arg(routine_min_dt->_M_memobj_out.lock(), 0);
		kernel_dt_store->set_arg(memobj_header, 1);

#if 1
		// test
		const unsigned int n1 = 35;
		unsigned int arr[n1];
		for (unsigned int i = 0; i < n1; ++i)
		{
			arr[i] = i+1;
		}

		OCL::RoutineArrayReduce<unsigned int> routine;
		routine._M_mgr = ocl;
		routine._M_prg = program;
		routine.init("k_sum_uint", n1, 4);
		routine.write(arr);

		routine.run();

		unsigned int res;

		routine._M_memobj_out.lock()->EnqueueRead(&res, sizeof(unsigned int));

		/*unsigned long arr0[n1];
		unsigned long arr1[1];

		routine._M_memobj0.lock()->EnqueueRead(&arr0, n1 * sizeof(unsigned long));

		routine._M_memobj1.lock()->EnqueueRead(&arr1, sizeof(unsigned long));*/

		printf("sum of 1 ... %u = %lu\n", n1, res);
#endif

		
	}
	void	shutdown()
	{
		memobj_bodies0.reset();
		memobj_bodies1.reset();
		memobj_pairs.reset();
		memobj_header.reset();

		kernel_calc_acc.reset();
		kernel_step.reset();
		kernel_dt_calc.reset();
		kernel_dt_store.reset();

		routine_min_dt.reset();

		ocl->shutdown();
	}

	std::shared_ptr<OCL::Manager>	ocl;

	std::shared_ptr<OCL::MemObj>	memobj_bodies0;
	std::shared_ptr<OCL::MemObj>	memobj_bodies1;
	std::shared_ptr<OCL::MemObj>	memobj_pairs;
	std::shared_ptr<OCL::MemObj>	memobj_header;

	std::shared_ptr<OCL::Kernel>	kernel_calc_acc;
	std::shared_ptr<OCL::Kernel>	kernel_step;
	std::shared_ptr<OCL::Kernel>	kernel_dt_calc;
	std::shared_ptr<OCL::Kernel>	kernel_dt_store;

	std::shared_ptr<OCL::RoutineArrayReduce<double>>	routine_min_dt;
};

#endif