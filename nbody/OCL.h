#ifndef OCL_H
#define OCL_H

#include <exception>
#include <vector>
#include <memory>
#include <iostream>

#include <CL/cl.h>

#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x10000)

namespace OCL
{
	class Error: public std::exception
	{

	};

	inline const char *getErrorString(cl_int error)
	{
		switch (error){
		case CL_OUT_OF_RESOURCES:
			return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:
			return "CL_OUT_OF_HOST_MEMORY";


			// run-time and JIT compiler errors
		case 0: return "CL_SUCCESS";
		case -1: return "CL_DEVICE_NOT_FOUND";
		case -2: return "CL_DEVICE_NOT_AVAILABLE";
		case -3: return "CL_COMPILER_NOT_AVAILABLE";
		case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		//case -5: return "CL_OUT_OF_RESOURCES";
		//case -6: return "CL_OUT_OF_HOST_MEMORY";
		case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case -8: return "CL_MEM_COPY_OVERLAP";
		case -9: return "CL_IMAGE_FORMAT_MISMATCH";
		case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case -11: return "CL_BUILD_PROGRAM_FAILURE";
		case -12: return "CL_MAP_FAILURE";
		case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case -15: return "CL_COMPILE_PROGRAM_FAILURE";
		case -16: return "CL_LINKER_NOT_AVAILABLE";
		case -17: return "CL_LINK_PROGRAM_FAILURE";
		case -18: return "CL_DEVICE_PARTITION_FAILED";
		case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

			// compile-time errors
		case -30: return "CL_INVALID_VALUE";
		case -31: return "CL_INVALID_DEVICE_TYPE";
		case -32: return "CL_INVALID_PLATFORM";
		case -33: return "CL_INVALID_DEVICE";
		case -34: return "CL_INVALID_CONTEXT";
		case -35: return "CL_INVALID_QUEUE_PROPERTIES";
		case -36: return "CL_INVALID_COMMAND_QUEUE";
		case -37: return "CL_INVALID_HOST_PTR";
		case -38: return "CL_INVALID_MEM_OBJECT";
		case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case -40: return "CL_INVALID_IMAGE_SIZE";
		case -41: return "CL_INVALID_SAMPLER";
		case -42: return "CL_INVALID_BINARY";
		case -43: return "CL_INVALID_BUILD_OPTIONS";
		case -44: return "CL_INVALID_PROGRAM";
		case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
		case -46: return "CL_INVALID_KERNEL_NAME";
		case -47: return "CL_INVALID_KERNEL_DEFINITION";
		case -48: return "CL_INVALID_KERNEL";
		case -49: return "CL_INVALID_ARG_INDEX";
		case -50: return "CL_INVALID_ARG_VALUE";
		case -51: return "CL_INVALID_ARG_SIZE";
		case -52: return "CL_INVALID_KERNEL_ARGS";
		case -53: return "CL_INVALID_WORK_DIMENSION";
		case -54: return "CL_INVALID_WORK_GROUP_SIZE";
		case -55: return "CL_INVALID_WORK_ITEM_SIZE";
		case -56: return "CL_INVALID_GLOBAL_OFFSET";
		case -57: return "CL_INVALID_EVENT_WAIT_LIST";
		case -58: return "CL_INVALID_EVENT";
		case -59: return "CL_INVALID_OPERATION";
		case -60: return "CL_INVALID_GL_OBJECT";
		case -61: return "CL_INVALID_BUFFER_SIZE";
		case -62: return "CL_INVALID_MIP_LEVEL";
		case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
		case -64: return "CL_INVALID_PROPERTY";
		case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
		case -66: return "CL_INVALID_COMPILER_OPTIONS";
		case -67: return "CL_INVALID_LINKER_OPTIONS";
		case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

			// extension errors
		case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
		case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
		case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
		case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
		case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
		case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
		default: return "Unknown OpenCL error";
		}
	}
	inline void errorcheck(const char * s, int result)
	{
		if (result != CL_SUCCESS){
			std::cerr << s << ":" << getErrorString(result) << std::endl;
			getchar();
			throw OCL::Error();
		}
	}

	class Manager;

	class MemObj
	{
	public:
		~MemObj();
		void EnqueueWrite(void * src, unsigned int size);
		void EnqueueRead(void * dest, unsigned int sz);

		std::weak_ptr<Manager> _M_mgr;
		cl_mem id;
	};

	class Kernel
	{
	public:
		~Kernel();
		
		void	set_arg(std::shared_ptr<MemObj> memobj, int arg);
		void	set_arg(int arg, unsigned int size);

		void	enqueue_ND_range_kernel(unsigned int global_size, unsigned int local_size);
		void	enqueue_ND_range_kernel();
		
		unsigned int	gs;
		unsigned int	ls;

		cl_kernel id;
		std::weak_ptr<Manager> _M_mgr;
		std::string name;
	};
	
	class Program
	{
	public:
		~Program()
		{
			_M_kernels.clear();

			cl_int ret;
			ret = clReleaseProgram(id);
		}
		std::shared_ptr<Kernel>		create_kernel(char const * kernel_name);

		std::vector<std::shared_ptr<Kernel>>		_M_kernels;

		cl_program id = NULL;
		std::weak_ptr<Manager> _M_mgr;
	};

	class Manager: public std::enable_shared_from_this<Manager>
	{
	public:
		void init();
		std::shared_ptr<MemObj>		create_buffer(cl_mem_flags mem_flags, unsigned int size);
		std::shared_ptr<Program>	create_program(char ** fileNames, int len, std::string args);

		void						flush();
		void						shutdown();
		
		void						test();

		cl_platform_id platform_id = NULL;
		cl_device_id device_id = NULL;
		cl_context context = NULL;
		cl_command_queue _M_command_queue = NULL;

		std::vector<std::shared_ptr<MemObj>> _M_memobj;
		std::vector<std::shared_ptr<Program>> _M_programs;

		unsigned int				device_max_work_group_size;
	};

	template<typename T>
	class RoutineArrayReduce
	{
	public:
		/*void								init(unsigned int len);
		void								write(double * arr);
		void								run();*/

		void								init(const char * s, unsigned int len, unsigned int N)
		{
			_M_len = len;

			std::shared_ptr<OCL::Manager> mgr = _M_mgr.lock();
			std::shared_ptr<OCL::Program> prg = _M_prg.lock();

			unsigned int l = len;

			unsigned int gs = next_power_of_two(l);
			unsigned int ls = std::min(gs, N);
			unsigned int gc = gs / ls;

			std::shared_ptr<OCL::MemObj> m0 = mgr->create_buffer(CL_MEM_READ_WRITE, l * sizeof(T));
			std::shared_ptr<OCL::MemObj> m1 = mgr->create_buffer(CL_MEM_READ_WRITE, gc * sizeof(T));

			_M_memobj0 = m0;
			_M_memobj1 = m1;

			do
			{
				printf("gs=%8u ls=%8u gc=%8u l=%8u\n", gs, ls, gc, l);

				auto k = prg->create_kernel(s);

				k->gs = gs;
				k->ls = ls;

				std::shared_ptr<OCL::MemObj> m = mgr->create_buffer(CL_MEM_READ_WRITE, sizeof(unsigned int));
				m->EnqueueWrite(&l, sizeof(unsigned int));

				k->set_arg(m0, 0);
				k->set_arg(m, 1);
				k->set_arg(m1, 2);
				k->set_arg(3, N * sizeof(T));

				_M_kernels.push_back(k);

				std::swap(m0, m1);

				l = gc;
				gs = gc;
				ls = std::min(gs, N);
				gc = gs / ls;
			} while (gs > 1);
			
			_M_memobj_out = m0;
		}
		void								write(T * arr)
		{
			std::shared_ptr<OCL::MemObj> m0 = _M_memobj0.lock();

			m0->EnqueueWrite(arr, _M_len);
		}
		void								run()
		{
			for (unsigned int i = 0; i < _M_kernels.size(); ++i)
			{
				std::shared_ptr<OCL::Kernel> k = _M_kernels[i].lock();
				k->enqueue_ND_range_kernel();
			}
		}

		unsigned int						_M_len;

		std::weak_ptr<Manager>				_M_mgr;
		std::weak_ptr<Program>				_M_prg;

		std::vector<std::weak_ptr<Kernel>>	_M_kernels;
		std::weak_ptr<MemObj>				_M_memobj0;
		std::weak_ptr<MemObj>				_M_memobj1;
		std::weak_ptr<MemObj>				_M_memobj_out;
	};

	void OCLtest(OCL::Manager& ocl);

}


#endif