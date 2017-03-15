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
			exit(0);
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

		void enqueue_ND_range_kernel(unsigned int global_size, unsigned int local_size);
		
		cl_kernel id;
		std::weak_ptr<Manager> _M_mgr;
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
		std::shared_ptr<MemObj>		create_buffer(cl_mem_flags mem_flags, unsigned int size)
		{
			printf("create buffer %u\n", size);

			int rc;
			std::shared_ptr<MemObj> ret = std::make_shared<MemObj>();
			ret->id = clCreateBuffer(context, mem_flags, size, NULL, &rc);
			ret->_M_mgr = shared_from_this();
			_M_memobj.push_back(ret);
			return ret;
		}
		std::shared_ptr<Program>	create_program(char ** fileNames, int len, std::string args)
		{
			
			char ** source_str = new char*[len];
			//char source_str[MAX_SOURCE_SIZE];
			size_t * source_size = new size_t[len];

			/* Load the source code containing the kernel*/

			for (int i = 0; i < len; ++i)
			{
				FILE *fp;
				printf("opening file %s\n", fileNames[i]);
				fopen_s(&fp, fileNames[i], "r");

				if (!fp) {
					printf("file load error\n");
					throw std::exception("Failed to load kernel");
					abort();
				}
				
				source_str[i] = (char*)malloc(MAX_SOURCE_SIZE);
				source_size[i] = fread(source_str[i], 1, MAX_SOURCE_SIZE, fp);
				fclose(fp);

				source_str[i][source_size[i]] = 0;

				//printf("source size = %u\n", source_size[i]);
				//printf("%s\n", source_str[i]);
				
			}

			cl_int ret;

			cl_program program_id;

			program_id = clCreateProgramWithSource(context, 1, (const char **)source_str, (const size_t *)source_size, &ret);

			errorcheck("clCreateProgramWithSource", ret);

			free(source_str);

			//printf("building program %s\n", fileNames);

			printf("define string: %s\n", args.c_str());

			/* Build Kernel Program */
			ret = clBuildProgram(program_id, 1, &device_id, args.c_str(), NULL, NULL);

			if (ret != CL_SUCCESS) {
				char* programLog;
				cl_build_status status;
				size_t logSize;

				// check build error and build status first
				clGetProgramBuildInfo(program_id, device_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);

				// check build log
				clGetProgramBuildInfo(program_id, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);

				programLog = (char*)calloc(logSize + 1, sizeof(char));

				clGetProgramBuildInfo(program_id, device_id, CL_PROGRAM_BUILD_LOG, logSize + 1, programLog, NULL);

				printf("Build failed; error=%d, status=%d, programLog:nn%s", ret, status, programLog);
				free(programLog);
				getchar(); exit(0);
			}

			

			printf("build successful\n");

			auto program = std::make_shared<Program>();
			program->id = program_id;
			program->_M_mgr = shared_from_this();

			_M_programs.push_back(program);

			return program;
		}
		void flush()
		{
			cl_int ret;
			ret = clFlush(_M_command_queue);
			OCL::errorcheck("clFlush", ret);
			ret = clFinish(_M_command_queue);
			OCL::errorcheck("clFinish", ret);
		}
		void shutdown()
		{
			cl_int ret;
			
			
			
			_M_programs.clear();
			_M_memobj.clear();

			
			ret = clReleaseCommandQueue(_M_command_queue);
			ret = clReleaseContext(context);
		}
		
		

		cl_platform_id platform_id = NULL;
		cl_device_id device_id = NULL;
		cl_context context = NULL;
		cl_command_queue _M_command_queue = NULL;

		std::vector<std::shared_ptr<MemObj>> _M_memobj;
		std::vector<std::shared_ptr<Program>> _M_programs;

		//cl_mem memobj1 = NULL;
		//cl_mem memobj2 = NULL;


		
		//cl_kernel kernel = NULL;
	};

	void OCLtest(OCL::Manager& ocl);

}


#endif