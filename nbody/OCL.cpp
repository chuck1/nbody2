#include "stdafx.h"

#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include "OCL.h"
#include "Misc.h"


void OCL::OCLtest(OCL::Manager& ocl)
{
	cl_int ret;

	char string[MEM_SIZE];

	
	// initialize memobj2
	unsigned int i = 0;
	//clEnqueueWriteBuffer(ocl.command_queue, ocl.memobj2, CL_FALSE, 0, sizeof(unsigned int), &i, 0, 0, 0);

	ocl._M_memobj[1]->EnqueueWrite(&i, sizeof(unsigned int));

	unsigned int global_size = 5;
	unsigned int local_size = 5;

	/* Execute OpenCL Kernel */
	//ret = clEnqueueTask(ocl.command_queue, ocl.kernel, 0, NULL, NULL);

	auto kernel = ocl._M_programs[0]->_M_kernels[0];

	ret = clEnqueueNDRangeKernel(ocl._M_command_queue, kernel->id, 1, 0, &global_size, &local_size, 0, 0, 0);
	OCL::errorcheck("clEnqueueNDRangeKernel", ret);

	/* Copy results from the memory buffer */
	ret = clEnqueueReadBuffer(ocl._M_command_queue, ocl._M_memobj[0]->id, CL_TRUE, 0, MEM_SIZE * sizeof(char), string, 0, NULL, NULL);

	/* Display Result */
	puts(string);

	

	//printf("");
	getchar();
}


void								OCL::Manager::init()
{
	cl_int ret;
	cl_uint ret_num_devices;
	cl_uint ret_num_platforms;

	/* Get Platform and Device Info */
	ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);

	


	/* Create OpenCL context */
	context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

	/* Create Command Queue */
	_M_command_queue = clCreateCommandQueue(context, device_id, 0, &ret);



	char device_version[512];
	char device_name[512];
	
	unsigned int device_max_work_item_dimensions;
	unsigned int * device_max_work_item_sizes;
	cl_ulong device_global_mem_size;
	cl_ulong device_local_mem_size;
	cl_ulong device_max_mem_alloc_size;
	cl_uint device_max_compute_units;

	clGetDeviceInfo(device_id, CL_DEVICE_VERSION, 512, device_version, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_NAME, 512, device_name, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(unsigned int), &device_max_work_group_size, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(unsigned int), &device_max_work_item_dimensions, NULL);
	
	device_max_work_item_sizes = new unsigned int[device_max_work_item_dimensions];

	clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(unsigned int)*device_max_work_item_dimensions, device_max_work_item_sizes, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &device_global_mem_size, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_local_mem_size, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &device_max_mem_alloc_size, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &device_max_compute_units, NULL);

	printf("CL_DEVICE_NAME:                %s\n", device_name);
	printf("CL_DEVICE_VERSION:             %s\n", device_version);
	printf("CL_DEVICE_MAX_WORK_GROUP_SIZE: %16u\n", device_max_work_group_size);
	printf("CL_DEVICE_MAX_WORK_ITEM_SIZES: %16u\n", device_max_work_item_sizes[0]);
	printf("CL_DEVICE_GLOBAL_MEM_SIZE:     %16llu\n", device_global_mem_size);
	printf("CL_DEVICE_LOCAL_MEM_SIZE:      %16llu\n", device_local_mem_size);
	printf("CL_DEVICE_MAX_MEM_ALLOC_SIZE:  %16llu\n", device_max_mem_alloc_size);
	printf("CL_DEVICE_MAX_COMPUTE_UNITS:   %16u\n", device_max_compute_units);

}
std::shared_ptr<OCL::MemObj>		OCL::Manager::create_buffer(cl_mem_flags mem_flags, unsigned int size)
{
	printf("create buffer %16u B %16u KB\n", size, size/1024);

	int rc;
	std::shared_ptr<OCL::MemObj> ret = std::make_shared<OCL::MemObj>();

	ret->id = clCreateBuffer(context, mem_flags, size, NULL, &rc);

	ret->_M_mgr = shared_from_this();
	_M_memobj.push_back(ret);
	return ret;
}
std::shared_ptr<OCL::Program>		OCL::Manager::create_program(char ** fileNames, int len, std::string args)
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

		throw std::exception();
	}

	printf("build successful\n");

	auto program = std::make_shared<Program>();
	program->id = program_id;
	program->_M_mgr = shared_from_this();

	_M_programs.push_back(program);

	return program;
}
void								OCL::Manager::test()
{
	init();

	char buf[64*64*8*16];

	TimeMeasurement tm;


	unsigned int S[4];
	unsigned int s = 64*64*8;

	std::shared_ptr<OCL::MemObj> memobj[4];

	for (unsigned int i = 0; i < 4; ++i)
	{
		S[i] = s;
		s *= 2;
		memobj[i] = create_buffer(CL_MEM_READ_WRITE, S[i]);
	}



	for (unsigned int i = 0; i < 4; ++i)
	{
		tm.start();
		memobj[i]->EnqueueRead(buf, S[i]);
		tm.stop();

		printf("transfered %u bytes in %f seconds. %f KiB/s\n", S[i], tm.dif, (double)S[i] / 1024 / tm.dif);
	}

	shutdown();
}
void								OCL::Manager::flush()
{
	cl_int ret;
	ret = clFlush(_M_command_queue);
	OCL::errorcheck("clFlush", ret);
	ret = clFinish(_M_command_queue);
	OCL::errorcheck("clFinish", ret);
}
void								OCL::Manager::shutdown()
{
	printf("OCL shutdown\n");

	cl_int ret;

	_M_programs.clear();
	_M_memobj.clear();

	ret = clReleaseCommandQueue(_M_command_queue);
	ret = clReleaseContext(context);
}


OCL::MemObj::~MemObj()
{
	printf("MemObj dtor\n");
	int ret;
	ret = clReleaseMemObject(id);
}
void	OCL::MemObj::EnqueueWrite(void * src, unsigned int size)
{
	//printf("EnqueueWrite %16x %8i\n", src, size);

	auto mgr = _M_mgr.lock();
	int ret = clEnqueueWriteBuffer(mgr->_M_command_queue, id, CL_FALSE, 0, size, src, 0, 0, 0);
	OCL::errorcheck("clEnqueueWriteBuffer", ret);
}
void	OCL::MemObj::EnqueueRead(void * dst, unsigned int size)
{
	//printf("EnqueueRead  %16x %8i\n", dst, size);

	auto mgr = _M_mgr.lock();
	cl_int ret = clEnqueueReadBuffer(mgr->_M_command_queue, id, CL_TRUE, 0, size, dst, 0, NULL, NULL);
	OCL::errorcheck("clEnqueueReadBuffer", ret);
}



OCL::Kernel::~Kernel()
{
	cl_int ret;
	ret = clReleaseKernel(id);
}
void	OCL::Kernel::set_arg(std::shared_ptr<MemObj> memobj, int arg)
{
	cl_int ret;
	ret = clSetKernelArg(id, arg, sizeof(cl_mem), (void *)&(memobj->id));

	char s[128];
	sprintf_s(s, "clSetKernelArg %i", arg);
	errorcheck(s, ret);
}
void	OCL::Kernel::set_arg(int arg, unsigned int size)
{
	cl_int ret;
	ret = clSetKernelArg(id, arg, size, NULL);

	char s[128];
	sprintf_s(s, "clSetKernelArg %i", arg);
	errorcheck(s, ret);
}
void	OCL::Kernel::enqueue_ND_range_kernel(unsigned int global_size, unsigned int local_size)
{
	std::shared_ptr<Manager> mgr = _M_mgr.lock();
	cl_int ret = clEnqueueNDRangeKernel(mgr->_M_command_queue, id, 1, 0, &global_size, &local_size, 0, 0, 0);

	char buffer[128];
	sprintf_s(buffer, "clEnqueueNDRangeKernel %s", name.c_str());
	OCL::errorcheck(buffer, ret);
}
void	OCL::Kernel::enqueue_ND_range_kernel()
{
	enqueue_ND_range_kernel(gs, ls);
}




std::shared_ptr<OCL::Kernel>		OCL::Program::create_kernel(char const * kernel_name)
{
	cl_int ret;

	cl_kernel kernel_id = clCreateKernel(id, kernel_name, &ret);

	char buffer[512];
	sprintf_s(buffer, "clCreateKernel %s", kernel_name);
	errorcheck(buffer, ret);

	auto kernel = std::make_shared<Kernel>();

	kernel->id = kernel_id;
	kernel->_M_mgr = _M_mgr;
	kernel->name = kernel_name;

	_M_kernels.push_back(kernel);



	return kernel;
}



//void								OCL::RoutineArrayReduce::init(unsigned int len)
//{
//	auto m = mgr.lock();
//
//	m->create_buffer(CL_MEM_READ_WRITE, len * sizeof(double));
//}
//void								OCL::RoutineArrayReduce::write(double * arr)
//{
//
//}
//void								OCL::RoutineArrayReduce::run()
//{
//
//}
