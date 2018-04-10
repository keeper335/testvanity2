// ocltest3.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

void printPlatformInfo(cl_platform_id platform_id) {
	char nbuf[256];
	size_t len;
	clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, 256, nbuf, &len);
	nbuf[len] = '\0';
	printf("%s\n", nbuf);
}

static void printDeviceInfo(cl_device_id device_id) {
	char name[256], vendor[256];
	cl_device_type type;
	cl_uint comp_units, freq;
	cl_ulong mem_size;
	size_t len;
	clGetDeviceInfo(device_id, CL_DEVICE_NAME, 256, name, &len);
	name[len] = '\0';
	clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, 256, vendor, &len);
	vendor[len] = '\0';

	clGetDeviceInfo(device_id, CL_DEVICE_TYPE, sizeof(type), &type, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(comp_units), &comp_units, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(freq), &freq, NULL);
	clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);

	printf("%s-%s -- %d units per %d MHz\n", vendor, name, comp_units, freq);
}

int main(int argc, char *argv[])
{

	const int LIST_SIZE = 128;
	const size_t globalWorkSize[] = { LIST_SIZE };
	cl_int res;
	cl_uint ret_num_devices, ret_num_platforms, alloc_num_devices = 0;
	cl_device_id *device_id = NULL, working_device;
	cl_platform_id *platform_ids = NULL;

	cl_context ctx = NULL;
	cl_command_queue que = NULL;
	cl_program program = NULL;
	cl_kernel kernel = NULL;
	cl_event ev;
	cl_mem mem_[] = { NULL, NULL, NULL };

	cl_uint i_ = 0, j_ = 0;
	int i;

    cl_uint device_index = CL_INT_MAX;
	char *arg;


	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') break;
		switch (argv[i][1]) {
		case 'd':
			i++;
			arg = _strdup(argv[i]);
			device_index = atoi(arg);
			free(arg);
			break;
		default:
			goto end_arg;
		}
		continue;
	end_arg: break;
	}

	res = clGetPlatformIDs(0, NULL, &ret_num_platforms);
	if (res != CL_SUCCESS) {
		printf("Error: clGetPlatformIDs returns %d\n", res);
		goto error;
	}

	platform_ids = (cl_platform_id *)malloc(sizeof(platform_ids) * ret_num_platforms);
	clGetPlatformIDs(ret_num_platforms, platform_ids, NULL);
	for (i_ = 0; i_ < ret_num_platforms; i_++) {
		printf("\nPlatform #%d::", i_);
		printPlatformInfo(platform_ids[i_]);

		res = clGetDeviceIDs(platform_ids[i_], CL_DEVICE_TYPE_ALL, 0, NULL, &ret_num_devices);
		if (res != CL_SUCCESS) {
			printf("Error: clGetDeviceIDs returns %d\n", res);
			goto error;
		}
		if (ret_num_devices< 1 || ret_num_devices > 50) {
			printf("Error: number of devices is out of range %u\n", ret_num_devices);
			goto error;
		}

		if (device_id)
			device_id = (cl_device_id *)realloc(device_id, sizeof(cl_device_id) * ret_num_devices + sizeof(cl_device_id) * alloc_num_devices);
		else
			device_id = (cl_device_id *)malloc(sizeof(cl_device_id) * ret_num_devices);

		res = clGetDeviceIDs(platform_ids[i_], CL_DEVICE_TYPE_ALL, ret_num_devices, &device_id[alloc_num_devices], NULL);

		for (j_ = alloc_num_devices; j_ < ret_num_devices + alloc_num_devices; j_++) {
			printf("\tDevice %d::", j_);
			printDeviceInfo(device_id[j_]);
		}
		alloc_num_devices += ret_num_devices;
	}
	
	ret_num_devices = alloc_num_devices;

	if (ret_num_devices > 1 && ret_num_devices < device_index) {
		printf("Select device number \n");
		goto error;
	}

	working_device = device_id[device_index];
	ctx = clCreateContext(0, 1, &working_device, NULL, NULL, &res);
	if (!ctx || res != CL_SUCCESS) {
		printf("Create context fails %d\n", res);
		goto error;
	}
		
	//que = clCreateCommandQueue(ctx, &device_id[0], 0, &res);
	que = clCreateCommandQueue(ctx, working_device, NULL, &res);
	if (!que || res != CL_SUCCESS) {
		printf("Failed to create queue\n");
		goto error;
	}

	program = CreateProgram(ctx, working_device, "HelloKernel.cl");
	if (!program)
	{
		printf("Failed to create program\n");
		goto error;
	}
	else
		printf("Program Ok\n");

	kernel = clCreateKernel(program, "hello_kernel", NULL);
	if (!kernel)
	{
		printf("Failed to create kernel\n");
		goto error;
	}
	else
		printf("Kernel Ok\n");


	start();

	res = clEnqueueNDRangeKernel(que, kernel, 1, NULL, globalWorkSize, NULL, 0, NULL, &ev);
	if (res != CL_SUCCESS) {
		printf("Failed to enqueue ND range\n");
		goto error;
	}

	clWaitForEvents(1, &ev);
	clReleaseEvent(ev);
	stop();
error:
	if (device_id) free(device_id);
	if (platform_ids) free(platform_ids);
	cleanupGPU(ctx, que, program, kernel, mem_);

	return 1;
}

void cleanupGPU(cl_context context, cl_command_queue commandQueue, cl_program program, cl_kernel kernel, cl_mem *memObjects) {
	if (commandQueue) {
		clFlush(commandQueue);
		clFinish(commandQueue);
		clReleaseCommandQueue(commandQueue);
	}
	if (kernel)
		clReleaseKernel(kernel);
	if (program)
		clReleaseProgram(program);
	if (context)
		clReleaseContext(context);

	for (int i = 0; i < 3; ++i)
	{
		if (memObjects[i])
		{
			clReleaseMemObject(memObjects[i]);
		}
	}
}

