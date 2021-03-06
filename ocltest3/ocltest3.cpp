// ocltest3.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "../x64/Debug/secp.c"
#include "../testvanity2/b58.h"
#define __global

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

	const int LIST_SIZE = 1;
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
	int i, j;

    cl_uint device_index = CL_INT_MAX;
	char *arg;
	bool print_info = true;


	//check ctx
	secp256k1_ecmult_gen_context sec_ctx;
	sec_ctx.prec = (secp256k1_ge_storage(*)[64][16])malloc(sizeof(*sec_ctx.prec));

	secp256k1_ecmult_gen_context_build(&sec_ctx);
	printf("secp256k1_ecmult_gen_context_init\n");
	//secp256k1_ecmult_gen_context_build(&sec_ctx);
	secp256k1_scalar scalar_key = { 0, }, scalar_key_one = { {1} };
	uint8_t *p = (uint8_t *)scalar_key.d;
	secp256k1_gej gej, gej2, gej_one;
	secp256k1_ge ge;
	uint8_t sha_block[64] = { 0, }, sha_block2[20 * STEP] = { 0, };
	unsigned char sha_text[128] = { 0, };
	char p_text[128] = { 0, };

	p[0] = 0x34; p[1] = 0x12;
	sprintf(p_text, "%08x%08x%08x%08x%08x%08x%08x%08x", scalar_key.d[7], scalar_key.d[6],scalar_key.d[5], scalar_key.d[4],\
			scalar_key.d[3], scalar_key.d[2], scalar_key.d[1], scalar_key.d[0]);
	printf("Private key %s\n", p_text);
	secp256k1_ecmult_gen(&sec_ctx, &gej, &scalar_key);
	secp256k1_ge_set_gej_var(&ge, &gej);
	sha_block[0] = (secp256k1_fe_is_odd(&ge.y) ? 0x03 : 0x02);
	secp256k1_fe_get_b32(sha_block + 1, &ge.x);
	for (i = 0; i < 33; i++)
		sprintf((char *)&sha_text[i * 2], "%02x", sha_block[i]);
	printf("Public key %s\n", sha_text);

	unsigned char rmd_block[21], chksum[32], rmd_text[64] = { 0, }, buffer[26] = { 0, };
	HASH160(rmd_block, sha_block);
	buffer[0] = 0;
	memcpy(buffer + 1, rmd_block, 20);
	HASH256(chksum, buffer);
	memcpy(buffer + 21, chksum, 4);
	b58enc(rmd_text, buffer, 25);
	printf("Public address %s\n", rmd_text);

	//return 1;

	//0x00-001234
	//0337a4aef1f8423ca076e4b7d99a8cabff40ddb8231f2a9f01081f15d7fa65c1ba
	//9c44f51e29bc5f8e9723c166cfb18b8ce2065c9c4c4a4bfbca141286e8451d20 - sha
	//ec8f7d48e2add2a176b420fc8c1eed6cf714b654 - hash160
	//1NZpN6D8vby5WpUaiSbGd53imC6Qfk2gyM

	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') break;
		switch (argv[i][1]) {
		case 'd':
			i++;
			arg = _strdup(argv[i]);
			device_index = atoi(arg);
			print_info = false;
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
		if (print_info) {
			printf("\nPlatform #%d::", i_);
			printPlatformInfo(platform_ids[i_]);
		}

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

		if (print_info) {
			for (j_ = alloc_num_devices; j_ < ret_num_devices + alloc_num_devices; j_++) {
				printf("\tDevice %d::", j_);
				printDeviceInfo(device_id[j_]);
			}
		}
		alloc_num_devices += ret_num_devices;
	}
	
	ret_num_devices = alloc_num_devices;

	if (ret_num_devices > 1 && ret_num_devices <= device_index) {printf("Select device number \n");	goto error;	}

	working_device = device_id[device_index];
	ctx = clCreateContext(0, 1, &working_device, NULL, NULL, &res);
	if (!ctx || res != CL_SUCCESS) {printf("Create context fails %d\n", res);goto error;}
	
	que = clCreateCommandQueue(ctx, working_device, NULL, &res);
	if (!que || res != CL_SUCCESS) {printf("Failed to create queue\n");	goto error;}

	program = CreateProgram(ctx, working_device, "HelloKernel.cl");
	if (!program) {	printf("Failed to create program\n");	goto error;	}

	kernel = clCreateKernel(program, "hello_kernel", NULL);
	if (!kernel){printf("Failed to create kernel\n");goto error;}

	start();

	secp256k1_ecmult_gen(&sec_ctx, &gej2, &scalar_key);
	secp256k1_ecmult_gen(&sec_ctx, &gej_one, &scalar_key_one);

	mem_[0] = clCreateBuffer(ctx, CL_MEM_READ_WRITE, sizeof(sha_block2), NULL, &res);
	mem_[1] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(gej2), &gej2, &res);
	mem_[2] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(gej_one), &gej_one, &res);

	res = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&mem_[0]);
	res = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&mem_[1]);
	res = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&mem_[2]);

	res = clEnqueueNDRangeKernel(que, kernel, 1, NULL, globalWorkSize, NULL, 0, NULL, &ev);
	if (res != CL_SUCCESS) {printf("Failed to enqueue ND range\n");	goto error;}

	//clWaitForEvents(1, &ev);
	//clReleaseEvent(ev);

	res = clEnqueueReadBuffer(que, mem_[0], CL_TRUE, 0, sizeof(sha_block2), sha_block2, 0, NULL, NULL);
	if (res != CL_SUCCESS) {printf("Failed to read mem buffer\n");	goto error;}

	stop();

	for (j = STEP-10; j < STEP; j++) {
		for (i = 0; i < 20; i++) {
			sprintf((char *)&sha_text[i * 2], "%02x", sha_block2[j*20 + i]);
		}
		printf("HASH160 %d -> %s\n", j, sha_text);
	}
		
	


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

