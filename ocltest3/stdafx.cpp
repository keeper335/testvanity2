// stdafx.cpp: исходный файл, содержащий только стандартные включаемые модули
// ocltest3.pch будет использоваться в качестве предкомпилированного заголовка
// stdafx.obj будет содержать предварительно откомпилированные сведения о типе

#include "stdafx.h"
#include "math.h"

// TODO: Установите ссылки на любые требующиеся дополнительные заголовки в файле STDAFX.H
// , а не в данном файле

double PCFreq = 0.0;
LONGLONG CounterStart = 0;

void start()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

void stop()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	printf("%f\n", double(li.QuadPart - CounterStart) / PCFreq);
}

char *getFileSource(const char *fileName) {
	char *fileBuffer;
	FILE *file = fopen(fileName, "r");
	if (!file)
	{
		fprintf(stderr, "Failed to open file %s\n", fileName);
		return NULL;
	}

	struct stat st;
	stat(fileName, &st);
	size_t fileSize = st.st_size;

	if (0 == fileSize)
	{
		fprintf(stderr, "%s source file is empty\n", fileName);
		return NULL;
	}

	fileBuffer = (char *)malloc(sizeof(char) * fileSize);
	size_t bytesRead = fread(fileBuffer, sizeof(char), fileSize, file);
	fclose(file);
	if (bytesRead != fileSize)
	{
		fprintf(stderr, "Failed to read complete source file %s, read %zd, needs %zd\n", fileName, bytesRead, fileSize);
		//free(fileBuffer);
		//return NULL;
		fileBuffer[bytesRead] = '\0';
		//printf("%s\n", fileBuffer);

	}
	return fileBuffer;
}

cl_program CreateProgram(cl_context context, cl_device_id device, const char *kernelFileName)
{
	cl_int errNum = 0;
	cl_program program = NULL;

	char *fileBuffers[1] = { NULL, };

	//fileBuffers[0] = getFileSource("secp256k1.cl");
	fileBuffers[0] = getFileSource(kernelFileName);
	if (fileBuffers[0]) {
		program = clCreateProgramWithSource(context, 1, (const char **)fileBuffers, NULL, NULL);
		free(fileBuffers[0]);
	}
//	free(fileBuffers[1]);

	if (!program)
	{
		fprintf(stderr, "Failed to create program from source\n");
		return NULL;
	}

	errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (CL_SUCCESS != errNum)
	{
		char buildLog[16384];
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
		fprintf(stderr, "Error when building program:\n%s\n", buildLog);
		clReleaseProgram(program);
		return NULL;
	}

	return program;
}
