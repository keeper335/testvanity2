// stdafx.h: включаемый файл для стандартных системных включаемых файлов
// или включаемых файлов для конкретного проекта, которые часто используются, но
// не часто изменяются
//

#pragma once
#pragma comment(lib, "OpenCL.lib")

#define _CRT_SECURE_NO_WARNINGS
#define __global

#include "targetver.h"

#include <stdio.h>
#include <stdint.h>
#include <tchar.h>
#include <CL/cl.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <windows.h>

void start();
void stop();

void cleanupGPU(cl_context context, cl_command_queue commandQueue, cl_program program, cl_kernel kernel, cl_mem *memObjects);

char * getFileSource(const char *fileName);
cl_program CreateProgram(cl_context context, cl_device_id device, const char *kernelFileName);

// TODO: Установите здесь ссылки на дополнительные заголовки, требующиеся для программы
