#pragma once
#pragma comment(lib, "WS2_32.lib")
#pragma comment(lib, "secp256k1.lib")
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include <memory.h>
#include <string.h>


/* include secp256k1 */
//#include "src/libsecp256k1-config.h"
#include <secp256k1.c>
//#include "secp256_ext.h"

#include "b58.h"
#include "utils.h"

#define SHA256_DIGEST_LENGTH 32
#define RIPEMD160_DIGEST_LENGTH 20
#define STEP 3072
#define THREADS_NUM 1
#define ADDRESS_NUM 1
#define ERROR_CREATE_THREAD     -11

#define SHA256(x,y,z) snprintf(z,y,"%s",x)
#define RIPEMD160(x,y,z) snprintf(z,y,"%s",x)
#define NELEM(array) (int)(sizeof(array)/sizeof(array[0]))

