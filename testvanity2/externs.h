#pragma once
#pragma comment(lib, "Ws2_32.lib")
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
#include "src/libsecp256k1-config.h"
//#include <secp256k1.c>
#include "msecp.h"

#include "b58.h"
#include "utils.h"

#define SHA256_DIGEST_LENGTH 32
#define RIPEMD160_DIGEST_LENGTH 20

#define THREADS_NUM 1
#define ADDRESS_NUM 1
#define ERROR_CREATE_THREAD     -11

#define SHA256(x,y,z) snprintf(z,y,"%s",x)
#define RIPEMD160(x,y,z) snprintf(z,y,"%s",x)
#define NELEM(array) (int)(sizeof(array)/sizeof(array[0]))

const char *adr_to_find[] = { "1LzhS3k3e9Ub8i2W1V8xQFdB8n2MYCHPCa", "17aPYR1m6pVAacXg1PTDDU7XafvK1dxvhi", "15c9mPGLku1HuW9LRtBf4jcHVpBUt8txKz",
"1Dn8NF8qDyyfHMktmuoQLGyjWmZXgvosXf", "1HAX2n9Uruu9YDt4cqRgYcvtGvZj1rbUyt", "1Kn5h2qpgw9mWE5jKpk8PP4qvvJ1QVy8su", "1AVJKwzs9AskraJLGHAZPiaZcrpDr1U6AB" };

const char fname_result[] = "found.txt";

