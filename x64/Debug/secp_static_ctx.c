#pragma once
#define USE_ECMULT_STATIC_PRECOMPUTATION 1
#ifndef NULL
//#define NULL ((void *)0)
#define NULL 0
#endif

#ifndef __constant
#define __constant
#endif
#ifndef __global
#define __global
#endif

#ifndef _STDINT
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long      uint64_t;
typedef long  int64_t;
typedef int   int32_t;
typedef short int16_t;
typedef char  int8_t;
#endif

typedef struct { uint32_t n[10]; } secp256k1_fe;
typedef struct { uint32_t n[8]; } secp256k1_fe_storage;
typedef struct { uint32_t d[8]; } secp256k1_scalar;
//typedef uint64_t uint128_t;
//typedef int64_t int128_t;

typedef struct {
	secp256k1_fe x;
	secp256k1_fe y;
	int infinity; /* whether this represents the point at infinity */
} secp256k1_ge;

typedef struct {
	secp256k1_fe x; /* actual X: x/z^2 */
	secp256k1_fe y; /* actual Y: y/z^3 */
	secp256k1_fe z;
	int infinity; /* whether this represents the point at infinity */
} secp256k1_gej;

typedef struct {
	secp256k1_fe_storage x;
	secp256k1_fe_storage y;
} secp256k1_ge_storage;

typedef struct {
	secp256k1_ge_storage(*prec)[64][16]; /* prec[j][i] = 16^j * i * G + U_i */
	secp256k1_scalar blind;
	secp256k1_gej initial;
} secp256k1_ecmult_gen_context;

#ifndef _VCRUNTIME_H
void memcpy(void* dst, void const* src, size_t size)
{
	size_t ret = 0;
	char *i_dst = (char *)dst;
	for (; ret <= size; ret++)
		i_dst[ret] = ((char *)src)[ret];
}

void memset(void* dst, int byte, size_t size)
{
	size_t ret = 0;
	int *i_dst = (int *)dst;
	for (; ret <= size; ret++)
		i_dst[ret] = byte;
}
#endif