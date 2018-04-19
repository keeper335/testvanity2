#pragma once
#define USE_ECMULT_STATIC_PRECOMPUTATION 1
#ifndef NULL
//#define NULL ((void *)0)
#define NULL 0
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
	secp256k1_ge_storage *prec; /* prec[j][i] = 16^j * i * G + U_i */
	secp256k1_scalar blind;
	secp256k1_gej initial;
} secp256k1_ecmult_gen_context;


void my_memcpy(void* dst, void const* src, size_t size)
{
	int ret = 0;
	unsigned char *i_dst = (unsigned char *)dst;
	for (; ret < size; ret++)
		i_dst[ret] = *((unsigned char *)src + ret);
}

void my_memset(void* dst, int byte, size_t size)
{
	int ret = 0;
	unsigned char *i_dst = (unsigned char *)dst;
	for (; ret < size; ret++)
		i_dst[ret] = byte;
}

//#ifndef _INC_STRING
#if 0
void* __cdecl memcpy(
	void* _Dst,
	void const* _Src,
	size_t      _Size
);

void* __cdecl memset(
	void*  _Dst,
	int    _Val,
	size_t _Size
);
#else
#define memset my_memset
#define memcpy my_memcpy
#endif