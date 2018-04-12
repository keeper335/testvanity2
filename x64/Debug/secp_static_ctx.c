#pragma once
#define USE_ECMULT_STATIC_PRECOMPUTATION 1
#ifndef NULL
//#define NULL ((void *)0)
#define NULL 0
#endif

typedef unsigned long long uint64_t;
typedef unsigned long uint32_t;
typedef unsigned int uint16_t;
//typedef unsigned int size_t;
typedef int ssize_t;
typedef long long int64_t;
typedef long int32_t;
typedef int int16_t;
typedef unsigned char uint8_t;
typedef struct { uint64_t n[5]; } secp256k1_fe;
typedef struct { uint64_t n[4]; } secp256k1_fe_storage;
typedef struct { uint64_t d[4]; } secp256k1_scalar;
typedef uint64_t uint128_t;
typedef int64_t int128_t;

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

#define SECP256K1_FE_CONST_INNER(d7, d6, d5, d4, d3, d2, d1, d0) { \
    (d0) | (((uint64_t)(d1) & 0xFFFFFUL) << 32), \
    ((uint64_t)(d1) >> 20) | (((uint64_t)(d2)) << 12) | (((uint64_t)(d3) & 0xFFUL) << 44), \
    ((uint64_t)(d3) >> 8) | (((uint64_t)(d4) & 0xFFFFFFFUL) << 24), \
    ((uint64_t)(d4) >> 28) | (((uint64_t)(d5)) << 4) | (((uint64_t)(d6) & 0xFFFFUL) << 36), \
    ((uint64_t)(d6) >> 16) | (((uint64_t)(d7)) << 16) \
}
#define SECP256K1_FE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {SECP256K1_FE_CONST_INNER((d7), (d6), (d5), (d4), (d3), (d2), (d1), (d0))}
#define SECP256K1_GE_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_CONST((i),(j),(k),(l),(m),(n),(o),(p)), 0}

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
