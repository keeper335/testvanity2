#pragma once
#define USE_ECMULT_STATIC_PRECOMPUTATION 1
#ifndef NULL
//#define NULL ((void *)0)
#define NULL 0
#endif

typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long      uint64_t;
typedef long  int64_t;
typedef int   int32_t;
typedef short int16_t;
typedef char  int8_t;

typedef struct { uint32_t n[10]; } secp256k1_fe;
typedef struct { uint32_t n[8] } secp256k1_fe_storage;
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

#define SECP256K1_FE_CONST_INNER(d7, d6, d5, d4, d3, d2, d1, d0) { \
    (d0) & 0x3FFFFFFUL, \
    (((uint32_t)d0) >> 26) | (((uint32_t)(d1) & 0xFFFFFUL) << 6), \
    (((uint32_t)d1) >> 20) | (((uint32_t)(d2) & 0x3FFFUL) << 12), \
    (((uint32_t)d2) >> 14) | (((uint32_t)(d3) & 0xFFUL) << 18), \
    (((uint32_t)d3) >> 8) | (((uint32_t)(d4) & 0x3UL) << 24), \
    (((uint32_t)d4) >> 2) & 0x3FFFFFFUL, \
    (((uint32_t)d4) >> 28) | (((uint32_t)(d5) & 0x3FFFFFUL) << 4), \
    (((uint32_t)d5) >> 22) | (((uint32_t)(d6) & 0xFFFFUL) << 10), \
    (((uint32_t)d6) >> 16) | (((uint32_t)(d7) & 0x3FFUL) << 16), \
    (((uint32_t)d7) >> 10) \
}
#define SECP256K1_FE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {SECP256K1_FE_CONST_INNER((d7), (d6), (d5), (d4), (d3), (d2), (d1), (d0))}
#define SECP256K1_FE_STORAGE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{ (d0), (d1), (d2), (d3), (d4), (d5), (d6), (d7) }}
#define SECP256K1_FE_STORAGE_CONST_GET(d) d.n[7], d.n[6], d.n[5], d.n[4],d.n[3], d.n[2], d.n[1], d.n[0]
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
