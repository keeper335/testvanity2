#include "secp_static_ctx.c"

/* Limbs of the secp256k1 order. */
#define SECP256K1_N_0 ((uint64_t)0xBFD25E8CD0364141ULL)
#define SECP256K1_N_1 ((uint64_t)0xBAAEDCE6AF48A03BULL)
#define SECP256K1_N_2 ((uint64_t)0xFFFFFFFFFFFFFFFEULL)
#define SECP256K1_N_3 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)

/* Limbs of 2^256 minus the secp256k1 order. */
#define SECP256K1_N_C_0 (~SECP256K1_N_0 + 1)
#define SECP256K1_N_C_1 (~SECP256K1_N_1)
#define SECP256K1_N_C_2 (1)

/* Limbs of half the secp256k1 order. */
#define SECP256K1_N_H_0 ((uint64_t)0xDFE92F46681B20A0ULL)
#define SECP256K1_N_H_1 ((uint64_t)0x5D576E7357A4501DULL)
#define SECP256K1_N_H_2 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)
#define SECP256K1_N_H_3 ((uint64_t)0x7FFFFFFFFFFFFFFFULL)

static void secp256k1_scalar_clear(secp256k1_scalar *r) {
	r->d[0] = 0;
	r->d[1] = 0;
	r->d[2] = 0;
	r->d[3] = 0;
}

static void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v) {
	r->d[0] = v;
	r->d[1] = 0;
	r->d[2] = 0;
	r->d[3] = 0;
}

static int secp256k1_scalar_check_overflow(const secp256k1_scalar *a) {
	int yes = 0;
	int no = 0;
	no |= (a->d[3] < SECP256K1_N_3); /* No need for a > check. */
	no |= (a->d[2] < SECP256K1_N_2);
	yes |= (a->d[2] > SECP256K1_N_2) & ~no;
	no |= (a->d[1] < SECP256K1_N_1);
	yes |= (a->d[1] > SECP256K1_N_1) & ~no;
	yes |= (a->d[0] >= SECP256K1_N_0) & ~no;
	return yes;
}

static int secp256k1_scalar_is_zero(const secp256k1_scalar *a) {
	return (a->d[0] | a->d[1] | a->d[2] | a->d[3]) == 0;
}

static unsigned int secp256k1_scalar_get_bits(const secp256k1_scalar *a, unsigned int offset, unsigned int count) {
	return (a->d[offset >> 6] >> (offset & 0x3F)) & ((((uint64_t)1) << count) - 1);
}

static int secp256k1_scalar_reduce(secp256k1_scalar *r, unsigned int overflow) {
	uint128_t t;
	t = (uint128_t)r->d[0] + overflow * SECP256K1_N_C_0;
	r->d[0] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint128_t)r->d[1] + overflow * SECP256K1_N_C_1;
	r->d[1] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint128_t)r->d[2] + overflow * SECP256K1_N_C_2;
	r->d[2] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint64_t)r->d[3];
	r->d[3] = t & 0xFFFFFFFFFFFFFFFFULL;
	return overflow;
}

static int secp256k1_scalar_add(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b) {
	int overflow;
	uint128_t t = (uint128_t)a->d[0] + b->d[0];
	r->d[0] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint128_t)a->d[1] + b->d[1];
	r->d[1] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint128_t)a->d[2] + b->d[2];
	r->d[2] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	t += (uint128_t)a->d[3] + b->d[3];
	r->d[3] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
	overflow = t + secp256k1_scalar_check_overflow(r);
	secp256k1_scalar_reduce(r, overflow);
	return overflow;
}

static void secp256k1_scalar_set_b32(secp256k1_scalar *r, const unsigned char *b32, int *overflow) {
	int over;
	r->d[0] = (uint64_t)b32[31] | (uint64_t)b32[30] << 8 | (uint64_t)b32[29] << 16 | (uint64_t)b32[28] << 24 | (uint64_t)b32[27] << 32 | (uint64_t)b32[26] << 40 | (uint64_t)b32[25] << 48 | (uint64_t)b32[24] << 56;
	r->d[1] = (uint64_t)b32[23] | (uint64_t)b32[22] << 8 | (uint64_t)b32[21] << 16 | (uint64_t)b32[20] << 24 | (uint64_t)b32[19] << 32 | (uint64_t)b32[18] << 40 | (uint64_t)b32[17] << 48 | (uint64_t)b32[16] << 56;
	r->d[2] = (uint64_t)b32[15] | (uint64_t)b32[14] << 8 | (uint64_t)b32[13] << 16 | (uint64_t)b32[12] << 24 | (uint64_t)b32[11] << 32 | (uint64_t)b32[10] << 40 | (uint64_t)b32[9] << 48 | (uint64_t)b32[8] << 56;
	r->d[3] = (uint64_t)b32[7] | (uint64_t)b32[6] << 8 | (uint64_t)b32[5] << 16 | (uint64_t)b32[4] << 24 | (uint64_t)b32[3] << 32 | (uint64_t)b32[2] << 40 | (uint64_t)b32[1] << 48 | (uint64_t)b32[0] << 56;
	over = secp256k1_scalar_reduce(r, secp256k1_scalar_check_overflow(r));
	if (overflow) {
		*overflow = over;
	}
}

static void secp256k1_scalar_get_b32(unsigned char *bin, const secp256k1_scalar* a) {
	bin[0] = a->d[3] >> 56; bin[1] = a->d[3] >> 48; bin[2] = a->d[3] >> 40; bin[3] = a->d[3] >> 32; bin[4] = a->d[3] >> 24; bin[5] = a->d[3] >> 16; bin[6] = a->d[3] >> 8; bin[7] = a->d[3];
	bin[8] = a->d[2] >> 56; bin[9] = a->d[2] >> 48; bin[10] = a->d[2] >> 40; bin[11] = a->d[2] >> 32; bin[12] = a->d[2] >> 24; bin[13] = a->d[2] >> 16; bin[14] = a->d[2] >> 8; bin[15] = a->d[2];
	bin[16] = a->d[1] >> 56; bin[17] = a->d[1] >> 48; bin[18] = a->d[1] >> 40; bin[19] = a->d[1] >> 32; bin[20] = a->d[1] >> 24; bin[21] = a->d[1] >> 16; bin[22] = a->d[1] >> 8; bin[23] = a->d[1];
	bin[24] = a->d[0] >> 56; bin[25] = a->d[0] >> 48; bin[26] = a->d[0] >> 40; bin[27] = a->d[0] >> 32; bin[28] = a->d[0] >> 24; bin[29] = a->d[0] >> 16; bin[30] = a->d[0] >> 8; bin[31] = a->d[0];
}

static void secp256k1_scalar_negate(secp256k1_scalar *r, const secp256k1_scalar *a) {
	uint64_t nonzero = 0xFFFFFFFFFFFFFFFFULL * (secp256k1_scalar_is_zero(a) == 0);
	uint128_t t = (uint128_t)(~a->d[0]) + SECP256K1_N_0 + 1;
	r->d[0] = t & nonzero; t >>= 64;
	t += (uint128_t)(~a->d[1]) + SECP256K1_N_1;
	r->d[1] = t & nonzero; t >>= 64;
	t += (uint128_t)(~a->d[2]) + SECP256K1_N_2;
	r->d[2] = t & nonzero; t >>= 64;
	t += (uint128_t)(~a->d[3]) + SECP256K1_N_3;
	r->d[3] = t & nonzero;
}
