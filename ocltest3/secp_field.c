#include "secp_static_ctx.c"

static void secp256k1_fe_mul_inner(uint64_t *r, const uint64_t *a, const uint64_t *b) {
	uint128_t c, d;
	uint64_t t3, t4, tx, u0;
	uint64_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	const uint64_t M = 0xFFFFFFFFFFFFFULL, R = 0x1000003D10ULL;

	d = (uint128_t)a0 * b[3] + (uint128_t)a1 * b[2]	+ (uint128_t)a2 * b[1]	+ (uint128_t)a3 * b[0];
	c = (uint128_t)a4 * b[4];
	d += (c & M) * R; c >>= 52;
	t3 = d & M; d >>= 52;
	d += (uint128_t)a0 * b[4]+ (uint128_t)a1 * b[3]+ (uint128_t)a2 * b[2]+ (uint128_t)a3 * b[1]+ (uint128_t)a4 * b[0];
	d += c * R;
	t4 = d & M; d >>= 52;
	tx = (t4 >> 48); t4 &= (M >> 4);
	c = (uint128_t)a0 * b[0];
	d += (uint128_t)a1 * b[4]+ (uint128_t)a2 * b[3]+ (uint128_t)a3 * b[2]+ (uint128_t)a4 * b[1];
	u0 = d & M; d >>= 52;
	u0 = (u0 << 4) | tx;
	c += (uint128_t)u0 * (R >> 4);
	r[0] = c & M; c >>= 52;
	c += (uint128_t)a0 * b[1]+ (uint128_t)a1 * b[0];
	d += (uint128_t)a2 * b[4]+ (uint128_t)a3 * b[3]	+ (uint128_t)a4 * b[2];
	c += (d & M) * R; d >>= 52;
	r[1] = c & M; c >>= 52;
	c += (uint128_t)a0 * b[2]+ (uint128_t)a1 * b[1]	+ (uint128_t)a2 * b[0];
	d += (uint128_t)a3 * b[4]+ (uint128_t)a4 * b[3];
	c += (d & M) * R; d >>= 52;
	r[2] = c & M; c >>= 52;
	c += d * R + t3;;
	r[3] = c & M; c >>= 52;
	c += t4;
	r[4] = c;
}

static void secp256k1_fe_sqr_inner(uint64_t *r, const uint64_t *a) {
	uint128_t c, d;
	uint64_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
	int64_t t3, t4, tx, u0;
	const uint64_t M = 0xFFFFFFFFFFFFFULL, R = 0x1000003D10ULL;

	d = (uint128_t)(a0 * 2) * a3 + (uint128_t)(a1 * 2) * a2;
	c = (uint128_t)a4 * a4;
	d += (c & M) * R; c >>= 52;
	t3 = d & M; d >>= 52;
	a4 *= 2;
	d += (uint128_t)a0 * a4	+ (uint128_t)(a1 * 2) * a3	+ (uint128_t)a2 * a2;
	d += c * R;
	t4 = d & M; d >>= 52;
	tx = (t4 >> 48); t4 &= (M >> 4);
	c = (uint128_t)a0 * a0;
	d += (uint128_t)a1 * a4	+ (uint128_t)(a2 * 2) * a3;
	u0 = d & M; d >>= 52;
	u0 = (u0 << 4) | tx;
	c += (uint128_t)u0 * (R >> 4);
	r[0] = c & M; c >>= 52;
	a0 *= 2;
	c += (uint128_t)a0 * a1;
	d += (uint128_t)a2 * a4	+ (uint128_t)a3 * a3;
	c += (d & M) * R; d >>= 52;
	r[1] = c & M; c >>= 52;
	c += (uint128_t)a0 * a2	+ (uint128_t)a1 * a1;
	d += (uint128_t)a3 * a4;
	c += (d & M) * R; d >>= 52;
	r[2] = c & M; c >>= 52;
	c += d * R + t3;;
	r[3] = c & M; c >>= 52;
	c += t4;
	r[4] = c;
}

static int secp256k1_fe_is_zero(const secp256k1_fe *a) {
	const uint64_t *t = a->n;
	return (t[0] | t[1] | t[2] | t[3] | t[4]) == 0;
}

static void secp256k1_fe_clear(secp256k1_fe *a) {
	int i;
	for (i = 0; i<5; i++) a->n[i] = 0;
}

static void secp256k1_fe_normalize(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t m;
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; m = t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; m &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; m &= t3;
	x = (t4 >> 48) | ((t4 == 0x0FFFFFFFFFFFFULL) & (m == 0xFFFFFFFFFFFFFULL) & (t0 >= 0xFFFFEFFFFFC2FULL));
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;
	t4 &= 0x0FFFFFFFFFFFFULL;
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
}

static void secp256k1_fe_normalize_weak(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];

	/* Reduce t4 at the start so there will be at most a single carry from the first pass */
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;

	/* The first pass ensures the magnitude is 1, ... */
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
}

static int secp256k1_fe_normalizes_to_zero(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t z0, z1;
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL; z0 = t0; z1 = t0 ^ 0x1000003D0ULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; z0 |= t1; z1 &= t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; z0 |= t2; z1 &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; z0 |= t3; z1 &= t3;
	z0 |= t4; z1 &= t4 ^ 0xF000000000000ULL;
	return (z0 == 0) | (z1 == 0xFFFFFFFFFFFFFULL);
}

static void secp256k1_fe_negate(secp256k1_fe *r, const secp256k1_fe *a, int m) {
	r->n[0] = 0xFFFFEFFFFFC2FULL * 2 * (m + 1) - a->n[0];
	r->n[1] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[1];
	r->n[2] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[2];
	r->n[3] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[3];
	r->n[4] = 0x0FFFFFFFFFFFFULL * 2 * (m + 1) - a->n[4];
}

static void secp256k1_fe_set_int(secp256k1_fe *r, int a) {
	r->n[0] = a;
	r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
}

static int secp256k1_fe_set_b32(secp256k1_fe *r, const unsigned char *a) {
	int i;
	r->n[0] = r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
	for (i = 0; i<32; i++) {
		int j;
		for (j = 0; j<2; j++) {
			int limb = (8 * i + 4 * j) / 52;
			int shift = (8 * i + 4 * j) % 52;
			r->n[limb] |= (uint64_t)((a[31 - i] >> (4 * j)) & 0xF) << shift;
		}
	}
	if (r->n[4] == 0x0FFFFFFFFFFFFULL && (r->n[3] & r->n[2] & r->n[1]) == 0xFFFFFFFFFFFFFULL && r->n[0] >= 0xFFFFEFFFFFC2FULL) {
		return 0;
	}
	return 1;
}


static void secp256k1_fe_sqr(secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe_sqr_inner(r->n, a->n);
}

static void secp256k1_fe_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b) {
	secp256k1_fe_mul_inner(r->n, a->n, b->n);
}

static void secp256k1_fe_cmov(secp256k1_fe *r, const secp256k1_fe *a, int flag) {
	uint64_t mask0, mask1;
	mask0 = flag + ~((uint64_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
	r->n[4] = (r->n[4] & mask0) | (a->n[4] & mask1);
}

static void secp256k1_fe_storage_cmov(secp256k1_fe_storage *r, const secp256k1_fe_storage *a, int flag) {
	uint64_t mask0, mask1;
	mask0 = flag + ~((uint64_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
}

static void secp256k1_fe_to_storage(secp256k1_fe_storage *r, const secp256k1_fe *a) {
	r->n[0] = a->n[0] | a->n[1] << 52;
	r->n[1] = a->n[1] >> 12 | a->n[2] << 40;
	r->n[2] = a->n[2] >> 24 | a->n[3] << 28;
	r->n[3] = a->n[3] >> 36 | a->n[4] << 16;
}

static void secp256k1_fe_from_storage(secp256k1_fe *r, const secp256k1_fe_storage *a) {
	r->n[0] = a->n[0] & 0xFFFFFFFFFFFFFULL;
	r->n[1] = a->n[0] >> 52 | ((a->n[1] << 12) & 0xFFFFFFFFFFFFFULL);
	r->n[2] = a->n[1] >> 40 | ((a->n[2] << 24) & 0xFFFFFFFFFFFFFULL);
	r->n[3] = a->n[2] >> 28 | ((a->n[3] << 36) & 0xFFFFFFFFFFFFFULL);
	r->n[4] = a->n[3] >> 16;
}

static void secp256k1_fe_add(secp256k1_fe *r, const secp256k1_fe *a) {
	r->n[0] += a->n[0];
	r->n[1] += a->n[1];
	r->n[2] += a->n[2];
	r->n[3] += a->n[3];
	r->n[4] += a->n[4];
}

static void secp256k1_fe_mul_int(secp256k1_fe *r, int a) {
	r->n[0] *= a;
	r->n[1] *= a;
	r->n[2] *= a;
	r->n[3] *= a;
	r->n[4] *= a;
}
