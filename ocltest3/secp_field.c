#include "secp_static_ctx.c"

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
