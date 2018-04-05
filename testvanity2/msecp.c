#include "msecp.h"

static void default_illegal_callback_fn(const char* str, void* data) {
	(void)data;
	fprintf(stderr, "[libsecp256k1] illegal argument: %s\n", str);
	abort();
}

static const secp256k1_callback default_illegal_callback = {
	default_illegal_callback_fn,
	NULL
};

static void default_error_callback_fn(const char* str, void* data) {
	(void)data;
	fprintf(stderr, "[libsecp256k1] internal consistency check failed: %s\n", str);
	abort();
}

static const secp256k1_callback default_error_callback = {
	default_error_callback_fn,
	NULL
};

void secp256k1_ecmult_context_init(secp256k1_ecmult_context *ctx) {
	ctx->pre_g = NULL;
}

static void secp256k1_ecmult_odd_multiples_table(int n, secp256k1_gej *prej, secp256k1_fe *zr, const secp256k1_gej *a) {
	secp256k1_gej d;
	secp256k1_ge a_ge, d_ge;
	int i;
	VERIFY_CHECK(!a->infinity);
	secp256k1_gej_double_var(&d, a, NULL);
	d_ge.x = d.x;
	d_ge.y = d.y;
	d_ge.infinity = 0;
	secp256k1_ge_set_gej_zinv(&a_ge, a, &d.z);
	prej[0].x = a_ge.x;
	prej[0].y = a_ge.y;
	prej[0].z = a->z;
	prej[0].infinity = 0;

	zr[0] = d.z;
	for (i = 1; i < n; i++)
		secp256k1_gej_add_ge_var(&prej[i], &prej[i - 1], &d_ge, &zr[i]);
	secp256k1_fe_mul(&prej[n - 1].z, &prej[n - 1].z, &d.z);
}

static void secp256k1_ecmult_odd_multiples_table_globalz_windowa(secp256k1_ge *pre, secp256k1_fe *globalz, const secp256k1_gej *a) {
	secp256k1_gej prej[ECMULT_TABLE_SIZE(WINDOW_A)];
	secp256k1_fe zr[ECMULT_TABLE_SIZE(WINDOW_A)];
	secp256k1_ecmult_odd_multiples_table(ECMULT_TABLE_SIZE(WINDOW_A), prej, zr, a);
	secp256k1_ge_globalz_set_table_gej(ECMULT_TABLE_SIZE(WINDOW_A), pre, globalz, prej, zr);
}

static void secp256k1_ecmult_odd_multiples_table_storage_var(int n, secp256k1_ge_storage *pre, const secp256k1_gej *a, const secp256k1_callback *cb) {
	secp256k1_gej *prej = (secp256k1_gej*)checked_malloc(cb, sizeof(secp256k1_gej) * n);
	secp256k1_ge *prea = (secp256k1_ge*)checked_malloc(cb, sizeof(secp256k1_ge) * n);
	secp256k1_fe *zr = (secp256k1_fe*)checked_malloc(cb, sizeof(secp256k1_fe) * n);
	int i;
	secp256k1_ecmult_odd_multiples_table(n, prej, zr, a);
	secp256k1_ge_set_table_gej_var(n, prea, prej, zr);
	for (i = 0; i < n; i++)
		secp256k1_ge_to_storage(&pre[i], &prea[i]);
	free(prea);
	free(prej);
	free(zr);
}


void secp256k1_ecmult_context_build(secp256k1_ecmult_context *ctx, const secp256k1_callback *cb) {
	secp256k1_gej gj;
	if (ctx->pre_g != NULL)
		return;
	secp256k1_gej_set_ge(&gj, &secp256k1_ge_const_g);
	ctx->pre_g = (secp256k1_ge_storage(*)[])checked_malloc(cb, sizeof((*ctx->pre_g)[0]) * ECMULT_TABLE_SIZE(WINDOW_G));
	secp256k1_ecmult_odd_multiples_table_storage_var(ECMULT_TABLE_SIZE(WINDOW_G), *ctx->pre_g, &gj, cb);
}

static void secp256k1_ecmult_context_clone(secp256k1_ecmult_context *dst,
	const secp256k1_ecmult_context *src, const secp256k1_callback *cb) {
	if (src->pre_g == NULL) {
		dst->pre_g = NULL;
	}
	else {
		size_t size = sizeof((*dst->pre_g)[0]) * ECMULT_TABLE_SIZE(WINDOW_G);
		dst->pre_g = (secp256k1_ge_storage(*)[])checked_malloc(cb, size);
		memcpy(dst->pre_g, src->pre_g, size);
	}
}

static int secp256k1_ecmult_context_is_built(const secp256k1_ecmult_context *ctx) {
	return ctx->pre_g != NULL;
}

static void secp256k1_ecmult_context_clear(secp256k1_ecmult_context *ctx) {
	free(ctx->pre_g);
	secp256k1_ecmult_context_init(ctx);
}

int secp256k1_ecmult_wnaf(int *wnaf, int len, const secp256k1_scalar *a, int w) {
	secp256k1_scalar s = *a;
	int last_set_bit = -1;
	int bit = 0;
	int sign = 1;
	int carry = 0;
	VERIFY_CHECK(wnaf != NULL);
	VERIFY_CHECK(0 <= len && len <= 256);
	VERIFY_CHECK(a != NULL);
	VERIFY_CHECK(2 <= w && w <= 31);
	memset(wnaf, 0, len * sizeof(wnaf[0]));
	if (secp256k1_scalar_get_bits(&s, 255, 1)) {
		secp256k1_scalar_negate(&s, &s);
		sign = -1;
	}
	while (bit < len) {
		int now;
		int word;
		if (secp256k1_scalar_get_bits(&s, bit, 1) == (unsigned int)carry) {
			bit++;
			continue;
		}
		now = w;
		if (now > len - bit) {
			now = len - bit;
		}
		word = secp256k1_scalar_get_bits_var(&s, bit, now) + carry;
		carry = (word >> (w - 1)) & 1;
		word -= carry << w;
		wnaf[bit] = sign * word;
		last_set_bit = bit;
		bit += now;
	}
	return last_set_bit + 1;
}

void secp256k1_ecmult(const secp256k1_ecmult_context *ctx, secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_scalar *na, const secp256k1_scalar *ng) {
	secp256k1_ge pre_a[ECMULT_TABLE_SIZE(WINDOW_A)];
	secp256k1_ge tmpa;
	secp256k1_fe Z;

	int wnaf_na[256];
	int bits_na;
	int wnaf_ng[256];
	int bits_ng;

	int i;
	int bits;
	/* build wnaf representation for na. */
	bits_na = secp256k1_ecmult_wnaf(wnaf_na, 256, na, WINDOW_A);
	bits = bits_na;
	secp256k1_ecmult_odd_multiples_table_globalz_windowa(pre_a, &Z, a);
	bits_ng = secp256k1_ecmult_wnaf(wnaf_ng, 256, ng, WINDOW_G);
	if (bits_ng > bits) {
		bits = bits_ng;
	}
	secp256k1_gej_set_infinity(r);

	for (i = bits - 1; i >= 0; i--) {
		int n;
		secp256k1_gej_double_var(r, r, NULL);
		if (i < bits_na && (n = wnaf_na[i])) {
			ECMULT_TABLE_GET_GE(&tmpa, pre_a, n, WINDOW_A);
			secp256k1_gej_add_ge_var(r, r, &tmpa, NULL);
		}
		if (i < bits_ng && (n = wnaf_ng[i])) {
			ECMULT_TABLE_GET_GE_STORAGE(&tmpa, *ctx->pre_g, n, WINDOW_G);
			secp256k1_gej_add_zinv_var(r, r, &tmpa, &Z);
		}
	}

	if (!r->infinity) {
		secp256k1_fe_mul(&r->z, &r->z, &Z);
	}
}

/* ecmult gen h */
void secp256k1_ecmult_gen_context_init(secp256k1_ecmult_gen_context *ctx) {
	ctx->prec = NULL;
}

void secp256k1_ecmult_gen_context_build(secp256k1_ecmult_gen_context *ctx, const secp256k1_callback* cb) {
	secp256k1_ge prec[1024];
	secp256k1_gej gj;
	secp256k1_gej nums_gej;
	int i, j;

	if (ctx->prec != NULL) {
		return;
	}

	ctx->prec = (secp256k1_ge_storage(*)[64][16])malloc(sizeof(*ctx->prec));

	/* get the generator */
	secp256k1_gej_set_ge(&gj, &secp256k1_ge_const_g);

	/* Construct a group element with no known corresponding scalar (nothing up my sleeve). */
	{
		static const unsigned char nums_b32[33] = "The scalar for this x is unknown";
		secp256k1_fe nums_x;
		secp256k1_ge nums_ge;
		VERIFY_CHECK(secp256k1_fe_set_b32(&nums_x, nums_b32));
		VERIFY_CHECK(secp256k1_ge_set_xo_var(&nums_ge, &nums_x, 0));
		secp256k1_gej_set_ge(&nums_gej, &nums_ge);
		/* Add G to make the bits in x uniformly distributed. */
		secp256k1_gej_add_ge_var(&nums_gej, &nums_gej, &secp256k1_ge_const_g, NULL);
	}

	/* compute prec. */
	{
		secp256k1_gej precj[1024]; /* Jacobian versions of prec. */
		secp256k1_gej gbase;
		secp256k1_gej numsbase;
		gbase = gj; /* 16^j * G */
		numsbase = nums_gej; /* 2^j * nums. */
		for (j = 0; j < 64; j++) {
			/* Set precj[j*16 .. j*16+15] to (numsbase, numsbase + gbase, ..., numsbase + 15*gbase). */
			precj[j * 16] = numsbase;
			for (i = 1; i < 16; i++) {
				secp256k1_gej_add_var(&precj[j * 16 + i], &precj[j * 16 + i - 1], &gbase, NULL);
			}
			/* Multiply gbase by 16. */
			for (i = 0; i < 4; i++) {
				secp256k1_gej_double_var(&gbase, &gbase, NULL);
			}
			/* Multiply numbase by 2. */
			secp256k1_gej_double_var(&numsbase, &numsbase, NULL);
			if (j == 62) {
				/* In the last iteration, numsbase is (1 - 2^j) * nums instead. */
				secp256k1_gej_neg(&numsbase, &numsbase);
				secp256k1_gej_add_var(&numsbase, &numsbase, &nums_gej, NULL);
			}
		}
		secp256k1_ge_set_all_gej_var(1024, prec, precj, cb);
	}
	for (j = 0; j < 64; j++) {
		for (i = 0; i < 16; i++) {
			secp256k1_ge_to_storage(&(*ctx->prec)[j][i], &prec[j * 16 + i]);
		}
	}

	secp256k1_ecmult_gen_blind(ctx, NULL);
}
void secp256k1_ecmult_gen_context_clear(secp256k1_ecmult_gen_context *ctx) {
	free(ctx->prec);
	secp256k1_scalar_clear(&ctx->blind);
	secp256k1_gej_clear(&ctx->initial);
	ctx->prec = NULL;
}

void secp256k1_ecmult_gen(const secp256k1_ecmult_gen_context *ctx, secp256k1_gej *r, const secp256k1_scalar *gn) {
	secp256k1_ge add;
	secp256k1_ge_storage adds;
	secp256k1_scalar gnb;
	int bits;
	int i, j;
	memset(&adds, 0, sizeof(adds));
	*r = ctx->initial;
	secp256k1_scalar_add(&gnb, gn, &ctx->blind);
	add.infinity = 0;
	for (j = 0; j < 64; j++) {
		bits = secp256k1_scalar_get_bits(&gnb, j * 4, 4);
		for (i = 0; i < 16; i++) {
			secp256k1_ge_storage_cmov(&adds, &(*ctx->prec)[j][i], i == bits);
		}
		secp256k1_ge_from_storage(&add, &adds);
		secp256k1_gej_add_ge(r, r, &add);
	}
	bits = 0;
	secp256k1_ge_clear(&add);
	secp256k1_scalar_clear(&gnb);
}

void secp256k1_ecmult_gen_blind(secp256k1_ecmult_gen_context *ctx, const unsigned char *seed32) {
	secp256k1_scalar b;
	secp256k1_gej gb;
	secp256k1_fe s;
	unsigned char nonce32[32];
	secp256k1_rfc6979_hmac_sha256_t rng;
	int retry;
	unsigned char keydata[64] = { 0 };
	if (seed32 == NULL) {
		/* When seed is NULL, reset the initial point and blinding value. */
		secp256k1_gej_set_ge(&ctx->initial, &secp256k1_ge_const_g);
		secp256k1_gej_neg(&ctx->initial, &ctx->initial);
		secp256k1_scalar_set_int(&ctx->blind, 1);
	}
	secp256k1_scalar_get_b32(nonce32, &ctx->blind);
	memcpy(keydata, nonce32, 32);
	if (seed32 != NULL) {
		memcpy(keydata + 32, seed32, 32);
	}
	secp256k1_rfc6979_hmac_sha256_initialize(&rng, keydata, seed32 ? 64 : 32);
	memset(keydata, 0, sizeof(keydata));
	/* Retry for out of range results to achieve uniformity. */
	do {
		secp256k1_rfc6979_hmac_sha256_generate(&rng, nonce32, 32);
		retry = !secp256k1_fe_set_b32(&s, nonce32);
		retry |= secp256k1_fe_is_zero(&s);
	} while (retry);
	/* Randomize the projection to defend against multiplier sidechannels. */
	secp256k1_gej_rescale(&ctx->initial, &s);
	secp256k1_fe_clear(&s);
	do {
		secp256k1_rfc6979_hmac_sha256_generate(&rng, nonce32, 32);
		secp256k1_scalar_set_b32(&b, nonce32, &retry);
		/* A blinding value of 0 works, but would undermine the projection hardening. */
		retry |= secp256k1_scalar_is_zero(&b);
	} while (retry);
	secp256k1_rfc6979_hmac_sha256_finalize(&rng);
	memset(nonce32, 0, 32);
	secp256k1_ecmult_gen(ctx, &gb, &b);
	secp256k1_scalar_negate(&b, &b);
	ctx->blind = b;
	ctx->initial = gb;
	secp256k1_scalar_clear(&b);
	secp256k1_gej_clear(&gb);
}


secp256k1_context* secp256k1_context_create(unsigned int flags) {
	secp256k1_context* ret = (secp256k1_context*)malloc(sizeof(secp256k1_context));
	ret->illegal_callback = default_illegal_callback;
	ret->error_callback = default_error_callback;

	secp256k1_ecmult_context_init(&ret->ecmult_ctx);
	secp256k1_ecmult_gen_context_init(&ret->ecmult_gen_ctx);

	if (flags & SECP256K1_CONTEXT_SIGN) {
		secp256k1_ecmult_gen_context_build(&ret->ecmult_gen_ctx, &ret->error_callback);
	}
	if (flags & SECP256K1_CONTEXT_VERIFY) {
		secp256k1_ecmult_context_build(&ret->ecmult_ctx, &ret->error_callback);
	}

	return ret;
}


/* field h */

void secp256k1_fe_verify(const secp256k1_fe *a) {
	(void)a;
}
void secp256k1_fe_normalize(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t m;
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; m = t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; m &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; m &= t3;
	VERIFY_CHECK(t4 >> 49 == 0);
	x = (t4 >> 48) | ((t4 == 0x0FFFFFFFFFFFFULL) & (m == 0xFFFFFFFFFFFFFULL) & (t0 >= 0xFFFFEFFFFFC2FULL));
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;
	VERIFY_CHECK(t4 >> 48 == x);
	t4 &= 0x0FFFFFFFFFFFFULL;
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
}

void secp256k1_fe_normalize_weak(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;
	VERIFY_CHECK(t4 >> 49 == 0);
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
}

void secp256k1_fe_normalize_var(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t m;
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; m = t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; m &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; m &= t3;
	VERIFY_CHECK(t4 >> 49 == 0);
	x = (t4 >> 48) | ((t4 == 0x0FFFFFFFFFFFFULL) & (m == 0xFFFFFFFFFFFFFULL) & (t0 >= 0xFFFFEFFFFFC2FULL));

	if (x) {
		t0 += 0x1000003D1ULL;
		t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
		t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
		t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
		t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;
		VERIFY_CHECK(t4 >> 48 == x);
		t4 &= 0x0FFFFFFFFFFFFULL;
	}

	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
}

int secp256k1_fe_normalizes_to_zero(secp256k1_fe *r) {
	uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];
	uint64_t z0, z1;
	uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
	t0 += x * 0x1000003D1ULL;
	t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL; z0 = t0; z1 = t0 ^ 0x1000003D0ULL;
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; z0 |= t1; z1 &= t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; z0 |= t2; z1 &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; z0 |= t3; z1 &= t3;
	z0 |= t4; z1 &= t4 ^ 0xF000000000000ULL;
	VERIFY_CHECK(t4 >> 49 == 0);
	return (z0 == 0) | (z1 == 0xFFFFFFFFFFFFFULL);
}

int secp256k1_fe_normalizes_to_zero_var(secp256k1_fe *r) {
	uint64_t t0, t1, t2, t3, t4;
	uint64_t z0, z1;
	uint64_t x;

	t0 = r->n[0];
	t4 = r->n[4];
	x = t4 >> 48;
	t0 += x * 0x1000003D1ULL;
	z0 = t0 & 0xFFFFFFFFFFFFFULL;
	z1 = z0 ^ 0x1000003D0ULL;
	if ((z0 != 0ULL) & (z1 != 0xFFFFFFFFFFFFFULL)) {
		return 0;
	}
	t1 = r->n[1];
	t2 = r->n[2];
	t3 = r->n[3];
	t4 &= 0x0FFFFFFFFFFFFULL;
	t1 += (t0 >> 52);
	t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; z0 |= t1; z1 &= t1;
	t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; z0 |= t2; z1 &= t2;
	t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; z0 |= t3; z1 &= t3;
	z0 |= t4; z1 &= t4 ^ 0xF000000000000ULL;
	VERIFY_CHECK(t4 >> 49 == 0);
	return (z0 == 0) | (z1 == 0xFFFFFFFFFFFFFULL);
}

 void inline secp256k1_fe_set_int(secp256k1_fe *r, int a) {
	r->n[0] = a;
	r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
}

int inline secp256k1_fe_is_zero(const secp256k1_fe *a) {
	const uint64_t *t = a->n;
	return (t[0] | t[1] | t[2] | t[3] | t[4]) == 0;
}

int inline secp256k1_fe_is_odd(const secp256k1_fe *a) {
	return a->n[0] & 1;
}

void inline secp256k1_fe_clear(secp256k1_fe *a) {
	for (int i = 0; i<5; i++)
		a->n[i] = 0;
}

int secp256k1_fe_cmp_var(const secp256k1_fe *a, const secp256k1_fe *b) {
	for (int i = 4; i >= 0; i--) {
		if (a->n[i] > b->n[i]) return 1;
		if (a->n[i] < b->n[i]) return -1;
	}
	return 0;
}

int secp256k1_fe_set_b32(secp256k1_fe *r, const unsigned char *a) {
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

/** Convert a field element to a 32-byte big endian value. Requires the input to be normalized */
void secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe *a) {
	int i;
	for (i = 0; i<32; i++) {
		int j;
		int c = 0;
		for (j = 0; j<2; j++) {
			int limb = (8 * i + 4 * j) / 52;
			int shift = (8 * i + 4 * j) % 52;
			c |= ((a->n[limb] >> shift) & 0xF) << (4 * j);
		}
		r[31 - i] = c;
	}
}

void inline secp256k1_fe_negate(secp256k1_fe *r, const secp256k1_fe *a, int m) {
	r->n[0] = 0xFFFFEFFFFFC2FULL * 2 * (m + 1) - a->n[0];
	r->n[1] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[1];
	r->n[2] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[2];
	r->n[3] = 0xFFFFFFFFFFFFFULL * 2 * (m + 1) - a->n[3];
	r->n[4] = 0x0FFFFFFFFFFFFULL * 2 * (m + 1) - a->n[4];
}

void inline secp256k1_fe_mul_int(secp256k1_fe *r, int a) {
	r->n[0] *= a;
	r->n[1] *= a;
	r->n[2] *= a;
	r->n[3] *= a;
	r->n[4] *= a;
}

void inline secp256k1_fe_add(secp256k1_fe *r, const secp256k1_fe *a) {
	r->n[0] += a->n[0];
	r->n[1] += a->n[1];
	r->n[2] += a->n[2];
	r->n[3] += a->n[3];
	r->n[4] += a->n[4];
}

void secp256k1_fe_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b) {
	secp256k1_fe_mul_inner(r->n, a->n, b->n);
}

void secp256k1_fe_sqr(secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe_sqr_inner(r->n, a->n);
}

void inline secp256k1_fe_cmov(secp256k1_fe *r, const secp256k1_fe *a, int flag) {
	uint64_t mask0, mask1;
	mask0 = flag + ~((uint64_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
	r->n[4] = (r->n[4] & mask0) | (a->n[4] & mask1);
}

void inline secp256k1_fe_storage_cmov(secp256k1_fe_storage *r, const secp256k1_fe_storage *a, int flag) {
	uint64_t mask0, mask1;
	mask0 = flag + ~((uint64_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
}

void secp256k1_fe_to_storage(secp256k1_fe_storage *r, const secp256k1_fe *a) {
	r->n[0] = a->n[0] | a->n[1] << 52;
	r->n[1] = a->n[1] >> 12 | a->n[2] << 40;
	r->n[2] = a->n[2] >> 24 | a->n[3] << 28;
	r->n[3] = a->n[3] >> 36 | a->n[4] << 16;
}

void inline secp256k1_fe_from_storage(secp256k1_fe *r, const secp256k1_fe_storage *a) {
	r->n[0] = a->n[0] & 0xFFFFFFFFFFFFFULL;
	r->n[1] = a->n[0] >> 52 | ((a->n[1] << 12) & 0xFFFFFFFFFFFFFULL);
	r->n[2] = a->n[1] >> 40 | ((a->n[2] << 24) & 0xFFFFFFFFFFFFFULL);
	r->n[3] = a->n[2] >> 28 | ((a->n[3] << 36) & 0xFFFFFFFFFFFFFULL);
	r->n[4] = a->n[3] >> 16;
}

/* group h */

static void secp256k1_ge_set_gej_zinv(secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_fe *zi) {
	secp256k1_fe zi2;
	secp256k1_fe zi3;
	secp256k1_fe_sqr(&zi2, zi);
	secp256k1_fe_mul(&zi3, &zi2, zi);
	secp256k1_fe_mul(&r->x, &a->x, &zi2);
	secp256k1_fe_mul(&r->y, &a->y, &zi3);
	r->infinity = a->infinity;
}

static void secp256k1_ge_set_xy(secp256k1_ge *r, const secp256k1_fe *x, const secp256k1_fe *y) {
	r->infinity = 0;
	r->x = *x;
	r->y = *y;
}

static int secp256k1_ge_is_infinity(const secp256k1_ge *a) {
	return a->infinity;
}

static void secp256k1_ge_neg(secp256k1_ge *r, const secp256k1_ge *a) {
	*r = *a;
	secp256k1_fe_normalize_weak(&r->y);
	secp256k1_fe_negate(&r->y, &r->y, 1);
}

static void secp256k1_ge_set_gej(secp256k1_ge *r, secp256k1_gej *a) {
	secp256k1_fe z2, z3;
	r->infinity = a->infinity;
	secp256k1_fe_inv(&a->z, &a->z);
	secp256k1_fe_sqr(&z2, &a->z);
	secp256k1_fe_mul(&z3, &a->z, &z2);
	secp256k1_fe_mul(&a->x, &a->x, &z2);
	secp256k1_fe_mul(&a->y, &a->y, &z3);
	secp256k1_fe_set_int(&a->z, 1);
	r->x = a->x;
	r->y = a->y;
}

void secp256k1_ge_set_gej_var(secp256k1_ge *r, secp256k1_gej *a) {
	secp256k1_fe z2, z3;
	r->infinity = a->infinity;
	if (a->infinity) {
		return;
	}
	secp256k1_fe_inv_var(&a->z, &a->z);
	secp256k1_fe_sqr(&z2, &a->z);
	secp256k1_fe_mul(&z3, &a->z, &z2);
	secp256k1_fe_mul(&a->x, &a->x, &z2);
	secp256k1_fe_mul(&a->y, &a->y, &z3);
	secp256k1_fe_set_int(&a->z, 1);
	r->x = a->x;
	r->y = a->y;
}

void secp256k1_ge_set_all_gej_var(size_t len, secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_callback *cb) {
	secp256k1_fe *az;
	secp256k1_fe *azi;
	size_t i;
	size_t count = 0;
	az = (secp256k1_fe *)checked_malloc(cb, sizeof(secp256k1_fe) * len);
	for (i = 0; i < len; i++) {
		if (!a[i].infinity) {
			az[count++] = a[i].z;
		}
	}

	azi = (secp256k1_fe *)checked_malloc(cb, sizeof(secp256k1_fe) * count);
	secp256k1_fe_inv_all_var(count, azi, az);
	free(az);

	count = 0;
	for (i = 0; i < len; i++) {
		r[i].infinity = a[i].infinity;
		if (!a[i].infinity) {
			secp256k1_ge_set_gej_zinv(&r[i], &a[i], &azi[count++]);
		}
	}
	free(azi);
}

void secp256k1_ge_set_table_gej_var(size_t len, secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_fe *zr) {
	size_t i = len - 1;
	secp256k1_fe zi;

	if (len > 0) {
		secp256k1_fe_inv(&zi, &a[i].z);
		secp256k1_ge_set_gej_zinv(&r[i], &a[i], &zi);
		while (i > 0) {
			secp256k1_fe_mul(&zi, &zi, &zr[i]);
			i--;
			secp256k1_ge_set_gej_zinv(&r[i], &a[i], &zi);
		}
	}
}

void secp256k1_ge_globalz_set_table_gej(size_t len, secp256k1_ge *r, secp256k1_fe *globalz, const secp256k1_gej *a, const secp256k1_fe *zr) {
	size_t i = len - 1;
	secp256k1_fe zs;

	if (len > 0) {
		/* The z of the final point gives us the "global Z" for the table. */
		r[i].x = a[i].x;
		r[i].y = a[i].y;
		*globalz = a[i].z;
		r[i].infinity = 0;
		zs = zr[i];
		while (i > 0) {
			if (i != len - 1) {
				secp256k1_fe_mul(&zs, &zs, &zr[i]);
			}
			i--;
			secp256k1_ge_set_gej_zinv(&r[i], &a[i], &zs);
		}
	}
}

void secp256k1_gej_set_infinity(secp256k1_gej *r) {
	r->infinity = 1;
	secp256k1_fe_set_int(&r->x, 0);
	secp256k1_fe_set_int(&r->y, 0);
	secp256k1_fe_set_int(&r->z, 0);
}

void secp256k1_gej_clear(secp256k1_gej *r) {
	r->infinity = 0;
	secp256k1_fe_clear(&r->x);
	secp256k1_fe_clear(&r->y);
	secp256k1_fe_clear(&r->z);
}

void secp256k1_ge_clear(secp256k1_ge *r) {
	r->infinity = 0;
	secp256k1_fe_clear(&r->x);
	secp256k1_fe_clear(&r->y);
}

int secp256k1_ge_set_xo_var(secp256k1_ge *r, const secp256k1_fe *x, int odd) {
	secp256k1_fe x2, x3, c;
	r->x = *x;
	secp256k1_fe_sqr(&x2, x);
	secp256k1_fe_mul(&x3, x, &x2);
	r->infinity = 0;
	secp256k1_fe_set_int(&c, 7);
	secp256k1_fe_add(&c, &x3);
	if (!secp256k1_fe_sqrt_var(&r->y, &c)) {
		return 0;
	}
	secp256k1_fe_normalize_var(&r->y);
	if (secp256k1_fe_is_odd(&r->y) != odd) {
		secp256k1_fe_negate(&r->y, &r->y, 1);
	}
	return 1;
}

void secp256k1_gej_set_ge(secp256k1_gej *r, const secp256k1_ge *a) {
	r->infinity = a->infinity;
	r->x = a->x;
	r->y = a->y;
	secp256k1_fe_set_int(&r->z, 1);
}

int secp256k1_gej_eq_x_var(const secp256k1_fe *x, const secp256k1_gej *a) {
	secp256k1_fe r, r2;
	VERIFY_CHECK(!a->infinity);
	secp256k1_fe_sqr(&r, &a->z); secp256k1_fe_mul(&r, &r, x);
	r2 = a->x; secp256k1_fe_normalize_weak(&r2);
	return secp256k1_fe_equal_var(&r, &r2);
}

void secp256k1_gej_neg(secp256k1_gej *r, const secp256k1_gej *a) {
	r->infinity = a->infinity;
	r->x = a->x;
	r->y = a->y;
	r->z = a->z;
	secp256k1_fe_normalize_weak(&r->y);
	secp256k1_fe_negate(&r->y, &r->y, 1);
}

int secp256k1_gej_is_infinity(const secp256k1_gej *a) {
	return a->infinity;
}

int secp256k1_gej_is_valid_var(const secp256k1_gej *a) {
	secp256k1_fe y2, x3, z2, z6;
	if (a->infinity) {
		return 0;
	}
	secp256k1_fe_sqr(&y2, &a->y);
	secp256k1_fe_sqr(&x3, &a->x); secp256k1_fe_mul(&x3, &x3, &a->x);
	secp256k1_fe_sqr(&z2, &a->z);
	secp256k1_fe_sqr(&z6, &z2); secp256k1_fe_mul(&z6, &z6, &z2);
	secp256k1_fe_mul_int(&z6, 7);
	secp256k1_fe_add(&x3, &z6);
	secp256k1_fe_normalize_weak(&x3);
	return secp256k1_fe_equal_var(&y2, &x3);
}

int secp256k1_ge_is_valid_var(const secp256k1_ge *a) {
	secp256k1_fe y2, x3, c;
	if (a->infinity) {
		return 0;
	}
	/* y^2 = x^3 + 7 */
	secp256k1_fe_sqr(&y2, &a->y);
	secp256k1_fe_sqr(&x3, &a->x); secp256k1_fe_mul(&x3, &x3, &a->x);
	secp256k1_fe_set_int(&c, 7);
	secp256k1_fe_add(&x3, &c);
	secp256k1_fe_normalize_weak(&x3);
	return secp256k1_fe_equal_var(&y2, &x3);
}

void secp256k1_gej_double_var(secp256k1_gej *r, const secp256k1_gej *a, secp256k1_fe *rzr) {
	secp256k1_fe t1, t2, t3, t4;
	r->infinity = a->infinity;
	if (r->infinity) {
		if (rzr != NULL) {
			secp256k1_fe_set_int(rzr, 1);
		}
		return;
	}
	if (rzr != NULL) {
		*rzr = a->y;
		secp256k1_fe_normalize_weak(rzr);
		secp256k1_fe_mul_int(rzr, 2);
	}
	secp256k1_fe_mul(&r->z, &a->z, &a->y);
	secp256k1_fe_mul_int(&r->z, 2);       /* Z' = 2*Y*Z (2) */
	secp256k1_fe_sqr(&t1, &a->x);
	secp256k1_fe_mul_int(&t1, 3);         /* T1 = 3*X^2 (3) */
	secp256k1_fe_sqr(&t2, &t1);           /* T2 = 9*X^4 (1) */
	secp256k1_fe_sqr(&t3, &a->y);
	secp256k1_fe_mul_int(&t3, 2);         /* T3 = 2*Y^2 (2) */
	secp256k1_fe_sqr(&t4, &t3);
	secp256k1_fe_mul_int(&t4, 2);         /* T4 = 8*Y^4 (2) */
	secp256k1_fe_mul(&t3, &t3, &a->x);    /* T3 = 2*X*Y^2 (1) */
	r->x = t3;
	secp256k1_fe_mul_int(&r->x, 4);       /* X' = 8*X*Y^2 (4) */
	secp256k1_fe_negate(&r->x, &r->x, 4); /* X' = -8*X*Y^2 (5) */
	secp256k1_fe_add(&r->x, &t2);         /* X' = 9*X^4 - 8*X*Y^2 (6) */
	secp256k1_fe_negate(&t2, &t2, 1);     /* T2 = -9*X^4 (2) */
	secp256k1_fe_mul_int(&t3, 6);         /* T3 = 12*X*Y^2 (6) */
	secp256k1_fe_add(&t3, &t2);           /* T3 = 12*X*Y^2 - 9*X^4 (8) */
	secp256k1_fe_mul(&r->y, &t1, &t3);    /* Y' = 36*X^3*Y^2 - 27*X^6 (1) */
	secp256k1_fe_negate(&t2, &t4, 2);     /* T2 = -8*Y^4 (3) */
	secp256k1_fe_add(&r->y, &t2);         /* Y' = 36*X^3*Y^2 - 27*X^6 - 8*Y^4 (4) */
}

void secp256k1_gej_double_nonzero(secp256k1_gej *r, const secp256k1_gej *a, secp256k1_fe *rzr) {
	VERIFY_CHECK(!secp256k1_gej_is_infinity(a));
	secp256k1_gej_double_var(r, a, rzr);
}

void secp256k1_gej_add_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_gej *b, secp256k1_fe *rzr) {
	secp256k1_fe z22, z12, u1, u2, s1, s2, h, i, i2, h2, h3, t;

	if (a->infinity) {
		VERIFY_CHECK(rzr == NULL);
		*r = *b;
		return;
	}

	if (b->infinity) {
		if (rzr != NULL) {
			secp256k1_fe_set_int(rzr, 1);
		}
		*r = *a;
		return;
	}

	r->infinity = 0;
	secp256k1_fe_sqr(&z22, &b->z);
	secp256k1_fe_sqr(&z12, &a->z);
	secp256k1_fe_mul(&u1, &a->x, &z22);
	secp256k1_fe_mul(&u2, &b->x, &z12);
	secp256k1_fe_mul(&s1, &a->y, &z22); secp256k1_fe_mul(&s1, &s1, &b->z);
	secp256k1_fe_mul(&s2, &b->y, &z12); secp256k1_fe_mul(&s2, &s2, &a->z);
	secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
	secp256k1_fe_negate(&i, &s1, 1); secp256k1_fe_add(&i, &s2);
	if (secp256k1_fe_normalizes_to_zero_var(&h)) {
		if (secp256k1_fe_normalizes_to_zero_var(&i)) {
			secp256k1_gej_double_var(r, a, rzr);
		}
		else {
			if (rzr != NULL) {
				secp256k1_fe_set_int(rzr, 0);
			}
			r->infinity = 1;
		}
		return;
	}
	secp256k1_fe_sqr(&i2, &i);
	secp256k1_fe_sqr(&h2, &h);
	secp256k1_fe_mul(&h3, &h, &h2);
	secp256k1_fe_mul(&h, &h, &b->z);
	if (rzr != NULL) {
		*rzr = h;
	}
	secp256k1_fe_mul(&r->z, &a->z, &h);
	secp256k1_fe_mul(&t, &u1, &h2);
	r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3); secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
	secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t); secp256k1_fe_mul(&r->y, &r->y, &i);
	secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
	secp256k1_fe_add(&r->y, &h3);
}

void secp256k1_gej_add_ge_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b, secp256k1_fe *rzr) {
	secp256k1_fe z12, u1, u2, s1, s2, h, i, i2, h2, h3, t;
	secp256k1_fe_sqr(&z12, &a->z);
	u1 = a->x; secp256k1_fe_normalize_weak(&u1);
	secp256k1_fe_mul(&u2, &b->x, &z12);
	s1 = a->y; secp256k1_fe_normalize_weak(&s1);
	secp256k1_fe_mul(&s2, &b->y, &z12); secp256k1_fe_mul(&s2, &s2, &a->z);
	secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
	secp256k1_fe_negate(&i, &s1, 1); secp256k1_fe_add(&i, &s2);
	secp256k1_fe_sqr(&i2, &i);
	secp256k1_fe_sqr(&h2, &h);
	secp256k1_fe_mul(&h3, &h, &h2);

	secp256k1_fe_mul(&r->z, &a->z, &h);
	secp256k1_fe_mul(&t, &u1, &h2);
	r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3); secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
	secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t); secp256k1_fe_mul(&r->y, &r->y, &i);
	secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
	secp256k1_fe_add(&r->y, &h3);
}

void secp256k1_gej_add_zinv_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b, const secp256k1_fe *bzinv) {
	secp256k1_fe az, z12, u1, u2, s1, s2, h, i, i2, h2, h3, t;
	if (b->infinity) {
		*r = *a;
		return;
	}
	if (a->infinity) {
		secp256k1_fe bzinv2, bzinv3;
		r->infinity = b->infinity;
		secp256k1_fe_sqr(&bzinv2, bzinv);
		secp256k1_fe_mul(&bzinv3, &bzinv2, bzinv);
		secp256k1_fe_mul(&r->x, &b->x, &bzinv2);
		secp256k1_fe_mul(&r->y, &b->y, &bzinv3);
		secp256k1_fe_set_int(&r->z, 1);
		return;
	}
	r->infinity = 0;
	secp256k1_fe_mul(&az, &a->z, bzinv);
	secp256k1_fe_sqr(&z12, &az);
	u1 = a->x; secp256k1_fe_normalize_weak(&u1);
	secp256k1_fe_mul(&u2, &b->x, &z12);
	s1 = a->y; secp256k1_fe_normalize_weak(&s1);
	secp256k1_fe_mul(&s2, &b->y, &z12); secp256k1_fe_mul(&s2, &s2, &az);
	secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
	secp256k1_fe_negate(&i, &s1, 1); secp256k1_fe_add(&i, &s2);
	if (secp256k1_fe_normalizes_to_zero_var(&h)) {
		if (secp256k1_fe_normalizes_to_zero_var(&i)) {
			secp256k1_gej_double_var(r, a, NULL);
		}
		else {
			r->infinity = 1;
		}
		return;
	}
	secp256k1_fe_sqr(&i2, &i);
	secp256k1_fe_sqr(&h2, &h);
	secp256k1_fe_mul(&h3, &h, &h2);
	r->z = a->z; secp256k1_fe_mul(&r->z, &r->z, &h);
	secp256k1_fe_mul(&t, &u1, &h2);
	r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3); secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
	secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t); secp256k1_fe_mul(&r->y, &r->y, &i);
	secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
	secp256k1_fe_add(&r->y, &h3);
}


void secp256k1_gej_add_ge(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b) {
	static const secp256k1_fe fe_1 = SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 1);
	secp256k1_fe zz, u1, u2, s1, s2, t, tt, m, n, q, rr;
	secp256k1_fe m_alt, rr_alt;
	int infinity, degenerate;
	VERIFY_CHECK(!b->infinity);
	VERIFY_CHECK(a->infinity == 0 || a->infinity == 1);
	secp256k1_fe_sqr(&zz, &a->z);                       /* z = Z1^2 */
	u1 = a->x; secp256k1_fe_normalize_weak(&u1);        /* u1 = U1 = X1*Z2^2 (1) */
	secp256k1_fe_mul(&u2, &b->x, &zz);                  /* u2 = U2 = X2*Z1^2 (1) */
	s1 = a->y; secp256k1_fe_normalize_weak(&s1);        /* s1 = S1 = Y1*Z2^3 (1) */
	secp256k1_fe_mul(&s2, &b->y, &zz);                  /* s2 = Y2*Z1^2 (1) */
	secp256k1_fe_mul(&s2, &s2, &a->z);                  /* s2 = S2 = Y2*Z1^3 (1) */
	t = u1; secp256k1_fe_add(&t, &u2);                  /* t = T = U1+U2 (2) */
	m = s1; secp256k1_fe_add(&m, &s2);                  /* m = M = S1+S2 (2) */
	secp256k1_fe_sqr(&rr, &t);                          /* rr = T^2 (1) */
	secp256k1_fe_negate(&m_alt, &u2, 1);                /* Malt = -X2*Z1^2 */
	secp256k1_fe_mul(&tt, &u1, &m_alt);                 /* tt = -U1*U2 (2) */
	secp256k1_fe_add(&rr, &tt);                         /* rr = R = T^2-U1*U2 (3) */
	degenerate = secp256k1_fe_normalizes_to_zero(&m) &
		secp256k1_fe_normalizes_to_zero(&rr);
	rr_alt = s1;
	secp256k1_fe_mul_int(&rr_alt, 2);       /* rr = Y1*Z2^3 - Y2*Z1^3 (2) */
	secp256k1_fe_add(&m_alt, &u1);          /* Malt = X1*Z2^2 - X2*Z1^2 */
	secp256k1_fe_cmov(&rr_alt, &rr, !degenerate);
	secp256k1_fe_cmov(&m_alt, &m, !degenerate);
	secp256k1_fe_sqr(&n, &m_alt);                       /* n = Malt^2 (1) */
	secp256k1_fe_mul(&q, &n, &t);                       /* q = Q = T*Malt^2 (1) */
	secp256k1_fe_sqr(&n, &n);
	secp256k1_fe_cmov(&n, &m, degenerate);              /* n = M^3 * Malt (2) */
	secp256k1_fe_sqr(&t, &rr_alt);                      /* t = Ralt^2 (1) */
	secp256k1_fe_mul(&r->z, &a->z, &m_alt);             /* r->z = Malt*Z (1) */
	infinity = secp256k1_fe_normalizes_to_zero(&r->z) * (1 - a->infinity);
	secp256k1_fe_mul_int(&r->z, 2);                     /* r->z = Z3 = 2*Malt*Z (2) */
	secp256k1_fe_negate(&q, &q, 1);                     /* q = -Q (2) */
	secp256k1_fe_add(&t, &q);                           /* t = Ralt^2-Q (3) */
	secp256k1_fe_normalize_weak(&t);
	r->x = t;                                           /* r->x = Ralt^2-Q (1) */
	secp256k1_fe_mul_int(&t, 2);                        /* t = 2*x3 (2) */
	secp256k1_fe_add(&t, &q);                           /* t = 2*x3 - Q: (4) */
	secp256k1_fe_mul(&t, &t, &rr_alt);                  /* t = Ralt*(2*x3 - Q) (1) */
	secp256k1_fe_add(&t, &n);                           /* t = Ralt*(2*x3 - Q) + M^3*Malt (3) */
	secp256k1_fe_negate(&r->y, &t, 3);                  /* r->y = Ralt*(Q - 2x3) - M^3*Malt (4) */
	secp256k1_fe_normalize_weak(&r->y);
	secp256k1_fe_mul_int(&r->x, 4);                     /* r->x = X3 = 4*(Ralt^2-Q) */
	secp256k1_fe_mul_int(&r->y, 4);                     /* r->y = Y3 = 4*Ralt*(Q - 2x3) - 4*M^3*Malt (4) */
	secp256k1_fe_cmov(&r->x, &b->x, a->infinity);
	secp256k1_fe_cmov(&r->y, &b->y, a->infinity);
	secp256k1_fe_cmov(&r->z, &fe_1, a->infinity);
	r->infinity = infinity;
}

void secp256k1_gej_rescale(secp256k1_gej *r, const secp256k1_fe *s) {
	/* Operations: 4 mul, 1 sqr */
	secp256k1_fe zz;
	VERIFY_CHECK(!secp256k1_fe_is_zero(s));
	secp256k1_fe_sqr(&zz, s);
	secp256k1_fe_mul(&r->x, &r->x, &zz);                /* r->x *= s^2 */
	secp256k1_fe_mul(&r->y, &r->y, &zz);
	secp256k1_fe_mul(&r->y, &r->y, s);                  /* r->y *= s^3 */
	secp256k1_fe_mul(&r->z, &r->z, s);                  /* r->z *= s   */
}

void secp256k1_ge_to_storage(secp256k1_ge_storage *r, const secp256k1_ge *a) {
	secp256k1_fe x, y;
	VERIFY_CHECK(!a->infinity);
	x = a->x;
	secp256k1_fe_normalize(&x);
	y = a->y;
	secp256k1_fe_normalize(&y);
	secp256k1_fe_to_storage(&r->x, &x);
	secp256k1_fe_to_storage(&r->y, &y);
}

void secp256k1_ge_from_storage(secp256k1_ge *r, const secp256k1_ge_storage *a) {
	secp256k1_fe_from_storage(&r->x, &a->x);
	secp256k1_fe_from_storage(&r->y, &a->y);
	r->infinity = 0;
}

void secp256k1_ge_storage_cmov(secp256k1_ge_storage *r, const secp256k1_ge_storage *a, int flag) {
	secp256k1_fe_storage_cmov(&r->x, &a->x, flag);
	secp256k1_fe_storage_cmov(&r->y, &a->y, flag);
}


/* scalar h */
void secp256k1_scalar_clear(secp256k1_scalar *r) {
	unsigned int *p = (unsigned int *)r->d;
	p[0] = 0; p[1] = 0; p[2] = 0; p[3] = 0; p[4] = 0; p[5] = 0; p[6] = 0; p[7] = 0;
}

void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v) {
	unsigned int *p = (unsigned int *)r->d;
	p[0] = v; p[1] = 0; p[2] = 0; p[3] = 0; p[4] = 0; p[5] = 0; p[6] = 0; p[7] = 0;
}
unsigned int secp256k1_scalar_get_bits(const secp256k1_scalar *a, unsigned int offset, unsigned int count) {
	VERIFY_CHECK((offset + count - 1) >> 5 == offset >> 5);
	return (a->d[offset >> 5] >> (offset & 0x1F)) & ((1 << count) - 1);
}

unsigned int secp256k1_scalar_get_bits_var(const secp256k1_scalar *a, unsigned int offset, unsigned int count) {
	VERIFY_CHECK(count < 32);
	VERIFY_CHECK(offset + count <= 256);
	if ((offset + count - 1) >> 5 == offset >> 5) {
		return secp256k1_scalar_get_bits(a, offset, count);
	}
	else {
		VERIFY_CHECK((offset >> 5) + 1 < 8);
		return ((a->d[offset >> 5] >> (offset & 0x1F)) | (a->d[(offset >> 5) + 1] << (32 - (offset & 0x1F)))) & ((((uint32_t)1) << count) - 1);
	}
}

int secp256k1_scalar_check_overflow(const secp256k1_scalar *a) {
	int yes = 0;
	int no = 0;
	no |= (a->d[7] < SECP256K1_N_7); /* No need for a > check. */
	no |= (a->d[6] < SECP256K1_N_6); /* No need for a > check. */
	no |= (a->d[5] < SECP256K1_N_5); /* No need for a > check. */
	no |= (a->d[4] < SECP256K1_N_4);
	yes |= (a->d[4] > SECP256K1_N_4) & ~no;
	no |= (a->d[3] < SECP256K1_N_3) & ~yes;
	yes |= (a->d[3] > SECP256K1_N_3) & ~no;
	no |= (a->d[2] < SECP256K1_N_2) & ~yes;
	yes |= (a->d[2] > SECP256K1_N_2) & ~no;
	no |= (a->d[1] < SECP256K1_N_1) & ~yes;
	yes |= (a->d[1] > SECP256K1_N_1) & ~no;
	yes |= (a->d[0] >= SECP256K1_N_0) & ~no;
	return yes;
}

int secp256k1_scalar_reduce(secp256k1_scalar *r, uint32_t overflow) {
	uint64_t t;
	VERIFY_CHECK(overflow <= 1);
	t = (uint64_t)r->d[0] + overflow * SECP256K1_N_C_0;
	r->d[0] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[1] + overflow * SECP256K1_N_C_1;
	r->d[1] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[2] + overflow * SECP256K1_N_C_2;
	r->d[2] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[3] + overflow * SECP256K1_N_C_3;
	r->d[3] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[4] + overflow * SECP256K1_N_C_4;
	r->d[4] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[5];
	r->d[5] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[6];
	r->d[6] = t & 0xFFFFFFFFUL; t >>= 32;
	t += (uint64_t)r->d[7];
	r->d[7] = t & 0xFFFFFFFFUL;
	return overflow;
}

int secp256k1_scalar_add(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b) {
	int overflow;
	uint64_t t = (uint64_t)a->d[0] + b->d[0];
	r->d[0] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[1] + b->d[1];
	r->d[1] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[2] + b->d[2];
	r->d[2] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[3] + b->d[3];
	r->d[3] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[4] + b->d[4];
	r->d[4] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[5] + b->d[5];
	r->d[5] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[6] + b->d[6];
	r->d[6] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)a->d[7] + b->d[7];
	r->d[7] = t & 0xFFFFFFFFULL; t >>= 32;
	overflow = t + secp256k1_scalar_check_overflow(r);
	VERIFY_CHECK(overflow == 0 || overflow == 1);
	secp256k1_scalar_reduce(r, overflow);
	return overflow;
}

void secp256k1_scalar_cadd_bit(secp256k1_scalar *r, unsigned int bit, int flag) {
	uint64_t t;
	VERIFY_CHECK(bit < 256);
	bit += ((uint32_t)flag - 1) & 0x100;  /* forcing (bit >> 5) > 7 makes this a noop */
	t = (uint64_t)r->d[0] + (((uint32_t)((bit >> 5) == 0)) << (bit & 0x1F));
	r->d[0] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[1] + (((uint32_t)((bit >> 5) == 1)) << (bit & 0x1F));
	r->d[1] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[2] + (((uint32_t)((bit >> 5) == 2)) << (bit & 0x1F));
	r->d[2] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[3] + (((uint32_t)((bit >> 5) == 3)) << (bit & 0x1F));
	r->d[3] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[4] + (((uint32_t)((bit >> 5) == 4)) << (bit & 0x1F));
	r->d[4] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[5] + (((uint32_t)((bit >> 5) == 5)) << (bit & 0x1F));
	r->d[5] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[6] + (((uint32_t)((bit >> 5) == 6)) << (bit & 0x1F));
	r->d[6] = t & 0xFFFFFFFFULL; t >>= 32;
	t += (uint64_t)r->d[7] + (((uint32_t)((bit >> 5) == 7)) << (bit & 0x1F));
	r->d[7] = t & 0xFFFFFFFFULL;
}

void secp256k1_scalar_set_b32(secp256k1_scalar *r, const unsigned char *b32, int *overflow) {
	int over;
	r->d[0] = (uint32_t)b32[31] | (uint32_t)b32[30] << 8 | (uint32_t)b32[29] << 16 | (uint32_t)b32[28] << 24;
	r->d[1] = (uint32_t)b32[27] | (uint32_t)b32[26] << 8 | (uint32_t)b32[25] << 16 | (uint32_t)b32[24] << 24;
	r->d[2] = (uint32_t)b32[23] | (uint32_t)b32[22] << 8 | (uint32_t)b32[21] << 16 | (uint32_t)b32[20] << 24;
	r->d[3] = (uint32_t)b32[19] | (uint32_t)b32[18] << 8 | (uint32_t)b32[17] << 16 | (uint32_t)b32[16] << 24;
	r->d[4] = (uint32_t)b32[15] | (uint32_t)b32[14] << 8 | (uint32_t)b32[13] << 16 | (uint32_t)b32[12] << 24;
	r->d[5] = (uint32_t)b32[11] | (uint32_t)b32[10] << 8 | (uint32_t)b32[9] << 16 | (uint32_t)b32[8] << 24;
	r->d[6] = (uint32_t)b32[7] | (uint32_t)b32[6] << 8 | (uint32_t)b32[5] << 16 | (uint32_t)b32[4] << 24;
	r->d[7] = (uint32_t)b32[3] | (uint32_t)b32[2] << 8 | (uint32_t)b32[1] << 16 | (uint32_t)b32[0] << 24;
	over = secp256k1_scalar_reduce(r, secp256k1_scalar_check_overflow(r));
	if (overflow) {
		*overflow = over;
	}
}

void secp256k1_scalar_get_b32(unsigned char *bin, const secp256k1_scalar* a) {
	bin[0] = a->d[7] >> 24; bin[1] = a->d[7] >> 16; bin[2] = a->d[7] >> 8; bin[3] = a->d[7];
	bin[4] = a->d[6] >> 24; bin[5] = a->d[6] >> 16; bin[6] = a->d[6] >> 8; bin[7] = a->d[6];
	bin[8] = a->d[5] >> 24; bin[9] = a->d[5] >> 16; bin[10] = a->d[5] >> 8; bin[11] = a->d[5];
	bin[12] = a->d[4] >> 24; bin[13] = a->d[4] >> 16; bin[14] = a->d[4] >> 8; bin[15] = a->d[4];
	bin[16] = a->d[3] >> 24; bin[17] = a->d[3] >> 16; bin[18] = a->d[3] >> 8; bin[19] = a->d[3];
	bin[20] = a->d[2] >> 24; bin[21] = a->d[2] >> 16; bin[22] = a->d[2] >> 8; bin[23] = a->d[2];
	bin[24] = a->d[1] >> 24; bin[25] = a->d[1] >> 16; bin[26] = a->d[1] >> 8; bin[27] = a->d[1];
	bin[28] = a->d[0] >> 24; bin[29] = a->d[0] >> 16; bin[30] = a->d[0] >> 8; bin[31] = a->d[0];
}

int secp256k1_scalar_is_zero(const secp256k1_scalar *a) {
	return (a->d[0] | a->d[1] | a->d[2] | a->d[3] | a->d[4] | a->d[5] | a->d[6] | a->d[7]) == 0;
}

void secp256k1_scalar_negate(secp256k1_scalar *r, const secp256k1_scalar *a) {
	uint32_t nonzero = 0xFFFFFFFFUL * (secp256k1_scalar_is_zero(a) == 0);
	uint64_t t = (uint64_t)(~a->d[0]) + SECP256K1_N_0 + 1;
	r->d[0] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[1]) + SECP256K1_N_1;
	r->d[1] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[2]) + SECP256K1_N_2;
	r->d[2] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[3]) + SECP256K1_N_3;
	r->d[3] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[4]) + SECP256K1_N_4;
	r->d[4] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[5]) + SECP256K1_N_5;
	r->d[5] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[6]) + SECP256K1_N_6;
	r->d[6] = t & nonzero; t >>= 32;
	t += (uint64_t)(~a->d[7]) + SECP256K1_N_7;
	r->d[7] = t & nonzero;
}

int secp256k1_scalar_is_one(const secp256k1_scalar *a) {
	return ((a->d[0] ^ 1) | a->d[1] | a->d[2] | a->d[3] | a->d[4] | a->d[5] | a->d[6] | a->d[7]) == 0;
}

int secp256k1_scalar_is_high(const secp256k1_scalar *a) {
	int yes = 0;
	int no = 0;
	no |= (a->d[7] < SECP256K1_N_H_7);
	yes |= (a->d[7] > SECP256K1_N_H_7) & ~no;
	no |= (a->d[6] < SECP256K1_N_H_6) & ~yes; /* No need for a > check. */
	no |= (a->d[5] < SECP256K1_N_H_5) & ~yes; /* No need for a > check. */
	no |= (a->d[4] < SECP256K1_N_H_4) & ~yes; /* No need for a > check. */
	no |= (a->d[3] < SECP256K1_N_H_3) & ~yes;
	yes |= (a->d[3] > SECP256K1_N_H_3) & ~no;
	no |= (a->d[2] < SECP256K1_N_H_2) & ~yes;
	yes |= (a->d[2] > SECP256K1_N_H_2) & ~no;
	no |= (a->d[1] < SECP256K1_N_H_1) & ~yes;
	yes |= (a->d[1] > SECP256K1_N_H_1) & ~no;
	yes |= (a->d[0] > SECP256K1_N_H_0) & ~no;
	return yes;
}

int secp256k1_scalar_cond_negate(secp256k1_scalar *r, int flag) {
	/* If we are flag = 0, mask = 00...00 and this is a no-op;
	* if we are flag = 1, mask = 11...11 and this is identical to secp256k1_scalar_negate */
	uint32_t mask = !flag - 1;
	uint32_t nonzero = 0xFFFFFFFFUL * (secp256k1_scalar_is_zero(r) == 0);
	uint64_t t = (uint64_t)(r->d[0] ^ mask) + ((SECP256K1_N_0 + 1) & mask);
	r->d[0] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[1] ^ mask) + (SECP256K1_N_1 & mask);
	r->d[1] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[2] ^ mask) + (SECP256K1_N_2 & mask);
	r->d[2] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[3] ^ mask) + (SECP256K1_N_3 & mask);
	r->d[3] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[4] ^ mask) + (SECP256K1_N_4 & mask);
	r->d[4] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[5] ^ mask) + (SECP256K1_N_5 & mask);
	r->d[5] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[6] ^ mask) + (SECP256K1_N_6 & mask);
	r->d[6] = t & nonzero; t >>= 32;
	t += (uint64_t)(r->d[7] ^ mask) + (SECP256K1_N_7 & mask);
	r->d[7] = t & nonzero;
	return 2 * (mask == 0) - 1;
}

void secp256k1_scalar_reduce_512(secp256k1_scalar *r, const uint32_t *l) {
	uint64_t c;
	uint32_t n0 = l[8], n1 = l[9], n2 = l[10], n3 = l[11], n4 = l[12], n5 = l[13], n6 = l[14], n7 = l[15];
	uint32_t m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12;
	uint32_t p0, p1, p2, p3, p4, p5, p6, p7, p8;

	/* 96 bit accumulator. */
	uint32_t c0, c1, c2;

	c0 = l[0]; c1 = 0; c2 = 0;
	muladd_fast(n0, SECP256K1_N_C_0);
	extract_fast(m0);
	sumadd_fast(l[1]);
	muladd(n1, SECP256K1_N_C_0);
	muladd(n0, SECP256K1_N_C_1);
	extract(m1);
	sumadd(l[2]);
	muladd(n2, SECP256K1_N_C_0);
	muladd(n1, SECP256K1_N_C_1);
	muladd(n0, SECP256K1_N_C_2);
	extract(m2);
	sumadd(l[3]);
	muladd(n3, SECP256K1_N_C_0);
	muladd(n2, SECP256K1_N_C_1);
	muladd(n1, SECP256K1_N_C_2);
	muladd(n0, SECP256K1_N_C_3);
	extract(m3);
	sumadd(l[4]);
	muladd(n4, SECP256K1_N_C_0);
	muladd(n3, SECP256K1_N_C_1);
	muladd(n2, SECP256K1_N_C_2);
	muladd(n1, SECP256K1_N_C_3);
	sumadd(n0);
	extract(m4);
	sumadd(l[5]);
	muladd(n5, SECP256K1_N_C_0);
	muladd(n4, SECP256K1_N_C_1);
	muladd(n3, SECP256K1_N_C_2);
	muladd(n2, SECP256K1_N_C_3);
	sumadd(n1);
	extract(m5);
	sumadd(l[6]);
	muladd(n6, SECP256K1_N_C_0);
	muladd(n5, SECP256K1_N_C_1);
	muladd(n4, SECP256K1_N_C_2);
	muladd(n3, SECP256K1_N_C_3);
	sumadd(n2);
	extract(m6);
	sumadd(l[7]);
	muladd(n7, SECP256K1_N_C_0);
	muladd(n6, SECP256K1_N_C_1);
	muladd(n5, SECP256K1_N_C_2);
	muladd(n4, SECP256K1_N_C_3);
	sumadd(n3);
	extract(m7);
	muladd(n7, SECP256K1_N_C_1);
	muladd(n6, SECP256K1_N_C_2);
	muladd(n5, SECP256K1_N_C_3);
	sumadd(n4);
	extract(m8);
	muladd(n7, SECP256K1_N_C_2);
	muladd(n6, SECP256K1_N_C_3);
	sumadd(n5);
	extract(m9);
	muladd(n7, SECP256K1_N_C_3);
	sumadd(n6);
	extract(m10);
	sumadd_fast(n7);
	extract_fast(m11);
	VERIFY_CHECK(c0 <= 1);
	m12 = c0;
	c0 = m0; c1 = 0; c2 = 0;
	muladd_fast(m8, SECP256K1_N_C_0);
	extract_fast(p0);
	sumadd_fast(m1);
	muladd(m9, SECP256K1_N_C_0);
	muladd(m8, SECP256K1_N_C_1);
	extract(p1);
	sumadd(m2);
	muladd(m10, SECP256K1_N_C_0);
	muladd(m9, SECP256K1_N_C_1);
	muladd(m8, SECP256K1_N_C_2);
	extract(p2);
	sumadd(m3);
	muladd(m11, SECP256K1_N_C_0);
	muladd(m10, SECP256K1_N_C_1);
	muladd(m9, SECP256K1_N_C_2);
	muladd(m8, SECP256K1_N_C_3);
	extract(p3);
	sumadd(m4);
	muladd(m12, SECP256K1_N_C_0);
	muladd(m11, SECP256K1_N_C_1);
	muladd(m10, SECP256K1_N_C_2);
	muladd(m9, SECP256K1_N_C_3);
	sumadd(m8);
	extract(p4);
	sumadd(m5);
	muladd(m12, SECP256K1_N_C_1);
	muladd(m11, SECP256K1_N_C_2);
	muladd(m10, SECP256K1_N_C_3);
	sumadd(m9);
	extract(p5);
	sumadd(m6);
	muladd(m12, SECP256K1_N_C_2);
	muladd(m11, SECP256K1_N_C_3);
	sumadd(m10);
	extract(p6);
	sumadd_fast(m7);
	muladd_fast(m12, SECP256K1_N_C_3);
	sumadd_fast(m11);
	extract_fast(p7);
	p8 = c0 + m12;
	VERIFY_CHECK(p8 <= 2);

	/* Reduce 258 bits into 256. */
	/* r[0..7] = p[0..7] + p[8] * SECP256K1_N_C. */
	c = p0 + (uint64_t)SECP256K1_N_C_0 * p8;
	r->d[0] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p1 + (uint64_t)SECP256K1_N_C_1 * p8;
	r->d[1] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p2 + (uint64_t)SECP256K1_N_C_2 * p8;
	r->d[2] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p3 + (uint64_t)SECP256K1_N_C_3 * p8;
	r->d[3] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p4 + (uint64_t)p8;
	r->d[4] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p5;
	r->d[5] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p6;
	r->d[6] = c & 0xFFFFFFFFUL; c >>= 32;
	c += p7;
	r->d[7] = c & 0xFFFFFFFFUL; c >>= 32;

	/* Final reduction of r. */
	secp256k1_scalar_reduce(r, c + secp256k1_scalar_check_overflow(r));
}

void secp256k1_scalar_mul_512(uint32_t *l, const secp256k1_scalar *a, const secp256k1_scalar *b) {
	/* 96 bit accumulator. */
	uint32_t c0 = 0, c1 = 0, c2 = 0;

	/* l[0..15] = a[0..7] * b[0..7]. */
	muladd_fast(a->d[0], b->d[0]);
	extract_fast(l[0]);
	muladd(a->d[0], b->d[1]);
	muladd(a->d[1], b->d[0]);
	extract(l[1]);
	muladd(a->d[0], b->d[2]);
	muladd(a->d[1], b->d[1]);
	muladd(a->d[2], b->d[0]);
	extract(l[2]);
	muladd(a->d[0], b->d[3]);
	muladd(a->d[1], b->d[2]);
	muladd(a->d[2], b->d[1]);
	muladd(a->d[3], b->d[0]);
	extract(l[3]);
	muladd(a->d[0], b->d[4]);
	muladd(a->d[1], b->d[3]);
	muladd(a->d[2], b->d[2]);
	muladd(a->d[3], b->d[1]);
	muladd(a->d[4], b->d[0]);
	extract(l[4]);
	muladd(a->d[0], b->d[5]);
	muladd(a->d[1], b->d[4]);
	muladd(a->d[2], b->d[3]);
	muladd(a->d[3], b->d[2]);
	muladd(a->d[4], b->d[1]);
	muladd(a->d[5], b->d[0]);
	extract(l[5]);
	muladd(a->d[0], b->d[6]);
	muladd(a->d[1], b->d[5]);
	muladd(a->d[2], b->d[4]);
	muladd(a->d[3], b->d[3]);
	muladd(a->d[4], b->d[2]);
	muladd(a->d[5], b->d[1]);
	muladd(a->d[6], b->d[0]);
	extract(l[6]);
	muladd(a->d[0], b->d[7]);
	muladd(a->d[1], b->d[6]);
	muladd(a->d[2], b->d[5]);
	muladd(a->d[3], b->d[4]);
	muladd(a->d[4], b->d[3]);
	muladd(a->d[5], b->d[2]);
	muladd(a->d[6], b->d[1]);
	muladd(a->d[7], b->d[0]);
	extract(l[7]);
	muladd(a->d[1], b->d[7]);
	muladd(a->d[2], b->d[6]);
	muladd(a->d[3], b->d[5]);
	muladd(a->d[4], b->d[4]);
	muladd(a->d[5], b->d[3]);
	muladd(a->d[6], b->d[2]);
	muladd(a->d[7], b->d[1]);
	extract(l[8]);
	muladd(a->d[2], b->d[7]);
	muladd(a->d[3], b->d[6]);
	muladd(a->d[4], b->d[5]);
	muladd(a->d[5], b->d[4]);
	muladd(a->d[6], b->d[3]);
	muladd(a->d[7], b->d[2]);
	extract(l[9]);
	muladd(a->d[3], b->d[7]);
	muladd(a->d[4], b->d[6]);
	muladd(a->d[5], b->d[5]);
	muladd(a->d[6], b->d[4]);
	muladd(a->d[7], b->d[3]);
	extract(l[10]);
	muladd(a->d[4], b->d[7]);
	muladd(a->d[5], b->d[6]);
	muladd(a->d[6], b->d[5]);
	muladd(a->d[7], b->d[4]);
	extract(l[11]);
	muladd(a->d[5], b->d[7]);
	muladd(a->d[6], b->d[6]);
	muladd(a->d[7], b->d[5]);
	extract(l[12]);
	muladd(a->d[6], b->d[7]);
	muladd(a->d[7], b->d[6]);
	extract(l[13]);
	muladd_fast(a->d[7], b->d[7]);
	extract_fast(l[14]);
	VERIFY_CHECK(c1 == 0);
	l[15] = c0;
}

void secp256k1_scalar_sqr_512(uint32_t *l, const secp256k1_scalar *a) {
	/* 96 bit accumulator. */
	uint32_t c0 = 0, c1 = 0, c2 = 0;

	/* l[0..15] = a[0..7]^2. */
	muladd_fast(a->d[0], a->d[0]);
	extract_fast(l[0]);
	muladd2(a->d[0], a->d[1]);
	extract(l[1]);
	muladd2(a->d[0], a->d[2]);
	muladd(a->d[1], a->d[1]);
	extract(l[2]);
	muladd2(a->d[0], a->d[3]);
	muladd2(a->d[1], a->d[2]);
	extract(l[3]);
	muladd2(a->d[0], a->d[4]);
	muladd2(a->d[1], a->d[3]);
	muladd(a->d[2], a->d[2]);
	extract(l[4]);
	muladd2(a->d[0], a->d[5]);
	muladd2(a->d[1], a->d[4]);
	muladd2(a->d[2], a->d[3]);
	extract(l[5]);
	muladd2(a->d[0], a->d[6]);
	muladd2(a->d[1], a->d[5]);
	muladd2(a->d[2], a->d[4]);
	muladd(a->d[3], a->d[3]);
	extract(l[6]);
	muladd2(a->d[0], a->d[7]);
	muladd2(a->d[1], a->d[6]);
	muladd2(a->d[2], a->d[5]);
	muladd2(a->d[3], a->d[4]);
	extract(l[7]);
	muladd2(a->d[1], a->d[7]);
	muladd2(a->d[2], a->d[6]);
	muladd2(a->d[3], a->d[5]);
	muladd(a->d[4], a->d[4]);
	extract(l[8]);
	muladd2(a->d[2], a->d[7]);
	muladd2(a->d[3], a->d[6]);
	muladd2(a->d[4], a->d[5]);
	extract(l[9]);
	muladd2(a->d[3], a->d[7]);
	muladd2(a->d[4], a->d[6]);
	muladd(a->d[5], a->d[5]);
	extract(l[10]);
	muladd2(a->d[4], a->d[7]);
	muladd2(a->d[5], a->d[6]);
	extract(l[11]);
	muladd2(a->d[5], a->d[7]);
	muladd(a->d[6], a->d[6]);
	extract(l[12]);
	muladd2(a->d[6], a->d[7]);
	extract(l[13]);
	muladd_fast(a->d[7], a->d[7]);
	extract_fast(l[14]);
	VERIFY_CHECK(c1 == 0);
	l[15] = c0;
}

void secp256k1_scalar_mul(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b) {
	uint32_t l[16];
	secp256k1_scalar_mul_512(l, a, b);
	secp256k1_scalar_reduce_512(r, l);
}

int secp256k1_scalar_shr_int(secp256k1_scalar *r, int n) {
	int ret;
	VERIFY_CHECK(n > 0);
	VERIFY_CHECK(n < 16);
	ret = r->d[0] & ((1 << n) - 1);
	r->d[0] = (r->d[0] >> n) + (r->d[1] << (32 - n));
	r->d[1] = (r->d[1] >> n) + (r->d[2] << (32 - n));
	r->d[2] = (r->d[2] >> n) + (r->d[3] << (32 - n));
	r->d[3] = (r->d[3] >> n) + (r->d[4] << (32 - n));
	r->d[4] = (r->d[4] >> n) + (r->d[5] << (32 - n));
	r->d[5] = (r->d[5] >> n) + (r->d[6] << (32 - n));
	r->d[6] = (r->d[6] >> n) + (r->d[7] << (32 - n));
	r->d[7] = (r->d[7] >> n);
	return ret;
}

void secp256k1_scalar_sqr(secp256k1_scalar *r, const secp256k1_scalar *a) {
	uint32_t l[16];
	secp256k1_scalar_sqr_512(l, a);
	secp256k1_scalar_reduce_512(r, l);
}

int secp256k1_scalar_eq(const secp256k1_scalar *a, const secp256k1_scalar *b) {
	return ((a->d[0] ^ b->d[0]) | (a->d[1] ^ b->d[1]) | (a->d[2] ^ b->d[2]) | (a->d[3] ^ b->d[3]) | (a->d[4] ^ b->d[4]) | (a->d[5] ^ b->d[5]) | (a->d[6] ^ b->d[6]) | (a->d[7] ^ b->d[7])) == 0;
}

void secp256k1_scalar_mul_shift_var(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b, unsigned int shift) {
	uint32_t l[16];
	unsigned int shiftlimbs;
	unsigned int shiftlow;
	unsigned int shifthigh;
	VERIFY_CHECK(shift >= 256);
	secp256k1_scalar_mul_512(l, a, b);
	shiftlimbs = shift >> 5;
	shiftlow = shift & 0x1F;
	shifthigh = 32 - shiftlow;
	r->d[0] = shift < 512 ? (l[0 + shiftlimbs] >> shiftlow | (shift < 480 && shiftlow ? (l[1 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[1] = shift < 480 ? (l[1 + shiftlimbs] >> shiftlow | (shift < 448 && shiftlow ? (l[2 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[2] = shift < 448 ? (l[2 + shiftlimbs] >> shiftlow | (shift < 416 && shiftlow ? (l[3 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[3] = shift < 416 ? (l[3 + shiftlimbs] >> shiftlow | (shift < 384 && shiftlow ? (l[4 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[4] = shift < 384 ? (l[4 + shiftlimbs] >> shiftlow | (shift < 352 && shiftlow ? (l[5 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[5] = shift < 352 ? (l[5 + shiftlimbs] >> shiftlow | (shift < 320 && shiftlow ? (l[6 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[6] = shift < 320 ? (l[6 + shiftlimbs] >> shiftlow | (shift < 288 && shiftlow ? (l[7 + shiftlimbs] << shifthigh) : 0)) : 0;
	r->d[7] = shift < 288 ? (l[7 + shiftlimbs] >> shiftlow) : 0;
	secp256k1_scalar_cadd_bit(r, 0, (l[(shift - 1) >> 5] >> ((shift - 1) & 0x1f)) & 1);
}
