#pragma once
//#include "static_ctx.c"
#include "secp_hash.c"
#include "secp_scalar.c"
#include "secp_group.c"
#include "secp_field.c"

static void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v) {
	r->d[0] = v;
	r->d[1] = 0;
	r->d[2] = 0;
	r->d[3] = 0;
}

static const secp256k1_ge secp256k1_ge_const_g = SECP256K1_GE_CONST(0x79BE667EUL, 0xF9DCBBACUL, 0x55A06295UL, 0xCE870B07UL,	0x029BFCDBUL, 0x2DCE28D9UL, 0x59F2815BUL, 0x16F81798UL,	0x483ADA77UL, 0x26A3C465UL, 0x5DA4FBFCUL, 0x0E1108A8UL,	0xFD17B448UL, 0xA6855419UL, 0x9C47D08FUL, 0xFB10D4B8UL);

static void secp256k1_ecmult_gen_context_init(secp256k1_ecmult_gen_context *ctx) {
	ctx->prec = NULL;
	secp256k1_ecmult_gen_context_build(ctx);
}

static void secp256k1_ecmult_gen_context_build(secp256k1_ecmult_gen_context *ctx) {
	if (ctx->prec != NULL) {
		return;
	}
	ctx->prec = (secp256k1_ge_storage(*)[64][16])secp256k1_ecmult_static_context;
	secp256k1_ecmult_gen_blind(ctx, NULL);
}

static void secp256k1_ecmult_gen_blind(secp256k1_ecmult_gen_context *ctx, const unsigned char *seed32) {
	secp256k1_scalar b;
	secp256k1_gej gb;
	secp256k1_fe s;
	unsigned char nonce32[32];
	secp256k1_rfc6979_hmac_sha256_t rng;
	int retry;
	unsigned char keydata[64] = { 0 };

	/* When seed is NULL, reset the initial point and blinding value. */
	secp256k1_gej_set_ge(&ctx->initial, &secp256k1_ge_const_g);
	secp256k1_gej_neg(&ctx->initial, &ctx->initial);
	secp256k1_scalar_set_int(&ctx->blind, 1);

	/* The prior blinding value (if not reset) is chained forward by including it in the hash. */
	secp256k1_scalar_get_b32(nonce32, &ctx->blind);
	/** Using a CSPRNG allows a failure free interface, avoids needing large amounts of random data,
	*   and guards against weak or adversarial seeds.  This is a simpler and safer interface than
	*   asking the caller for blinding values directly and expecting them to retry on failure.
	*/
	memcpy(keydata, nonce32, 32);

	secp256k1_rfc6979_hmac_sha256_initialize(&rng, keydata, 32);
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