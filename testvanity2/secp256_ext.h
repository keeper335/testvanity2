#pragma once
#include "externs.h"
typedef struct {
	uint32_t d[8];
} secp256k1_scalar;

typedef struct {
	uint64_t n[5];
} secp256k1_fe;

typedef struct {
	uint64_t n[4];
} secp256k1_fe_storage;

typedef struct {
	secp256k1_fe x; /* actual X: x/z^2 */
	secp256k1_fe y; /* actual Y: y/z^3 */
	secp256k1_fe z;
	int infinity; /* whether this represents the point at infinity */
} secp256k1_gej;
typedef struct {
	secp256k1_fe x;
	secp256k1_fe y;
	int infinity; /* whether this represents the point at infinity */
} secp256k1_ge;

typedef struct {
	secp256k1_fe_storage x;
	secp256k1_fe_storage y;
} secp256k1_ge_storage;

typedef void* secp256k1_ecmult_context;

typedef void* secp256k1_ecmult_gen_context;

typedef struct {
	void(*fn)(const char *text, void* data);
	const void* data;
} secp256k1_callback;

struct secp256k1_context_struct {
	secp256k1_ecmult_context ecmult_ctx;
	secp256k1_ecmult_gen_context ecmult_gen_ctx;
	secp256k1_callback illegal_callback;
	secp256k1_callback error_callback;
};
typedef struct secp256k1_context_struct secp256k1_context;

#define SECP256K1_CONTEXT_VERIFY (1 << 0)
#define SECP256K1_CONTEXT_SIGN   (1 << 1)

#ifdef __cplusplus
extern   "C" {
#endif


	void secp256k1_fe_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b);
	void secp256k1_fe_mul_int(secp256k1_fe *r, int a);
	void secp256k1_fe_inv_var(secp256k1_fe *r, const secp256k1_fe *a);
	void secp256k1_fe_sqr(secp256k1_fe *r, const secp256k1_fe *a);
	void secp256k1_fe_normalize_weak(secp256k1_fe *r);
	void secp256k1_fe_negate(secp256k1_fe *r, const secp256k1_fe *a, int m);
	void secp256k1_fe_add(secp256k1_fe *r, const secp256k1_fe *a);

	void secp256k1_ge_set_gej_zinv(secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_fe *zi);


	secp256k1_context* secp256k1_context_create(unsigned int flags);
	void secp256k1_context_destroy(secp256k1_context* ctx);
	void secp256k1_scalar_set_b32(secp256k1_scalar *r, const unsigned char *b32, int *overflow);
	void secp256k1_scalar_get_b32(unsigned char *bin, const secp256k1_scalar* a);
	void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v);
	int secp256k1_scalar_add(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b);

	void secp256k1_ecmult_gen(const secp256k1_ecmult_gen_context *ctx, secp256k1_gej *r, const secp256k1_scalar *gn);
	void secp256k1_ge_set_gej_var(secp256k1_ge *r, secp256k1_gej *a);
	int inline secp256k1_fe_is_odd(const secp256k1_fe *a);
	void secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe *a);

#ifdef __cplusplus
}
#endif

