#pragma once
#include <stdint.h>
#include <stdio.h>

//#include "secp256k1/src/hash.h"

#define SECP256K1_CONTEXT_VERIFY (1 << 0)
#define SECP256K1_CONTEXT_SIGN   (1 << 1)
#define SECP256K1_INLINE __inline

#define ARG_CHECK(cond) do { \
    if (EXPECT(!(cond), 0)) { \
        secp256k1_callback_call(&ctx->illegal_callback, #cond); \
        return 0; \
    } \
} while(0)
#define VERIFY_CHECK(cond) do { (void)(cond); } while(0)
#define VERIFY_SETUP(stmt)


#define SECP256K1_FE_CONST_INNER(d7, d6, d5, d4, d3, d2, d1, d0) {(d0) | (((uint64_t)(d1) & 0xFFFFFUL) << 32), ((uint64_t)(d1) >> 20) | (((uint64_t)(d2)) << 12) | (((uint64_t)(d3) & 0xFFUL) << 44), ((uint64_t)(d3) >> 8) | (((uint64_t)(d4) & 0xFFFFFFFUL) << 24), ((uint64_t)(d4) >> 28) | (((uint64_t)(d5)) << 4) | (((uint64_t)(d6) & 0xFFFFUL) << 36), ((uint64_t)(d6) >> 16) | (((uint64_t)(d7)) << 16)}
#define SECP256K1_FE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {SECP256K1_FE_CONST_INNER((d7), (d6), (d5), (d4), (d3), (d2), (d1), (d0))}
#define SECP256K1_FE_STORAGE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{(d0) | (((uint64_t)(d1)) << 32), (d2) | (((uint64_t)(d3)) << 32), (d4) | (((uint64_t)(d5)) << 32), (d6) | (((uint64_t)(d7)) << 32)}}
#define SECP256K1_GE_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_CONST((i),(j),(k),(l),(m),(n),(o),(p)), 0}
#define SECP256K1_GE_CONST_INFINITY {SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), 1}
#define SECP256K1_GEJ_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_CONST((i),(j),(k),(l),(m),(n),(o),(p)), SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 1), 0}
#define SECP256K1_GEJ_CONST_INFINITY {SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), 1}

#define WINDOW_A 5
#define WINDOW_G 16
#define ECMULT_TABLE_SIZE(w) (1 << ((w)-2))
#define checked_malloc(x,y) malloc(y)


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

typedef struct {
	/* For accelerating the computation of a*P + b*G: */
	secp256k1_ge_storage(*pre_g)[];    /* odd multiples of the generator */
#ifdef USE_ENDOMORPHISM
	secp256k1_ge_storage(*pre_g_128)[]; /* odd multiples of 2^128*generator */
#endif
} secp256k1_ecmult_context;

typedef struct {
	secp256k1_ge_storage(*prec)[64][16]; /* prec[j][i] = 16^j * i * G + U_i */
	secp256k1_scalar blind;
	secp256k1_gej initial;
} secp256k1_ecmult_gen_context;

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

static const secp256k1_ge secp256k1_ge_const_g = SECP256K1_GE_CONST(
	0x79BE667EUL, 0xF9DCBBACUL, 0x55A06295UL, 0xCE870B07UL,
	0x029BFCDBUL, 0x2DCE28D9UL, 0x59F2815BUL, 0x16F81798UL,
	0x483ADA77UL, 0x26A3C465UL, 0x5DA4FBFCUL, 0x0E1108A8UL,
	0xFD17B448UL, 0xA6855419UL, 0x9C47D08FUL, 0xFB10D4B8UL
);

secp256k1_context* secp256k1_context_create(unsigned int flags);
void secp256k1_context_destroy(secp256k1_context* ctx);

/* ecmult h */
void secp256k1_ecmult_context_init(secp256k1_ecmult_context *ctx);
void secp256k1_context_destroy(secp256k1_context* ctx);
void secp256k1_ecmult_context_clear(secp256k1_ecmult_context *ctx);
int secp256k1_ecmult_context_is_built(const secp256k1_ecmult_context *ctx);
void secp256k1_ecmult_context_build(secp256k1_ecmult_context *ctx, const secp256k1_callback *cb);
void secp256k1_ecmult_context_clone(secp256k1_ecmult_context *dst, const secp256k1_ecmult_context *src, const secp256k1_callback *cb);
int secp256k1_ecmult_wnaf(int *wnaf, int len, const secp256k1_scalar *a, int w);
void secp256k1_ecmult(const secp256k1_ecmult_context *ctx, secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_scalar *na, const secp256k1_scalar *ng);


/* ecmult gen h */
void secp256k1_ecmult_gen_context_init(secp256k1_ecmult_gen_context *ctx);
void secp256k1_ecmult_gen(const secp256k1_ecmult_gen_context *ctx, secp256k1_gej *r, const secp256k1_scalar *gn);
void secp256k1_ecmult_gen_blind(secp256k1_ecmult_gen_context *ctx, const unsigned char *seed32);
void secp256k1_ecmult_gen_context_clear(secp256k1_ecmult_gen_context *ctx);
void secp256k1_ecmult_gen_context_build(secp256k1_ecmult_gen_context *ctx, const secp256k1_callback* cb);


/* field h */

void secp256k1_fe_verify(const secp256k1_fe *a);
void secp256k1_fe_normalize(secp256k1_fe *r);
void secp256k1_fe_normalize_weak(secp256k1_fe *r);
void secp256k1_fe_normalize_var(secp256k1_fe *r);
int secp256k1_fe_normalizes_to_zero(secp256k1_fe *r);
int secp256k1_fe_normalizes_to_zero_var(secp256k1_fe *r);
void inline secp256k1_fe_set_int(secp256k1_fe *r, int a);
int inline secp256k1_fe_is_zero(const secp256k1_fe *a);
int inline secp256k1_fe_is_odd(const secp256k1_fe *a);
void inline secp256k1_fe_clear(secp256k1_fe *a);
int secp256k1_fe_cmp_var(const secp256k1_fe *a, const secp256k1_fe *b);
int secp256k1_fe_set_b32(secp256k1_fe *r, const unsigned char *a);
void secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe *a);
void inline secp256k1_fe_negate(secp256k1_fe *r, const secp256k1_fe *a, int m);
void inline secp256k1_fe_mul_int(secp256k1_fe *r, int a);
void inline secp256k1_fe_add(secp256k1_fe *r, const secp256k1_fe *a);
void secp256k1_fe_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b);
void secp256k1_fe_sqr(secp256k1_fe *r, const secp256k1_fe *a);
void inline secp256k1_fe_cmov(secp256k1_fe *r, const secp256k1_fe *a, int flag);
void inline secp256k1_fe_storage_cmov(secp256k1_fe_storage *r, const secp256k1_fe_storage *a, int flag);
void secp256k1_fe_to_storage(secp256k1_fe_storage *r, const secp256k1_fe *a);
void inline secp256k1_fe_from_storage(secp256k1_fe *r, const secp256k1_fe_storage *a);
int secp256k1_fe_equal_var(const secp256k1_fe *a, const secp256k1_fe *b);
int secp256k1_fe_sqrt_var(secp256k1_fe *r, const secp256k1_fe *a);
void secp256k1_fe_inv(secp256k1_fe *r, const secp256k1_fe *a);
void secp256k1_fe_inv_var(secp256k1_fe *r, const secp256k1_fe *a);
void secp256k1_fe_inv_all_var(size_t len, secp256k1_fe *r, const secp256k1_fe *a);
/* group h */

void secp256k1_ge_set_gej_zinv(secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_fe *zi);
void secp256k1_ge_set_xy(secp256k1_ge *r, const secp256k1_fe *x, const secp256k1_fe *y);
int secp256k1_ge_is_infinity(const secp256k1_ge *a);
void secp256k1_ge_neg(secp256k1_ge *r, const secp256k1_ge *a);
void secp256k1_ge_set_gej(secp256k1_ge *r, secp256k1_gej *a);
void secp256k1_ge_set_gej_var(secp256k1_ge *r, secp256k1_gej *a);
void secp256k1_ge_set_all_gej_var(size_t len, secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_callback *cb);
void secp256k1_ge_set_table_gej_var(size_t len, secp256k1_ge *r, const secp256k1_gej *a, const secp256k1_fe *zr);
void secp256k1_ge_globalz_set_table_gej(size_t len, secp256k1_ge *r, secp256k1_fe *globalz, const secp256k1_gej *a, const secp256k1_fe *zr);
void secp256k1_gej_set_infinity(secp256k1_gej *r);
void secp256k1_gej_clear(secp256k1_gej *r);
void secp256k1_ge_clear(secp256k1_ge *r);
int secp256k1_ge_set_xo_var(secp256k1_ge *r, const secp256k1_fe *x, int odd);
void secp256k1_gej_set_ge(secp256k1_gej *r, const secp256k1_ge *a);
int secp256k1_gej_eq_x_var(const secp256k1_fe *x, const secp256k1_gej *a);
void secp256k1_gej_neg(secp256k1_gej *r, const secp256k1_gej *a);
int secp256k1_gej_is_infinity(const secp256k1_gej *a);
int secp256k1_gej_is_valid_var(const secp256k1_gej *a);
int secp256k1_ge_is_valid_var(const secp256k1_ge *a);
void secp256k1_gej_double_var(secp256k1_gej *r, const secp256k1_gej *a, secp256k1_fe *rzr);
void secp256k1_gej_double_nonzero(secp256k1_gej *r, const secp256k1_gej *a, secp256k1_fe *rzr);
void secp256k1_gej_add_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_gej *b, secp256k1_fe *rzr);
void secp256k1_gej_add_ge_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b, secp256k1_fe *rzr);
void secp256k1_gej_add_zinv_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b, const secp256k1_fe *bzinv);
void secp256k1_gej_add_ge(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b);
void secp256k1_gej_rescale(secp256k1_gej *r, const secp256k1_fe *s);
void secp256k1_ge_to_storage(secp256k1_ge_storage *r, const secp256k1_ge *a);
void secp256k1_ge_from_storage(secp256k1_ge *r, const secp256k1_ge_storage *a);
void secp256k1_ge_storage_cmov(secp256k1_ge_storage *r, const secp256k1_ge_storage *a, int flag);

/* scalar h */
#define SECP256K1_N_0 ((uint32_t)0xD0364141UL)
#define SECP256K1_N_1 ((uint32_t)0xBFD25E8CUL)
#define SECP256K1_N_2 ((uint32_t)0xAF48A03BUL)
#define SECP256K1_N_3 ((uint32_t)0xBAAEDCE6UL)
#define SECP256K1_N_4 ((uint32_t)0xFFFFFFFEUL)
#define SECP256K1_N_5 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_6 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_7 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_C_0 (~SECP256K1_N_0 + 1)
#define SECP256K1_N_C_1 (~SECP256K1_N_1)
#define SECP256K1_N_C_2 (~SECP256K1_N_2)
#define SECP256K1_N_C_3 (~SECP256K1_N_3)
#define SECP256K1_N_C_4 (1)
#define SECP256K1_N_H_0 ((uint32_t)0x681B20A0UL)
#define SECP256K1_N_H_1 ((uint32_t)0xDFE92F46UL)
#define SECP256K1_N_H_2 ((uint32_t)0x57A4501DUL)
#define SECP256K1_N_H_3 ((uint32_t)0x5D576E73UL)
#define SECP256K1_N_H_4 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_H_5 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_H_6 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_H_7 ((uint32_t)0x7FFFFFFFUL)
#define muladd(a,b) { \
    uint32_t tl, th; \
    { \
        uint64_t t = (uint64_t)a * b; \
        th = t >> 32;         /* at most 0xFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
    VERIFY_CHECK((c1 >= th) || (c2 != 0)); \
}
#define muladd_fast(a,b) { \
    uint32_t tl, th; \
    { \
        uint64_t t = (uint64_t)a * b; \
        th = t >> 32;         /* at most 0xFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFF */ \
    c1 += th;                 /* never overflows by contract (verified in the next line) */ \
    VERIFY_CHECK(c1 >= th); \
}
#define muladd2(a,b) { \
    uint32_t tl, th, th2, tl2; \
    { \
        uint64_t t = (uint64_t)a * b; \
        th = t >> 32;               /* at most 0xFFFFFFFE */ \
        tl = t; \
    } \
    th2 = th + th;                  /* at most 0xFFFFFFFE (in case th was 0x7FFFFFFF) */ \
    c2 += (th2 < th) ? 1 : 0;       /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((th2 >= th) || (c2 != 0)); \
    tl2 = tl + tl;                  /* at most 0xFFFFFFFE (in case the lowest 63 bits of tl were 0x7FFFFFFF) */ \
    th2 += (tl2 < tl) ? 1 : 0;      /* at most 0xFFFFFFFF */ \
    c0 += tl2;                      /* overflow is handled on the next line */ \
    th2 += (c0 < tl2) ? 1 : 0;      /* second overflow is handled on the next line */ \
    c2 += (c0 < tl2) & (th2 == 0);  /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c0 >= tl2) || (th2 != 0) || (c2 != 0)); \
    c1 += th2;                      /* overflow is handled on the next line */ \
    c2 += (c1 < th2) ? 1 : 0;       /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c1 >= th2) || (c2 != 0)); \
}
#define sumadd(a) { \
    unsigned int over; \
    c0 += (a);                  /* overflow is handled on the next line */ \
    over = (c0 < (a)) ? 1 : 0; \
    c1 += over;                 /* overflow is handled on the next line */ \
    c2 += (c1 < over) ? 1 : 0;  /* never overflows by contract */ \
}
#define sumadd_fast(a) { \
    c0 += (a);                 /* overflow is handled on the next line */ \
    c1 += (c0 < (a)) ? 1 : 0;  /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c1 != 0) | (c0 >= (a))); \
    VERIFY_CHECK(c2 == 0); \
}
#define extract(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = c2; \
    c2 = 0; \
}
#define extract_fast(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = 0; \
    VERIFY_CHECK(c2 == 0); \
}
void secp256k1_scalar_clear(secp256k1_scalar *r);
void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v);
unsigned int secp256k1_scalar_get_bits(const secp256k1_scalar *a, unsigned int offset, unsigned int count);
unsigned int secp256k1_scalar_get_bits_var(const secp256k1_scalar *a, unsigned int offset, unsigned int count);
int secp256k1_scalar_check_overflow(const secp256k1_scalar *a);
int secp256k1_scalar_reduce(secp256k1_scalar *r, uint32_t overflow);
int secp256k1_scalar_add(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b);
void secp256k1_scalar_cadd_bit(secp256k1_scalar *r, unsigned int bit, int flag);
void secp256k1_scalar_set_b32(secp256k1_scalar *r, const unsigned char *b32, int *overflow);
void secp256k1_scalar_get_b32(unsigned char *bin, const secp256k1_scalar* a);
int secp256k1_scalar_is_zero(const secp256k1_scalar *a);
void secp256k1_scalar_negate(secp256k1_scalar *r, const secp256k1_scalar *a);
int secp256k1_scalar_is_one(const secp256k1_scalar *a);
int secp256k1_scalar_is_high(const secp256k1_scalar *a);
int secp256k1_scalar_cond_negate(secp256k1_scalar *r, int flag);
void secp256k1_scalar_reduce_512(secp256k1_scalar *r, const uint32_t *l);
void secp256k1_scalar_mul_512(uint32_t *l, const secp256k1_scalar *a, const secp256k1_scalar *b);
void secp256k1_scalar_sqr_512(uint32_t *l, const secp256k1_scalar *a);
void secp256k1_scalar_mul(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b);
int secp256k1_scalar_shr_int(secp256k1_scalar *r, int n);
void secp256k1_scalar_sqr(secp256k1_scalar *r, const secp256k1_scalar *a);
int secp256k1_scalar_eq(const secp256k1_scalar *a, const secp256k1_scalar *b);
void secp256k1_scalar_mul_shift_var(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b, unsigned int shift);


#define ECMULT_TABLE_GET_GE(r,pre,n,w) do { \
    VERIFY_CHECK(((n) & 1) == 1); \
    VERIFY_CHECK((n) >= -((1 << ((w)-1)) - 1)); \
    VERIFY_CHECK((n) <=  ((1 << ((w)-1)) - 1)); \
    if ((n) > 0) { \
        *(r) = (pre)[((n)-1)/2]; \
    } else { \
        secp256k1_ge_neg((r), &(pre)[(-(n)-1)/2]); \
    } \
} while(0)

#define ECMULT_TABLE_GET_GE_STORAGE(r,pre,n,w) do { \
    VERIFY_CHECK(((n) & 1) == 1); \
    VERIFY_CHECK((n) >= -((1 << ((w)-1)) - 1)); \
    VERIFY_CHECK((n) <=  ((1 << ((w)-1)) - 1)); \
    if ((n) > 0) { \
        secp256k1_ge_from_storage((r), &(pre)[((n)-1)/2]); \
    } else { \
        secp256k1_ge_from_storage((r), &(pre)[(-(n)-1)/2]); \
        secp256k1_ge_neg((r), (r)); \
    } \
} while(0)

