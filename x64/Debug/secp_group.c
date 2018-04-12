#include "secp_static_ctx.c"
#ifndef __constant
#define __constant
#endif

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
#define SECP256K1_GE_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_CONST((i),(j),(k),(l),(m),(n),(o),(p)), 0}
__constant static const secp256k1_ge secp256k1_ge_const_g = SECP256K1_GE_CONST(0x79BE667EUL, 0xF9DCBBACUL, 0x55A06295UL, 0xCE870B07UL, 0x029BFCDBUL, 0x2DCE28D9UL, 0x59F2815BUL, 0x16F81798UL, 0x483ADA77UL, 0x26A3C465UL, 0x5DA4FBFCUL, 0x0E1108A8UL, 0xFD17B448UL, 0xA6855419UL, 0x9C47D08FUL, 0xFB10D4B8UL);

static void secp256k1_gej_clear(secp256k1_gej *r) {
	r->infinity = 0;
	secp256k1_fe_clear(&r->x);
	secp256k1_fe_clear(&r->y);
	secp256k1_fe_clear(&r->z);
}

static void secp256k1_ge_clear(secp256k1_ge *r) {
	r->infinity = 0;
	secp256k1_fe_clear(&r->x);
	secp256k1_fe_clear(&r->y);
}

static void secp256k1_gej_set_ge(secp256k1_gej *r, const secp256k1_ge *a) {
	r->infinity = a->infinity;
	r->x = a->x;
	r->y = a->y;
	secp256k1_fe_set_int(&r->z, 1);
}

static void secp256k1_gej_set_initial(secp256k1_gej *r) {
	const secp256k1_ge a = secp256k1_ge_const_g;
	r->infinity = a.infinity;
	memcpy(&r->x, &a.x, sizeof(a.x));
	memcpy(&r->y, &a.y, sizeof(a.y));
	secp256k1_fe_set_int(&r->z, 1);
}

static void secp256k1_gej_neg(secp256k1_gej *r, const secp256k1_gej *a) {
	r->infinity = a->infinity;
	r->x = a->x;
	r->y = a->y;
	r->z = a->z;
	secp256k1_fe_normalize_weak(&r->y);
	secp256k1_fe_negate(&r->y, &r->y, 1);
}

static void secp256k1_gej_rescale(secp256k1_gej *r, const secp256k1_fe *s) {
	/* Operations: 4 mul, 1 sqr */
	secp256k1_fe zz;
	secp256k1_fe_sqr(&zz, s);
	secp256k1_fe_mul(&r->x, &r->x, &zz);                /* r->x *= s^2 */
	secp256k1_fe_mul(&r->y, &r->y, &zz);
	secp256k1_fe_mul(&r->y, &r->y, s);                  /* r->y *= s^3 */
	secp256k1_fe_mul(&r->z, &r->z, s);                  /* r->z *= s   */
}

static void secp256k1_ge_storage_cmov(secp256k1_ge_storage *r, const secp256k1_ge_storage *a, int flag) {
	secp256k1_fe_storage_cmov(&r->x, &a->x, flag);
	secp256k1_fe_storage_cmov(&r->y, &a->y, flag);
}

static void secp256k1_ge_to_storage(secp256k1_ge_storage *r, const secp256k1_ge *a) {
	secp256k1_fe x, y;
	x = a->x;
	secp256k1_fe_normalize(&x);
	y = a->y;
	secp256k1_fe_normalize(&y);
	secp256k1_fe_to_storage(&r->x, &x);
	secp256k1_fe_to_storage(&r->y, &y);
}

static void secp256k1_ge_from_storage(secp256k1_ge *r, const secp256k1_ge_storage *a) {
	secp256k1_fe_from_storage(&r->x, &a->x);
	secp256k1_fe_from_storage(&r->y, &a->y);
	r->infinity = 0;
}

static void secp256k1_gej_add_ge_var(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b, secp256k1_fe *rzr) {
	/* 8 mul, 3 sqr, 4 normalize, 12 mul_int/add/negate */
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

static void secp256k1_gej_add_ge(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_ge *b) {
	/* Operations: 7 mul, 5 sqr, 4 normalize, 21 mul_int/add/negate/cmov */
	const secp256k1_fe fe_1 = SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 1);
	secp256k1_fe zz, u1, u2, s1, s2, t, tt, m, n, q, rr;
	secp256k1_fe m_alt, rr_alt;
	int infinity, degenerate;

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
														/** If lambda = R/M = 0/0 we have a problem (except in the "trivial"
														*  case that Z = z1z2 = 0, and this is special-cased later on). */
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
