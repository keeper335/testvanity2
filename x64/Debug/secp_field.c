#include "secp_static_ctx.c"
//#define SECP256K1_RESTRICT __restrict
#define SECP256K1_RESTRICT

void secp256k1_fe_mul_inner(uint32_t *r, const uint32_t *a, const uint32_t * SECP256K1_RESTRICT b) {
	uint64_t c, d;
	uint64_t u0, u1, u2, u3, u4, u5, u6, u7, u8;
	uint32_t t9, t1, t0, t2, t3, t4, t5, t6, t7;
	uint32_t M = 0x3FFFFFFUL, R0 = 0x3D10UL, R1 = 0x400UL;

	d = (uint64_t)a[0] * b[9]	+ (uint64_t)a[1] * b[8]	+ (uint64_t)a[2] * b[7]	+ (uint64_t)a[3] * b[6]	+ (uint64_t)a[4] * b[5]	
		+ (uint64_t)a[5] * b[4] + (uint64_t)a[6] * b[3]	+ (uint64_t)a[7] * b[2]	+ (uint64_t)a[8] * b[1]	+ (uint64_t)a[9] * b[0];
	t9 = d & M; d >>= 26;
	c = (uint64_t)a[0] * b[0];
	d += (uint64_t)a[1] * b[9]	+ (uint64_t)a[2] * b[8]	+ (uint64_t)a[3] * b[7]	+ (uint64_t)a[4] * b[6]	+ (uint64_t)a[5] * b[5]
		+ (uint64_t)a[6] * b[4]	+ (uint64_t)a[7] * b[3]	+ (uint64_t)a[8] * b[2]	+ (uint64_t)a[9] * b[1];
	u0 = d & M; d >>= 26; c += u0 * R0;
	t0 = c & M; c >>= 26; c += u0 * R1;
	c += (uint64_t)a[0] * b[1]	+ (uint64_t)a[1] * b[0];
	d += (uint64_t)a[2] * b[9]	+ (uint64_t)a[3] * b[8]	+ (uint64_t)a[4] * b[7]	+ (uint64_t)a[5] * b[6]	+ (uint64_t)a[6] * b[5]
		+ (uint64_t)a[7] * b[4]	+ (uint64_t)a[8] * b[3]	+ (uint64_t)a[9] * b[2];
	u1 = d & M; d >>= 26; c += u1 * R0;
	t1 = c & M; c >>= 26; c += u1 * R1;
	c += (uint64_t)a[0] * b[2]	+ (uint64_t)a[1] * b[1]	+ (uint64_t)a[2] * b[0];
	d += (uint64_t)a[3] * b[9]	+ (uint64_t)a[4] * b[8]	+ (uint64_t)a[5] * b[7]	+ (uint64_t)a[6] * b[6]	+ (uint64_t)a[7] * b[5]
		+ (uint64_t)a[8] * b[4]	+ (uint64_t)a[9] * b[3];
	u2 = d & M; d >>= 26; c += u2 * R0;
	t2 = c & M; c >>= 26; c += u2 * R1;
	c += (uint64_t)a[0] * b[3]	+ (uint64_t)a[1] * b[2]	+ (uint64_t)a[2] * b[1]	+ (uint64_t)a[3] * b[0];
	d += (uint64_t)a[4] * b[9]	+ (uint64_t)a[5] * b[8]	+ (uint64_t)a[6] * b[7]	+ (uint64_t)a[7] * b[6]	+ (uint64_t)a[8] * b[5]	+ (uint64_t)a[9] * b[4];
	u3 = d & M; d >>= 26; c += u3 * R0;
	t3 = c & M; c >>= 26; c += u3 * R1;
	c += (uint64_t)a[0] * b[4]	+ (uint64_t)a[1] * b[3]	+ (uint64_t)a[2] * b[2]	+ (uint64_t)a[3] * b[1]	+ (uint64_t)a[4] * b[0];
	d += (uint64_t)a[5] * b[9]	+ (uint64_t)a[6] * b[8]	+ (uint64_t)a[7] * b[7]	+ (uint64_t)a[8] * b[6]	+ (uint64_t)a[9] * b[5];
	u4 = d & M; d >>= 26; c += u4 * R0;
	t4 = c & M; c >>= 26; c += u4 * R1;
	c += (uint64_t)a[0] * b[5]	+ (uint64_t)a[1] * b[4]	+ (uint64_t)a[2] * b[3]	+ (uint64_t)a[3] * b[2]	+ (uint64_t)a[4] * b[1]	+ (uint64_t)a[5] * b[0];
	d += (uint64_t)a[6] * b[9]	+ (uint64_t)a[7] * b[8]	+ (uint64_t)a[8] * b[7]	+ (uint64_t)a[9] * b[6];
	u5 = d & M; d >>= 26; c += u5 * R0;
	t5 = c & M; c >>= 26; c += u5 * R1;
	c += (uint64_t)a[0] * b[6]	+ (uint64_t)a[1] * b[5]	+ (uint64_t)a[2] * b[4]	+ (uint64_t)a[3] * b[3]	+ (uint64_t)a[4] * b[2]
		+ (uint64_t)a[5] * b[1]	+ (uint64_t)a[6] * b[0];
	d += (uint64_t)a[7] * b[9]	+ (uint64_t)a[8] * b[8]	+ (uint64_t)a[9] * b[7];
	u6 = d & M; d >>= 26; c += u6 * R0;
	t6 = c & M; c >>= 26; c += u6 * R1;
	c += (uint64_t)a[0] * b[7]	+ (uint64_t)a[1] * b[6]	+ (uint64_t)a[2] * b[5]	+ (uint64_t)a[3] * b[4]	+ (uint64_t)a[4] * b[3]
		+ (uint64_t)a[5] * b[2]	+ (uint64_t)a[6] * b[1]	+ (uint64_t)a[7] * b[0];
	d += (uint64_t)a[8] * b[9]	+ (uint64_t)a[9] * b[8];
	u7 = d & M; d >>= 26; c += u7 * R0;
	t7 = c & M; c >>= 26; c += u7 * R1;
	c += (uint64_t)a[0] * b[8]	+ (uint64_t)a[1] * b[7]	+ (uint64_t)a[2] * b[6]	+ (uint64_t)a[3] * b[5]	+ (uint64_t)a[4] * b[4]
		+ (uint64_t)a[5] * b[3]	+ (uint64_t)a[6] * b[2]	+ (uint64_t)a[7] * b[1]	+ (uint64_t)a[8] * b[0];
	d += (uint64_t)a[9] * b[9];
	u8 = d & M; d >>= 26; c += u8 * R0;
	r[3] = t3;	r[4] = t4;	r[5] = t5;	r[6] = t6;	r[7] = t7;
	r[8] = c & M; c >>= 26; c += u8 * R1;
	c += d * R0 + t9;
	r[9] = c & (M >> 4); c >>= 22; c += d * (R1 << 4);
	d = c * (R0 >> 4) + t0;
	r[0] = d & M; d >>= 26;
	d += c * (R1 >> 4) + t1;
	r[1] = d & M; d >>= 26;
	d += t2;
	r[2] = (uint32_t)d;
}

void secp256k1_fe_sqr_inner(uint32_t *r, const uint32_t *a) {
	uint64_t c, d;
	uint64_t u0, u1, u2, u3, u4, u5, u6, u7, u8;
	uint32_t t9, t0, t1, t2, t3, t4, t5, t6, t7;
	uint32_t M = 0x3FFFFFFUL, R0 = 0x3D10UL, R1 = 0x400UL;

	d = (uint64_t)(a[0] * 2) * a[9]	+ (uint64_t)(a[1] * 2) * a[8]	+ (uint64_t)(a[2] * 2) * a[7]	+ (uint64_t)(a[3] * 2) * a[6]+ (uint64_t)(a[4] * 2) * a[5];
	t9 = d & M; d >>= 26;
	c = (uint64_t)a[0] * a[0];
	d += (uint64_t)(a[1] * 2) * a[9]+ (uint64_t)(a[2] * 2) * a[8]+ (uint64_t)(a[3] * 2) * a[7]+ (uint64_t)(a[4] * 2) * a[6]+ (uint64_t)a[5] * a[5];
	u0 = d & M; d >>= 26; c += u0 * R0;
	t0 = c & M; c >>= 26; c += u0 * R1;
	c += (uint64_t)(a[0] * 2) * a[1];
	d += (uint64_t)(a[2] * 2) * a[9]+ (uint64_t)(a[3] * 2) * a[8]+ (uint64_t)(a[4] * 2) * a[7]+ (uint64_t)(a[5] * 2) * a[6];
	u1 = d & M; d >>= 26; c += u1 * R0;
	t1 = c & M; c >>= 26; c += u1 * R1;
	c += (uint64_t)(a[0] * 2) * a[2]+ (uint64_t)a[1] * a[1];
	d += (uint64_t)(a[3] * 2) * a[9]+ (uint64_t)(a[4] * 2) * a[8]+ (uint64_t)(a[5] * 2) * a[7]+ (uint64_t)a[6] * a[6];
	u2 = d & M; d >>= 26; c += u2 * R0;
	t2 = c & M; c >>= 26; c += u2 * R1;
	c += (uint64_t)(a[0] * 2) * a[3]+ (uint64_t)(a[1] * 2) * a[2];
	d += (uint64_t)(a[4] * 2) * a[9]+ (uint64_t)(a[5] * 2) * a[8]+ (uint64_t)(a[6] * 2) * a[7];
	u3 = d & M; d >>= 26; c += u3 * R0;
	t3 = c & M; c >>= 26; c += u3 * R1;
	c += (uint64_t)(a[0] * 2) * a[4]+ (uint64_t)(a[1] * 2) * a[3]+ (uint64_t)a[2] * a[2];
	d += (uint64_t)(a[5] * 2) * a[9]+ (uint64_t)(a[6] * 2) * a[8]+ (uint64_t)a[7] * a[7];
	u4 = d & M; d >>= 26; c += u4 * R0;
	t4 = c & M; c >>= 26; c += u4 * R1;
	c += (uint64_t)(a[0] * 2) * a[5]+ (uint64_t)(a[1] * 2) * a[4]+ (uint64_t)(a[2] * 2) * a[3];
	d += (uint64_t)(a[6] * 2) * a[9]+ (uint64_t)(a[7] * 2) * a[8];
	u5 = d & M; d >>= 26; c += u5 * R0;
	t5 = c & M; c >>= 26; c += u5 * R1;
	c += (uint64_t)(a[0] * 2) * a[6]+ (uint64_t)(a[1] * 2) * a[5]+ (uint64_t)(a[2] * 2) * a[4]+ (uint64_t)a[3] * a[3];
	d += (uint64_t)(a[7] * 2) * a[9]+ (uint64_t)a[8] * a[8];
	u6 = d & M; d >>= 26; c += u6 * R0;
	t6 = c & M; c >>= 26; c += u6 * R1;
	c += (uint64_t)(a[0] * 2) * a[7]+ (uint64_t)(a[1] * 2) * a[6]+ (uint64_t)(a[2] * 2) * a[5]+ (uint64_t)(a[3] * 2) * a[4];
	d += (uint64_t)(a[8] * 2) * a[9];
	u7 = d & M; d >>= 26; c += u7 * R0;
	t7 = c & M; c >>= 26; c += u7 * R1;
	c += (uint64_t)(a[0] * 2) * a[8]+ (uint64_t)(a[1] * 2) * a[7]+ (uint64_t)(a[2] * 2) * a[6]+ (uint64_t)(a[3] * 2) * a[5]+ (uint64_t)a[4] * a[4];
	d += (uint64_t)a[9] * a[9];
	u8 = d & M; d >>= 26; c += u8 * R0;
	r[3] = t3;	r[4] = t4;	r[5] = t5;	r[6] = t6;	r[7] = t7;
	r[8] = c & M; c >>= 26; c += u8 * R1;
	c += d * R0 + t9;
	r[9] = c & (M >> 4); c >>= 22; c += d * (R1 << 4);
	d = c * (R0 >> 4) + t0;
	r[0] = d & M; d >>= 26;
	d += c * (R1 >> 4) + t1;
	r[1] = d & M; d >>= 26;
	d += t2;
	r[2] = (uint32_t)d;
}

int secp256k1_fe_is_zero(secp256k1_fe *a) {
	uint32_t *t = a->n;
	return (t[0] | t[1] | t[2] | t[3] | t[4] | t[5] | t[6] | t[7] | t[8] | t[9]) == 0;
}

void secp256k1_fe_clear(secp256k1_fe *a) {
	int i;
	for (i = 0; i<10; i++) {a->n[i] = 0;}
}

void secp256k1_fe_normalize(secp256k1_fe *r) {
	uint32_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4],
		t5 = r->n[5], t6 = r->n[6], t7 = r->n[7], t8 = r->n[8], t9 = r->n[9];
	uint32_t m;
	uint32_t x = t9 >> 22; t9 &= 0x03FFFFFUL;
	t0 += x * 0x3D1UL; t1 += (x << 6);
	t1 += (t0 >> 26); t0 &= 0x3FFFFFFUL;
	t2 += (t1 >> 26); t1 &= 0x3FFFFFFUL;
	t3 += (t2 >> 26); t2 &= 0x3FFFFFFUL; m = t2;
	t4 += (t3 >> 26); t3 &= 0x3FFFFFFUL; m &= t3;
	t5 += (t4 >> 26); t4 &= 0x3FFFFFFUL; m &= t4;
	t6 += (t5 >> 26); t5 &= 0x3FFFFFFUL; m &= t5;
	t7 += (t6 >> 26); t6 &= 0x3FFFFFFUL; m &= t6;
	t8 += (t7 >> 26); t7 &= 0x3FFFFFFUL; m &= t7;
	t9 += (t8 >> 26); t8 &= 0x3FFFFFFUL; m &= t8;
	x = (t9 >> 22) | ((t9 == 0x03FFFFFUL) & (m == 0x3FFFFFFUL)	& ((t1 + 0x40UL + ((t0 + 0x3D1UL) >> 26)) > 0x3FFFFFFUL));
	t0 += x * 0x3D1UL; t1 += (x << 6);
	t1 += (t0 >> 26); t0 &= 0x3FFFFFFUL;
	t2 += (t1 >> 26); t1 &= 0x3FFFFFFUL;
	t3 += (t2 >> 26); t2 &= 0x3FFFFFFUL;
	t4 += (t3 >> 26); t3 &= 0x3FFFFFFUL;
	t5 += (t4 >> 26); t4 &= 0x3FFFFFFUL;
	t6 += (t5 >> 26); t5 &= 0x3FFFFFFUL;
	t7 += (t6 >> 26); t6 &= 0x3FFFFFFUL;
	t8 += (t7 >> 26); t7 &= 0x3FFFFFFUL;
	t9 += (t8 >> 26); t8 &= 0x3FFFFFFUL;
	t9 &= 0x03FFFFFUL;
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
	r->n[5] = t5; r->n[6] = t6; r->n[7] = t7; r->n[8] = t8; r->n[9] = t9;
}

void secp256k1_fe_normalize_weak(secp256k1_fe *r) {
	uint32_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4],
		t5 = r->n[5], t6 = r->n[6], t7 = r->n[7], t8 = r->n[8], t9 = r->n[9];
	uint32_t x = t9 >> 22; t9 &= 0x03FFFFFUL;
	t0 += x * 0x3D1UL; t1 += (x << 6);
	t1 += (t0 >> 26); t0 &= 0x3FFFFFFUL;
	t2 += (t1 >> 26); t1 &= 0x3FFFFFFUL;
	t3 += (t2 >> 26); t2 &= 0x3FFFFFFUL;
	t4 += (t3 >> 26); t3 &= 0x3FFFFFFUL;
	t5 += (t4 >> 26); t4 &= 0x3FFFFFFUL;
	t6 += (t5 >> 26); t5 &= 0x3FFFFFFUL;
	t7 += (t6 >> 26); t6 &= 0x3FFFFFFUL;
	t8 += (t7 >> 26); t7 &= 0x3FFFFFFUL;
	t9 += (t8 >> 26); t8 &= 0x3FFFFFFUL;
	r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
	r->n[5] = t5; r->n[6] = t6; r->n[7] = t7; r->n[8] = t8; r->n[9] = t9;
}

int secp256k1_fe_normalizes_to_zero(secp256k1_fe *r) {
    uint32_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4],
             t5 = r->n[5], t6 = r->n[6], t7 = r->n[7], t8 = r->n[8], t9 = r->n[9];
    uint32_t z0, z1;
    uint32_t x = t9 >> 22; t9 &= 0x03FFFFFUL;
    t0 += x * 0x3D1UL; t1 += (x << 6);
    t1 += (t0 >> 26); t0 &= 0x3FFFFFFUL; z0  = t0; z1  = t0 ^ 0x3D0UL;
    t2 += (t1 >> 26); t1 &= 0x3FFFFFFUL; z0 |= t1; z1 &= t1 ^ 0x40UL;
    t3 += (t2 >> 26); t2 &= 0x3FFFFFFUL; z0 |= t2; z1 &= t2;
    t4 += (t3 >> 26); t3 &= 0x3FFFFFFUL; z0 |= t3; z1 &= t3;
    t5 += (t4 >> 26); t4 &= 0x3FFFFFFUL; z0 |= t4; z1 &= t4;
    t6 += (t5 >> 26); t5 &= 0x3FFFFFFUL; z0 |= t5; z1 &= t5;
    t7 += (t6 >> 26); t6 &= 0x3FFFFFFUL; z0 |= t6; z1 &= t6;
    t8 += (t7 >> 26); t7 &= 0x3FFFFFFUL; z0 |= t7; z1 &= t7;
    t9 += (t8 >> 26); t8 &= 0x3FFFFFFUL; z0 |= t8; z1 &= t8;
                                         z0 |= t9; z1 &= t9 ^ 0x3C00000UL;
    return (z0 == 0) | (z1 == 0x3FFFFFFUL);
}

int secp256k1_fe_normalizes_to_zero_var(secp256k1_fe *r) {
	uint32_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
	uint32_t z0, z1;
	uint32_t x;

	t0 = r->n[0];
	t9 = r->n[9];
	x = t9 >> 22;
	t0 += x * 0x3D1UL;
	z0 = t0 & 0x3FFFFFFUL;
	z1 = z0 ^ 0x3D0UL;
	if ((z0 != 0UL) & (z1 != 0x3FFFFFFUL))	return 0;

	t1 = r->n[1];	t2 = r->n[2];	t3 = r->n[3];	t4 = r->n[4];
	t5 = r->n[5];	t6 = r->n[6];	t7 = r->n[7];	t8 = r->n[8];

	t9 &= 0x03FFFFFUL;
	t1 += (x << 6);

	t1 += (t0 >> 26);
	t2 += (t1 >> 26); t1 &= 0x3FFFFFFUL; z0 |= t1; z1 &= t1 ^ 0x40UL;
	t3 += (t2 >> 26); t2 &= 0x3FFFFFFUL; z0 |= t2; z1 &= t2;
	t4 += (t3 >> 26); t3 &= 0x3FFFFFFUL; z0 |= t3; z1 &= t3;
	t5 += (t4 >> 26); t4 &= 0x3FFFFFFUL; z0 |= t4; z1 &= t4;
	t6 += (t5 >> 26); t5 &= 0x3FFFFFFUL; z0 |= t5; z1 &= t5;
	t7 += (t6 >> 26); t6 &= 0x3FFFFFFUL; z0 |= t6; z1 &= t6;
	t8 += (t7 >> 26); t7 &= 0x3FFFFFFUL; z0 |= t7; z1 &= t7;
	t9 += (t8 >> 26); t8 &= 0x3FFFFFFUL; z0 |= t8; z1 &= t8;
	z0 |= t9; z1 &= t9 ^ 0x3C00000UL;

	return (z0 == 0) | (z1 == 0x3FFFFFFUL);
}

void secp256k1_fe_negate(secp256k1_fe *r, const secp256k1_fe *a, int m) {
	r->n[0] = 0x3FFFC2FUL * 2 * (m + 1) - a->n[0];
	r->n[1] = 0x3FFFFBFUL * 2 * (m + 1) - a->n[1];
	r->n[2] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[2];
	r->n[3] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[3];
	r->n[4] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[4];
	r->n[5] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[5];
	r->n[6] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[6];
	r->n[7] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[7];
	r->n[8] = 0x3FFFFFFUL * 2 * (m + 1) - a->n[8];
	r->n[9] = 0x03FFFFFUL * 2 * (m + 1) - a->n[9];
}

void secp256k1_fe_set_int(secp256k1_fe *r, int a) {
	r->n[0] = a;
	r->n[1] = r->n[2] = r->n[3] = r->n[4] = r->n[5] = r->n[6] = r->n[7] = r->n[8] = r->n[9] = 0;
}

int secp256k1_fe_set_b32(secp256k1_fe *r, const unsigned char *a) {
	int i;
	r->n[0] = r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
	r->n[5] = r->n[6] = r->n[7] = r->n[8] = r->n[9] = 0;
	for (i = 0; i<32; i++) {
		int j;
		for (j = 0; j<4; j++) {
			int limb = (8 * i + 2 * j) / 26;
			int shift = (8 * i + 2 * j) % 26;
			r->n[limb] |= (uint32_t)((a[31 - i] >> (2 * j)) & 0x3) << shift;
		}
	}
	if (r->n[9] == 0x3FFFFFUL && (r->n[8] & r->n[7] & r->n[6] & r->n[5] & r->n[4] & r->n[3] & r->n[2]) == 0x3FFFFFFUL && (r->n[1] + 0x40UL + ((r->n[0] + 0x3D1UL) >> 26)) > 0x3FFFFFFUL) {
		return 0;
	}
	return 1;
}

void secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe *a) {
	int i;
	for (i = 0; i<32; i++) {
		int j;
		int c = 0;
		for (j = 0; j<4; j++) {
			int limb = (8 * i + 2 * j) / 26;
			int shift = (8 * i + 2 * j) % 26;
			c |= ((a->n[limb] >> shift) & 0x3) << (2 * j);
		}
		r[31 - i] = c;
	}
}


void secp256k1_fe_sqr(secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe_sqr_inner(r->n, a->n);
}

void secp256k1_fe_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b) {
	secp256k1_fe_mul_inner(r->n, a->n, b->n);
}

void secp256k1_fe_cmov(secp256k1_fe *r, const secp256k1_fe *a, int flag) {
	uint32_t mask0, mask1;
	mask0 = flag + ~((uint32_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
	r->n[4] = (r->n[4] & mask0) | (a->n[4] & mask1);
	r->n[5] = (r->n[5] & mask0) | (a->n[5] & mask1);
	r->n[6] = (r->n[6] & mask0) | (a->n[6] & mask1);
	r->n[7] = (r->n[7] & mask0) | (a->n[7] & mask1);
	r->n[8] = (r->n[8] & mask0) | (a->n[8] & mask1);
	r->n[9] = (r->n[9] & mask0) | (a->n[9] & mask1);
}

void secp256k1_fe_storage_cmov(secp256k1_fe_storage *r, const secp256k1_fe_storage *a, int flag) {
	uint32_t mask0, mask1;
	mask0 = flag + ~((uint32_t)0);
	mask1 = ~mask0;
	r->n[0] = (r->n[0] & mask0) | (a->n[0] & mask1);
	r->n[1] = (r->n[1] & mask0) | (a->n[1] & mask1);
	r->n[2] = (r->n[2] & mask0) | (a->n[2] & mask1);
	r->n[3] = (r->n[3] & mask0) | (a->n[3] & mask1);
	r->n[4] = (r->n[4] & mask0) | (a->n[4] & mask1);
	r->n[5] = (r->n[5] & mask0) | (a->n[5] & mask1);
	r->n[6] = (r->n[6] & mask0) | (a->n[6] & mask1);
	r->n[7] = (r->n[7] & mask0) | (a->n[7] & mask1);
}

void secp256k1_fe_to_storage(secp256k1_fe_storage *r, const secp256k1_fe *a) {
	r->n[0] = a->n[0] | a->n[1] << 26;
	r->n[1] = a->n[1] >> 6 | a->n[2] << 20;
	r->n[2] = a->n[2] >> 12 | a->n[3] << 14;
	r->n[3] = a->n[3] >> 18 | a->n[4] << 8;
	r->n[4] = a->n[4] >> 24 | a->n[5] << 2 | a->n[6] << 28;
	r->n[5] = a->n[6] >> 4 | a->n[7] << 22;
	r->n[6] = a->n[7] >> 10 | a->n[8] << 16;
	r->n[7] = a->n[8] >> 16 | a->n[9] << 10;
}

void secp256k1_fe_from_storage(secp256k1_fe *r, const secp256k1_fe_storage *a) {
	r->n[0] = a->n[0] & 0x3FFFFFFUL;
	r->n[1] = a->n[0] >> 26 | ((a->n[1] << 6) & 0x3FFFFFFUL);
	r->n[2] = a->n[1] >> 20 | ((a->n[2] << 12) & 0x3FFFFFFUL);
	r->n[3] = a->n[2] >> 14 | ((a->n[3] << 18) & 0x3FFFFFFUL);
	r->n[4] = a->n[3] >> 8 | ((a->n[4] << 24) & 0x3FFFFFFUL);
	r->n[5] = (a->n[4] >> 2) & 0x3FFFFFFUL;
	r->n[6] = a->n[4] >> 28 | ((a->n[5] << 4) & 0x3FFFFFFUL);
	r->n[7] = a->n[5] >> 22 | ((a->n[6] << 10) & 0x3FFFFFFUL);
	r->n[8] = a->n[6] >> 16 | ((a->n[7] << 16) & 0x3FFFFFFUL);
	r->n[9] = a->n[7] >> 10;
}

void secp256k1_fe_add(secp256k1_fe *r, const secp256k1_fe *a) {
	r->n[0] += a->n[0];
	r->n[1] += a->n[1];
	r->n[2] += a->n[2];
	r->n[3] += a->n[3];
	r->n[4] += a->n[4];
	r->n[5] += a->n[5];
	r->n[6] += a->n[6];
	r->n[7] += a->n[7];
	r->n[8] += a->n[8];
	r->n[9] += a->n[9];
}

void secp256k1_fe_mul_int(secp256k1_fe *r, int a) {
	r->n[0] *= a;
	r->n[1] *= a;
	r->n[2] *= a;
	r->n[3] *= a;
	r->n[4] *= a;
	r->n[5] *= a;
	r->n[6] *= a;
	r->n[7] *= a;
	r->n[8] *= a;
	r->n[9] *= a;
}

int secp256k1_fe_is_odd(secp256k1_fe *a) {
	return a->n[0] & 1;
}


void secp256k1_fe_inv(secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe x2, x3, x6, x9, x11, x22, x44, x88, x176, x220, x223, t1;
	int j;

	/** The binary representation of (p - 2) has 5 blocks of 1s, with lengths in
	*  { 1, 2, 22, 223 }. Use an addition chain to calculate 2^n - 1 for each block:
	*  [1], [2], 3, 6, 9, 11, [22], 44, 88, 176, 220, [223]
	*/

	secp256k1_fe_sqr(&x2, a);
	secp256k1_fe_mul(&x2, &x2, a);

	secp256k1_fe_sqr(&x3, &x2);
	secp256k1_fe_mul(&x3, &x3, a);

	x6 = x3;
	for (j = 0; j<3; j++) {
		secp256k1_fe_sqr(&x6, &x6);
	}
	secp256k1_fe_mul(&x6, &x6, &x3);

	x9 = x6;
	for (j = 0; j<3; j++) {
		secp256k1_fe_sqr(&x9, &x9);
	}
	secp256k1_fe_mul(&x9, &x9, &x3);

	x11 = x9;
	for (j = 0; j<2; j++) {
		secp256k1_fe_sqr(&x11, &x11);
	}
	secp256k1_fe_mul(&x11, &x11, &x2);

	x22 = x11;
	for (j = 0; j<11; j++) {
		secp256k1_fe_sqr(&x22, &x22);
	}
	secp256k1_fe_mul(&x22, &x22, &x11);

	x44 = x22;
	for (j = 0; j<22; j++) {
		secp256k1_fe_sqr(&x44, &x44);
	}
	secp256k1_fe_mul(&x44, &x44, &x22);

	x88 = x44;
	for (j = 0; j<44; j++) {
		secp256k1_fe_sqr(&x88, &x88);
	}
	secp256k1_fe_mul(&x88, &x88, &x44);

	x176 = x88;
	for (j = 0; j<88; j++) {
		secp256k1_fe_sqr(&x176, &x176);
	}
	secp256k1_fe_mul(&x176, &x176, &x88);

	x220 = x176;
	for (j = 0; j<44; j++) {
		secp256k1_fe_sqr(&x220, &x220);
	}
	secp256k1_fe_mul(&x220, &x220, &x44);

	x223 = x220;
	for (j = 0; j<3; j++) {
		secp256k1_fe_sqr(&x223, &x223);
	}
	secp256k1_fe_mul(&x223, &x223, &x3);

	/* The final result is then assembled using a sliding window over the blocks. */

	t1 = x223;
	for (j = 0; j<23; j++) {
		secp256k1_fe_sqr(&t1, &t1);
	}
	secp256k1_fe_mul(&t1, &t1, &x22);
	for (j = 0; j<5; j++) {
		secp256k1_fe_sqr(&t1, &t1);
	}
	secp256k1_fe_mul(&t1, &t1, a);
	for (j = 0; j<3; j++) {
		secp256k1_fe_sqr(&t1, &t1);
	}
	secp256k1_fe_mul(&t1, &t1, &x2);
	for (j = 0; j<2; j++) {
		secp256k1_fe_sqr(&t1, &t1);
	}
	secp256k1_fe_mul(r, a, &t1);
}

void secp256k1_fe_inv_var(secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe_inv(r, a);
}

void secp256k1_fe_inv_all_var(size_t len, secp256k1_fe *r, const secp256k1_fe *a) {
	secp256k1_fe u;
	size_t i;
	if (len < 1) return;
	r[0] = a[0];

	i = 0;
	while (++i < len) secp256k1_fe_mul(&r[i], &r[i - 1], &a[i]);

	secp256k1_fe_inv_var(&u, &r[--i]);

	while (i > 0) {
		size_t j = i--;
		secp256k1_fe_mul(&r[j], &r[i], &u);
		secp256k1_fe_mul(&u, &u, &a[j]);
	}

	r[0] = u;
}
