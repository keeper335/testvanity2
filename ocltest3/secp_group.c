#include "secp_static_ctx.c"
#include "secp_field.c"

static void secp256k1_gej_set_ge(secp256k1_gej *r, const secp256k1_ge *a) {
	r->infinity = a->infinity;
	r->x = a->x;
	r->y = a->y;
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
