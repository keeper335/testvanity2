#define STEP 3072

#include "secp.c"

void copy_ge_from_global(secp256k1_gej *dst, __global secp256k1_gej *src) {
	int i;
	for (i = 0; i < 10; i++) {
		dst->x.n[i] = src->x.n[i];
		dst->y.n[i] = src->y.n[i];
		dst->z.n[i] = src->z.n[i];
	}
	dst->infinity = src->infinity;
}


__kernel void hello_kernel(__global uchar* out, __global secp256k1_gej *gej_g, __global unsigned int *out_i) {
    int i = get_global_id(0), j;

	//secp256k1_ecmult_gen_context sec_ctx;
	//secp256k1_scalar scalar_key;
	//uint8_t *p = (uint8_t *)scalar_key.d;
	//memset(p, 0, sizeof(secp256k1_scalar));
	secp256k1_gej gej;
	secp256k1_ge ge;
	uint8_t sha_block[64];
	copy_ge_from_global(&gej, gej_g);

	//uint8_t *sha_block = (uint8_t *) out;
//	p[0] = 0x34; p[1] = 0x12;

	secp256k1_ge_set_gej_var(&ge, &gej);

	out_i[0] = gej.x.n[0];
	out_i[1] = gej.y.n[0];
	out_i[2] = gej.z.n[0];

	sha_block[0] = (secp256k1_fe_is_odd(&ge.y) ? 0x03 : 0x02);
	secp256k1_fe_get_b32(sha_block + 1, &ge.x);

	for (j=0; j < 32; j++)
		out[j] = sha_block[j];

}

void temp_kernel (__global uchar* out) {

	
	
//secp256k1_ecmult_gen_context_init(&sec_ctx);

//secp256k1_ecmult_gen(&sec_ctx, &gej, &scalar_key);
	//secp256k1_ge_set_gej_var(&ge, &gej);
	//sha_block[0] = (secp256k1_fe_is_odd(&ge.y) ? 0x03 : 0x02);
	//secp256k1_fe_get_b32(sha_block + 1, &ge.x);
}