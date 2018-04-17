#define STEP 3072

#include "secp.c"

__kernel void hello_kernel(__global uchar* out) {
    int i = get_global_id(0), j;
	
	out[0] = 'H';
	out[1] = 'e';
	out[2] = 'l';
	out[3] = 'l';
	out[4] = 'O';
	out[5] = '\0';

	secp256k1_ecmult_gen_context sec_ctx;
	secp256k1_scalar scalar_key = { 0, };
	uint8_t *p = (uint8_t *)scalar_key.d;
	secp256k1_gej gej;
	secp256k1_ge ge;
	uint8_t *sha_block = (uint8_t *) out;
	p[0] = 0x34; p[1] = 0x12;
}

__kernel void temp_kernel (__global uchar* out) {

	
	
//secp256k1_ecmult_gen_context_init(&sec_ctx);

//secp256k1_ecmult_gen(&sec_ctx, &gej, &scalar_key);
	//secp256k1_ge_set_gej_var(&ge, &gej);
	//sha_block[0] = (secp256k1_fe_is_odd(&ge.y) ? 0x03 : 0x02);
	//secp256k1_fe_get_b32(sha_block + 1, &ge.x);
}