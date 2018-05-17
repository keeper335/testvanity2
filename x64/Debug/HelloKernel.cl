//#include "secp.c"
#include "secp_static_ctx.c"
//#include "secp_hash.c"
#include "secp_field.c"
//#include "secp_group.c"

#define GE_TO_GE(dst, src) do{int i__;for(i__=0;i__<10;i__++) {(dst)->x.n[i__] = (src)->x.n[i__];(dst)->y.n[i__] = (src)->y.n[i__];(dst)->z.n[i__] = (src)->z.n[i__];}(dst)->infinity = (src)->infinity;} while(0)

void secp256k1_fe_get_b32_global(__global unsigned char *r, secp256k1_fe *a) {
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

__kernel void hello_kernel(__global unsigned int* out, __global secp256k1_gej *gej_g, __global secp256k1_gej *gej_g_one) {
    int i = get_global_id(0), j, k;

	//secp256k1_gej gej, gej_one;
	//secp256k1_ge ge_offset;//, rslt;
	//uint8_t sha_block[64], hash160[21];
	
	//secp256k1_gej base[STEP];
	//secp256k1_ge rslt[STEP];

	secp256k1_fe z2, z3, z, x, y;

	for (j=0; j<10; j++) {
		x.n[j] = gej_g->x.n[j];
		y.n[j] = gej_g->y.n[j];
		z.n[j] = gej_g->z.n[j];
	}

	//GE_TO_GE(&gej, gej_g);
	//GE_TO_GE(&gej_one, gej_g_one);

	//secp256k1_ge_set_gej_var(&ge_offset, &gej_one);

	do {
		//secp256k1_gej_set_gej(&base[0], &gej);
		//for(k=1;k < STEP;k++) {
//			secp256k1_gej_add_ge(&base[k], &base[k-1], &ge_offset);
	//	}
		//my_secp256k1_ge_set_all_gej_var(rslt, base);

		//for(k=0;k < STEP;k++) {
		  //sha_block[0]=(secp256k1_fe_is_odd(&rslt[k].y) ? 0x03 : 0x02);
		  //secp256k1_fe_get_b32(sha_block+1, &rslt[k].x);
		  //HASH160(hash160, sha_block);
		  
		  /* gej to ge */


		  //for (j=0; j< 10; j++) out[j] = z.n[j];
		  secp256k1_fe_sqr(&z, &z);
		  for (j=0; j< 10; j++) out[j] = z.n[j];
		  //secp256k1_fe_inv_var(&z, &z);
		  
		  //secp256k1_fe_sqr(&z2, &z);
		  //secp256k1_fe_mul(&z3, &z, &z2);

		  //secp256k1_fe_mul(&x, &x, &z2);
		  //for (j=0; j< 10; j++) out[j] = x.n[j];
		  //secp256k1_fe_mul(&y, &y, &z3);
		  //secp256k1_fe_set_int(&z, 1);
		  //x = &gej.x;
		  //y = &gej.y;
		  //out[0]=(secp256k1_fe_is_odd(y) ? 0x03 : 0x02);
		  //secp256k1_fe_get_b32_global(out+1, x);
		  

		  //for (j=0; j < 33; j++) out[j] = sha_block[j];
		//}
	} while (0);

}


__kernel void hello_kernel2(__global uchar* out) {
    int i = get_global_id(0), j, k;

	for(k=0;k < STEP*20;k++) {
		out[k] = 'a';
	}

}