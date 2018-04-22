//#include "secp.c"
#include "secp_static_ctx.c"
#include "secp_hash.c"
#include "secp_field.c"
#include "secp_group.c"

#define GE_TO_GE(dst, src) do{int i__;for(i__=0;i__<10;i__++) {(dst)->x.n[i__] = (src)->x.n[i__];(dst)->y.n[i__] = (src)->y.n[i__];(dst)->z.n[i__] = (src)->z.n[i__];}(dst)->infinity = (src)->infinity;} while(0)

void secp_fe_get_b32(unsigned char r[32], secp256k1_fe *a) {
	int i,j,c,limb,shift;
	for (i = 0; i<32; i++) {
		c = 0;
		for (j = 0; j<4; j++) {
			limb = (8 * i + 2 * j) / 26;
			shift = (8 * i + 2 * j) % 26;
			c |= ((a->n[limb] >> shift) & 0x3) << (2 * j);
		}
		r[31 - i] = c;
	}
}


__kernel void hello_kernel(__global uchar* out, __global secp256k1_gej *gej_g, __global secp256k1_gej *gej_g_one) {
    int i = get_global_id(0), j, k, c;

	secp256k1_gej gej, gej_one;
	secp256k1_ge ge_offset;//, rslt;
	unsigned char sha_block[33], hash160[21];
	
	secp256k1_gej base[STEP];
	secp256k1_ge rslt[STEP];

	GE_TO_GE(&gej, gej_g);
	GE_TO_GE(&gej_one, gej_g_one);

	secp256k1_ge_set_gej_var(&ge_offset, &gej_one);

	do {
		secp256k1_gej_set_gej(&base[0], &gej);
		//for(k=1;k < STEP;k++) {
			//secp256k1_gej_add_ge(&base[k], &base[k-1], &ge_offset);
		//}
		//my_secp256k1_ge_set_all_gej_var(rslt, base);

		secp256k1_ge_set_gej_var(&rslt[0], &base[0]);

		sha_block[0]=(secp256k1_fe_is_odd(&rslt[0].y) ? 0x03 : 0x02);
		  //secp_fe_get_b32(sha_block+1, &rslt[0].x);
		  for(i=0; i<32;i++) {
			c=0;
			for (j = 0; j<4; j++) {
				int limb = (8 * i + 2 * j) / 26;
				int shift = (8 * i + 2 * j) % 26;
				c |= ((rslt[0].x.n[limb] >> shift) & 0x3) << (2 * j);
			}
			//sha_block[31 - i + 1] = c;
		  }
	//	  HASH160(hash160, sha_block);

		for(k=0;k < STEP;k++) {
		//  sha_block[0]=(secp256k1_fe_is_odd(&rslt[k].y) ? 0x03 : 0x02);
		//  secp256k1_fe_get_b32(sha_block+1, &rslt[k].x);
		//  HASH160(hash160, sha_block);
		  
		  for (j=0; j < 20; j++) out[k*20 + j] = hash160[j];
		}
	} while (0);
	//for(k=0;k < STEP*20;k++) {
	//	out[k] = 'a';
	//}
}


__kernel void hello_kernel2(__global uchar* out) {
    int i = get_global_id(0), j, k;

	for(k=0;k < STEP*20;k++) {
		out[k] = 'a';
	}

}