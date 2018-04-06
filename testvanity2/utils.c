//#include "externs.h"
#include "utils.h"
#include "b58.h"

/**** libsecp256k1 Overrides *************************************************/

void my_secp256k1_fe_inv_all_gej_var(secp256k1_fe *r,
	const secp256k1_gej *a)
{
	secp256k1_fe u;
	int i;

	r[0] = a[0].z;

	for (i = 1;i < STEP;i++)
		secp256k1_fe_mul(&r[i], &r[i - 1], &a[i].z);

	secp256k1_fe_inv_var(&u, &r[--i]);

	for (;i > 0;i--) {
		secp256k1_fe_mul(&r[i], &r[i - 1], &u);
		secp256k1_fe_mul(&u, &u, &a[i].z);
	}

	r[0] = u;
}

void my_secp256k1_ge_set_all_gej_var(secp256k1_ge *r,
	const secp256k1_gej *a)
{
	static secp256k1_fe azi[STEP];
	int i;

	my_secp256k1_fe_inv_all_gej_var(azi, a);

	for (i = 0;i < STEP;i++)
		secp256k1_ge_set_gej_zinv(&r[i], &a[i], &azi[i]);
}

void my_secp256k1_gej_add_ge_var(secp256k1_gej *r,
	const secp256k1_gej *a,
	const secp256k1_ge *b)
{
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
	r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3);
	secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
	secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t);
	secp256k1_fe_mul(&r->y, &r->y, &i);
	secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
	secp256k1_fe_add(&r->y, &h3);
}

void scalartohex(unsigned char *buf, secp256k1_scalar *scalar)
{
	uint32_t *p = (uint32_t *)scalar->d;
	sprintf(buf + strlen(buf), "%016lx%016lx%016lx%016lx", p[3], p[2], p[1], p[0]);
}

void randScalar7Bytes(secp256k1_scalar *scalar, uint8_t setByte)
{
	uint8_t *p = (uint8_t *)scalar->d;
	memset(scalar, 0, sizeof(secp256k1_scalar));
	p[4] = rand() & 0xFF;
	p[5] = rand() & 0xFF;
	p[6] = setByte;

	unsigned char buf[128];
	scalartohex(buf, scalar);
	printf("Scalar hex %s\n", buf);
}

int verify_key(const uint8_t result[52])
{
	secp256k1_context *sec_ctx;
	secp256k1_scalar scalar;
	secp256k1_gej gej;
	secp256k1_ge ge;
	uint8_t sha_block[SHA256_DIGEST_LENGTH + 1], rmd_block[SHA256_DIGEST_LENGTH], pubkey[20];
	int ret, overflow;

	/* Initialize the secp256k1 context */
	sec_ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN);

	/* Copy private key to secp256k1 scalar format */
	secp256k1_scalar_set_b32(&scalar, result, &overflow);
	if (overflow) {
		secp256k1_context_destroy(sec_ctx);
		return 0;  /* Invalid private key */
	}

	/* Create a group element for the private key we're verifying */
	secp256k1_ecmult_gen(&sec_ctx->ecmult_gen_ctx, &gej, &scalar);

	/* Convert to affine coordinates */
	secp256k1_ge_set_gej_var(&ge, &gej);

	/* Extract the 33-byte compressed public key from the group element */
	sha_block[0] = (secp256k1_fe_is_odd(&ge.y) ? 0x03 : 0x02);
	secp256k1_fe_get_b32(sha_block + 1, &ge.x);

	/* Hash public key */
	SHA256(sha_block, sizeof(sha_block), rmd_block);
	RIPEMD160(rmd_block, sizeof(rmd_block), pubkey);


	/* Verify that the hashed public key matches the result */
	ret = !memcmp(pubkey, result + 32, 20);

	secp256k1_context_destroy(sec_ctx);
	return ret;
}

int add_prefix2(const char *prefix, uint8_t *pattern)
{
	/* Determine range of matching public keys */
	size_t pattern_sz = 25;
	size_t b58sz = strlen(prefix);
	uint8_t pattern1[32];
	int j;

	if (!b58tobin(pattern1, &pattern_sz, prefix, b58sz)) {
		fprintf(stderr, "Error: Address '%s' contains an invalid character.\n",
			prefix);
		return 0;
	}

	printf("add prefix %s and its pattern ", prefix);
	for (j = 1;j < 21;j++) printf("%02x", pattern1[j]);
	printf("\n");
	memcpy(pattern, pattern1 + 1, 20);

	return 1;
}

int socketpair(int af, int type, int protocol, SOCKET *pair) {
	struct sockaddr_in s0, s1;
	int af_ = AF_INET;
	if ((pair[0] = socket(af_, type, 0)) == INVALID_SOCKET)
	{
		printf("Could not create socket : %d", WSAGetLastError());
	}

	if ((pair[0] = socket(af_, type, IPPROTO_UDP)) == INVALID_SOCKET)
	{
		printf("Could not create socket : %d", WSAGetLastError());
	}

	memset((char *)&s0, 0, sizeof(s0));
	memset((char *)&s1, 0, sizeof(s1));
	s0.sin_family = af_;
	s0.sin_addr.s_addr = INADDR_ANY;
	s0.sin_port = htons(PORT);
	s1.sin_family = af_;
	s1.sin_port = htons(PORT);
	s1.sin_addr.S_un.S_addr = inet_addr("127.0.0.1");
	if (bind(pair[0], (struct sockaddr *)&s0, sizeof(s0)) == SOCKET_ERROR)
	{
		printf("Bind failed with error code : %d", WSAGetLastError());
		return 1;
	}
	return 0;
}
