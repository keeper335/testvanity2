#pragma once
static void my_secp256k1_ge_set_all_gej_var(secp256k1_ge *r,
	const secp256k1_gej *a);
static void my_secp256k1_gej_add_ge_var(secp256k1_gej *r,
	const secp256k1_gej *a,
	const secp256k1_ge *b);

void scalartohex(unsigned char *buf, secp256k1_scalar *scalar);
void randScalar7Bytes(secp256k1_scalar *scalar, uint8_t setByte);
static int verify_key(const uint8_t result[52]);
int add_prefix2(const char *prefix, uint8_t *pattern);
#ifndef socketpair
#define PORT 19214
int socketpair(_In_ int af, _In_ int type, _In_ int protocol, SOCKET *pair);
#endif