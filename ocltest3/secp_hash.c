#include "secp_static_ctx.c"

typedef struct {
	uint32_t s[32];
	uint32_t buf[16]; /* In big endian */
	size_t bytes;
} secp256k1_sha256_t;

typedef struct {
	secp256k1_sha256_t inner, outer;
} secp256k1_hmac_sha256_t;

typedef struct {
	unsigned char v[32];
	unsigned char k[32];
	int retry;
} secp256k1_rfc6979_hmac_sha256_t;

