//#include "externs.h"
//int b58tobin(void *bin, size_t *binszp, const char *b58, size_t b58sz);
//int b58enc(char *b58, const void *data, size_t binsz);
#include <windows.h>

#define B58_IN_SIZE 25
#define B58_OUT_SIZE (B58_IN_SIZE + 3)/4
static const char b58digits_map[] = {
	-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
	-1, 0, 1, 2, 3, 4, 5, 6,  7, 8,-1,-1,-1,-1,-1,-1,
	-1, 9,10,11,12,13,14,15, 16,-1,17,18,19,20,21,-1,
	22,23,24,25,26,27,28,29, 30,31,32,-1,-1,-1,-1,-1,
	-1,33,34,35,36,37,38,39, 40,41,42,43,-1,44,45,46,
	47,48,49,50,51,52,53,54, 55,56,57,-1,-1,-1,-1,-1,
};

int b58tobin(void *bin, size_t *binszp, const char *b58, size_t b58sz)
{
	unsigned int binsz = B58_IN_SIZE;
	const unsigned char *b58u = (unsigned char *)b58;
	unsigned char *binu = (unsigned char *)bin;
	unsigned int outi[B58_OUT_SIZE];
	unsigned long t;
	unsigned int c;
	unsigned int i, j;
	unsigned char bytesleft = B58_IN_SIZE % 4;
	unsigned int zeromask = bytesleft ? (0xffffffff << (bytesleft * 8)) : 0;
	unsigned zerocount = 0;

	if (!b58sz)
		b58sz = strlen(b58);

	memset(outi, 0, B58_OUT_SIZE * sizeof(*outi));

	// Leading zeros, just count
	for (i = 0; i < b58sz && b58u[i] == '1'; ++i)
		++zerocount;

	for (; i < b58sz; ++i)
	{
		if (b58u[i] & 0x80)
			// High-bit set on invalid digit
			return 0;
		if (b58digits_map[b58u[i]] == -1)
			// Invalid base58 digit
			return 0;
		c = (unsigned)b58digits_map[b58u[i]];
		for (j = B58_OUT_SIZE; j--; )
		{
			t = ((unsigned long)outi[j]) * 58 + c;
			c = (t & 0x3f00000000) >> 32;
			outi[j] = t & 0xffffffff;
		}
		if (c)
			// Output number too big (carry to the next int32)
			return 0;
		if (outi[0] & zeromask)
			// Output number too big (last int32 filled too far)
			return 0;
	}

	j = 0;
	switch (bytesleft) {
	case 3:
		*(binu++) = (outi[0] & 0xff0000) >> 16;
	case 2:
		*(binu++) = (outi[0] & 0xff00) >> 8;
	case 1:
		*(binu++) = (outi[0] & 0xff);
		++j;
	default:
		break;
	}

	for (; j < B58_OUT_SIZE; ++j)
	{
		*(binu++) = (outi[j] >> 0x18) & 0xff;
		*(binu++) = (outi[j] >> 0x10) & 0xff;
		*(binu++) = (outi[j] >> 8) & 0xff;
		*(binu++) = (outi[j] >> 0) & 0xff;
	}

	// Count canonical base58 byte count
	binu = (unsigned char *)bin;
	for (i = 0; i < B58_IN_SIZE; ++i)
	{
		if (binu[i])
			break;
		--*binszp;
	}
	*binszp += zerocount;

	return 1;
}

static const char b58digits_ordered[] = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

typedef long long ssize_t;
#define B58_ENC_SIZE 26 * 138 / 100 + 1
int b58enc(unsigned char *b58, const void *data, unsigned int binsz)
{
	const unsigned char *bin = (unsigned char *)data;
	unsigned char buf[B58_ENC_SIZE] = { 0, };
	int carry;
	unsigned int i, j, high, zcount = 0, size;

	while (zcount < binsz && !bin[zcount]) ++zcount;

	size = (binsz - zcount) * 138 / 100 + 1;
	for (i = zcount, high = size - 1; i < binsz; ++i, high = j)
	{
		for (carry = bin[i], j = size - 1; (j > high) || carry; --j)
		{
			carry += 256 * buf[j];
			buf[j] = carry % 58;
			carry /= 58;
		}
	}

	for (j = 0; j < size && !buf[j]; ++j);

	if (zcount) memset(b58, '1', zcount);
	for (i = zcount; j < size; ++i, ++j) b58[i] = b58digits_ordered[buf[j]];
	b58[i] = '\0';

	return 1;
}
