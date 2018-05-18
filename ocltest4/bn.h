#pragma once
#include <stdio.h>

typedef struct {
	unsigned int words[20];
	unsigned int length;
	unsigned int negative;
} BN_t;

/* predefines */
const BN_t BN_n = { {0x364141, 0x97a334, 0x203bbfd, 0x39abd22, 0x2baaedc, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3fffff}, 10, 0 };
const BN_t BN_nc = { {0x3c9bebf, 0x3685ccb, 0x1fc4402, 0x6542dd, 0x1455123, 0, 0}, 5, 0};
const BN_t BN_nh = { {0x1b20a0, 0x24bd19a, 0x101ddfe, 0x1cd5e91, 0x35d576e, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x1fffff }, 10, 0 };
const BN_t BN_p = { {0x3fffc2f, 0x3ffffbf, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3ffffff, 0x3fffff }, 10, 0 };
const BN_t BN_psn = { {0x3c9baee, 0x3685c8b, 0x1fc4402, 0x6542dd, 0x1455123, 0x0, 0x0, 0x0, 0x0, 0x0 }, 5, 0 };

void BN_strip(BN_t* bn);
void BN_fromNumber(BN_t *bn, int num);
void BN_fromBuffer(BN_t *bn, unsigned char *buf, int buf_len);
void BN_clone(BN_t* to, BN_t* from);
void BN_normSign(BN_t *bn);
bool BN_isEven(BN_t *bn);
bool BN_isOdd(BN_t *bn);
bool BN_isZero(BN_t *bn);
int BN_ucmp(BN_t *bn1, BN_t *bn2);
bool BN_gtOne(BN_t *bn);
bool BN_isOverflow(BN_t *bn);
bool BN_isHigh(BN_t *bn);
bool BN_bitLengthGT256(BN_t *bn);
void BN_iuaddn(BN_t *bn, unsigned int num);
void BN_iadd(BN_t *bn, BN_t *num);
void BN_isub(BN_t *bn, BN_t *num);

void BN_fromNumber(BN_t *bn, int num) {
	bn->length = 1;
	if (num < 0) {
		bn->negative = 1;
		bn->words[0] = -num & 0x03ffffff;
	}
	else {
		bn->negative = 0;
		bn->words[0] = num & 0x03ffffff;
	}
}

void BN_fromBuffer(BN_t *bn, unsigned char *buf, int buf_len) {
	unsigned char b32[32] = { 0, };
	for (int i = 31; i >= 0, buf_len > 0; i--, buf_len--) b32[i] = buf[buf_len];
	bn->words[0] = (b32[28] & 0x03) << 24 | b32[29] << 16 | b32[30] << 8 | b32[31];
	bn->words[1] = (b32[25] & 0x0F) << 22 | b32[26] << 14 | b32[27] << 6 | b32[28] >> 2;
	bn->words[2] = (b32[22] & 0x3F) << 20 | b32[23] << 12 | b32[24] << 4 | b32[25] >> 4;
	bn->words[3] = (b32[19] & 0xFF) << 18 | b32[20] << 10 | b32[21] << 2 | b32[22] >> 6;
	bn->words[4] = (b32[15] & 0x03) << 24 | b32[16] << 16 | b32[17] << 8 | b32[18];
	bn->words[5] = (b32[12] & 0x0F) << 22 | b32[13] << 14 | b32[14] << 6 | b32[15] >> 2;
	bn->words[6] = (b32[9] & 0x3F) << 20 | b32[10] << 12 | b32[11] << 4 | b32[12] >> 4;
	bn->words[7] = (b32[6] & 0xFF) << 18 | b32[7] << 10 | b32[8] << 2 | b32[9] >> 6;
	bn->words[8] = (b32[2] & 0x03) << 24 | b32[3] << 16 | b32[4] << 8 | b32[5];
	bn->words[9] = b32[0] << 14 | b32[1] << 6 | b32[2] >> 2;
	bn->length = 10;
	BN_strip(bn);
}

void BN_clone(BN_t* to, BN_t* from) {
	to->length = from->length;
	to->negative = from->negative;
	for (unsigned int i = 0; i < from->length; i++)
		to->words[i] = from->words[i];
}

void BN_strip(BN_t* bn) {
	int i = 0;
	bn->length = 0;
	for (i = 0; i < 20; i++) {
		if (bn->words[i] == 0)
			break;
		bn->length++;
	}
	for (; i < 20; i++) bn->words[i] = 0;
}

void BN_normSign(BN_t *bn) {
	if (bn->length == 1 && bn->words[0] == 0)
		bn->negative = 0;
}

bool BN_isEven(BN_t *bn) {
	return ((bn->words[0] & 1) == 0);
}

bool BN_isOdd(BN_t *bn) {
	return ((bn->words[0] & 1) == 1);
}

bool BN_isZero(BN_t *bn) {
	return (bn->length == 1 && bn->words[0] == 0);
}

int BN_ucmp(BN_t *bn1, BN_t *bn2) {
	if (bn1->length > bn2->length) return 1;
	else if (bn1->length < bn2->length) return -1;
	else {
		for (int i = bn1->length - 1; i >= 0; --i) {
			if (bn1->words[i] > bn2->words[i]) return 1;
			else if (bn1->words[i] < bn2->words[i]) return -1;
		}
	}
	return 0;
}

bool BN_gtOne(BN_t *bn) {
	return (bn->length > 1 || bn->words[0] > 1);
}

bool BN_isOverflow(BN_t *bn) {
	return 1;//BN_ucmp(bn, BN_n) >= 0;
}

bool BN_isHigh(BN_t *bn) {
	return 1;//BN_ucmp(bn, BN_nh) == 1;
}

bool BN_bitLengthGT256(BN_t *bn) {
	return (bn->length > 10 || (bn->length && bn->words[9] > 0x003fffff));
}

void BN_iuaddn(BN_t *bn, unsigned int num) {
	int i = 0;
	bn->words[0] += num;
	for (; bn->words[i] > 0x03FFFFFF && i < bn->length; ++i) {
		bn->words[i] -= 0x04000000;
		bn->words[i + 1] += 1;
	}
	if (i == bn->length) {
		bn->words[i] = 1;
		bn->length++;
	}
}

void BN_iadd(BN_t *bn, BN_t *num) {
	if (bn->negative != num->negative) {
		if (bn->negative != 0) {
			bn->negative = 0;
			BN_isub(bn, num);
			bn->negative = 1;
		}
		else {
			num->negative = 0;
			BN_isub(bn, num);
			num->negative = 1;
		}
		BN_normSign(bn);
		return;
	}

	BN_t *a, *b;
	unsigned int i, carry, word;
	if (bn->length > num->length) {
		a = bn; b = num;
	}
	else {
		a = num; b = bn;
	}
	for (i = 0, carry = 0; i < b->length; ++i) {
		word = a->words[i] + b->words[i] + carry;
		bn->words[i] = word & 0x03ffffff;
		carry = word >> 26;
	}
	for (; carry != 0 && i < a->length; ++i) {
		word = a->words[i] + carry;
		bn->words[i] = word & 0x03ffffff;
		carry = word >> 26;
	}
	bn->length = a->length;
	if (carry != 0) {
		bn->words[bn->length] = carry;
		bn->length++;
	}
	else if (a != bn) {
		for (; i < a->length; ++i)
			bn->words[i] = a->words[i];
	}
	
}

void BN_isub(BN_t *bn, BN_t *num) {
	if (bn->negative != num->negative) {
		if (bn->negative != 0) {
			bn->negative = 0;
			BN_iadd(bn, num);
			bn->negative = 1;
		}
		else {
			num->negative = 0;
			BN_iadd(bn, num);
			num->negative = 1;
		}
		BN_normSign(bn);
		return;
	}

	int cmp = BN_ucmp(bn, num);
	BN_t *a, *b;
	unsigned int i, carry, word;
	if (cmp == 0) {
		bn->negative = 0;
		bn->length = 1;
		bn->words[0] = 0;
		return;
	}
	else if (bn->length > num->length) {
		a = bn; b = num;
	}
	else {
		a = num; b = bn;
	}
	for (i = 0, carry = 0; i < b->length; ++i) {
		word = a->words[i] - b->words[i] + carry;
		carry = word >> 26;
		bn->words[i] = word & 0x03FFFFFF;
	}
	for (; carry != 0 && i < a->length; ++i) {
		word = a->words[i] + carry;
		carry = word >> 26;
		a->words[i] = word & 0x03FFFFFF;
	}

	if (carry == 0 && i < a->length && a != bn) {
		for (; i < a->length;++i) {
			bn->words[i] = a->words[i];
		}
	}

	if (i > bn->length)	bn->length = i;
	if (a != bn) bn->negative ^= 1;
	BN_strip(bn);
	BN_normSign(bn);
}

unsigned int max_(unsigned int a, unsigned int b) {
	return a > b ? a : b;
}
unsigned int min_(unsigned int a, unsigned int b) {
	return a < b ? a : b;
}

void BN_umulTo(BN_t *num1, BN_t *num2, BN_t *out) {
	out->length = num1->length + num2->length - 1;

	unsigned long long a1 = num1->words[0], b1 = num2->words[0];
	//TODO extend r1;
	unsigned long long r1 = a1 * b1;

	unsigned int carry = (r1 / 0x04000000) | 0, k, maxK, j, maxJ;
	out->words[0] = r1 & 0x03ffffff;
	for (k = 1, maxK = out->length; k < maxK; k++) {
		unsigned int ncarry = carry >> 26;
		unsigned int rword = carry & 0x03ffffff;
		
		for (j = max_(0, k - num1->length + 1), maxJ = min_(k, num2->length - 1); j <= maxJ; j++) {
			unsigned int i__ = k - j;
			unsigned long long a = num1->words[i__];
			unsigned long long b = num2->words[j];
			unsigned long long r = a * b + rword;
			ncarry += (r / 0x04000000) | 0;
			rword = r & 0x03ffffff;
		}
		out->words[k] = rword;
		carry = ncarry;
	}

	if (carry != 0) out->words[out->length++] = carry;
	BN_strip(out);
}

void BN_umulnTo(BN_t *num, unsigned int k, BN_t *out) {
	if (k == 0) {
		out->words[0] = 0;
		out->length = 1;
		return;
	}

	unsigned int i, carry, r;
	unsigned long long r;
	for (i = 0, carry = 0; i < num->length; ++i) {
		r = num->words[i];
		r = r * k + carry;
		out->words[i] = r & 0x03ffffff;
		carry = (r / 0x04000000) | 0;
	}

	if (carry > 0) {
		out->words[i] = carry;
		out->length = num->length + 1;
	}
	else {
		out->length = num->length;
	}
}

void BN_umul(BN_t *bn, BN_t *num, BN_t *out) {
	if (bn->length == 1) {
		BN_umulnTo(num, bn->words[0], out);
	}
	else if (num->length == 1) {
		BN_umulnTo(bn, num->words[0], out);
	}
	else {
		BN_umulTo(bn, num, out);
	}
}