#pragma once


typedef struct {
	int negative;
	int words[10];
	int length;
} BN_t;

void BN_fromNumber(BN_t *bn, int num) {
	bn->length = 1;
	if (num < 0) {
		bn->negative = 1;
		bn->words[0] = -num;
	}
	else {
		bn->negative = 0;
		bn->words[0] = num;
	}
}

void BN_fromBuffer(BN_t *bn, char *buf, int buf_len) {
	int len = 0, i = buf_len-1;
	
	bn->words[0] = (int)buf[3] << 24 | (int)buf[2] << 16 | (int)buf[1] << 8 | (int)buf[0];

	while (i >= 0) {
		bn->words[len] = buf[]
	}

}