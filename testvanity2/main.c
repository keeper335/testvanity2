
#include "externs.h"

const char *adr_to_find[] = { "1LzhS3k3e9Ub8i2W1V8xQFdB8n2MYCHPCa", "17aPYR1m6pVAacXg1PTDDU7XafvK1dxvhi", "15c9mPGLku1HuW9LRtBf4jcHVpBUt8txKz",
"1Dn8NF8qDyyfHMktmuoQLGyjWmZXgvosXf", "1HAX2n9Uruu9YDt4cqRgYcvtGvZj1rbUyt", "1Kn5h2qpgw9mWE5jKpk8PP4qvvJ1QVy8su", "1AVJKwzs9AskraJLGHAZPiaZcrpDr1U6AB" };

const char fname_result[] = "found.txt";

typedef struct {
	uint64_t *thread_count;
	int thread_num;
	SOCKET sock;
	uint8_t **pattern_to_find;
} thread_data_t;

uint64_t *thread_count;
static void manager_loop(int threads, SOCKET sock);
static void announce_result(int found, const uint8_t result[52]);
DWORD WINAPI threadEngine(void *args_);

static void announce_result(int found, const uint8_t result[52])
{
	uint8_t pub_block[RIPEMD160_DIGEST_LENGTH + 5] = { 0, }, checksum[SHA256_DIGEST_LENGTH], wif[35];
	int j;
	char buf[512];
	FILE *fp;

	memset(buf, 0, 256);

	printf("\n");

	/* Display matching keys in hexadecimal */
	sprintf(buf + strlen(buf), "Private match: ");
	for (j = 0;j < 32;j++)
		sprintf(buf + strlen(buf), "%02x", result[j]);
	sprintf(buf + strlen(buf), "\nPublic match:  ");
	for (j = 0;j < 20;j++)
		sprintf(buf + strlen(buf), "%02x", result[j + 32]);

	/* Convert Public Key to Compressed WIF */
	memcpy(pub_block + 1, result + 32, 20);
	/* Compute checksum and copy first 4-bytes to end of public key */
	SHA256(pub_block, RIPEMD160_DIGEST_LENGTH + 1, checksum);
	SHA256(checksum, SHA256_DIGEST_LENGTH, checksum);
	memcpy(pub_block + 21, checksum, 4);
	b58enc(wif, pub_block, sizeof(pub_block));

	sprintf(buf + strlen(buf), "\nAddress:       %s\n---\n", wif);
	printf("%s", buf);

	if ((fp = fopen(fname_result, "a+")) == NULL) {
		printf("Cannot open file %s\n", fname_result);
		exit(1);
	}
	fprintf(fp, "%s", buf);
	fclose(fp);
}

int main(int argc, char *argv[])
{
	int i, status;
	SOCKET sock[2];
	thread_data_t *args[THREADS_NUM];
	HANDLE  hThreadArray[THREADS_NUM];
	uint8_t *pattern_to_find[ADDRESS_NUM];
	
	WSADATA	wsadata;
	WSAStartup(MAKEWORD(2, 2), &wsadata);
	if (socketpair(AF_UNIX, SOCK_DGRAM, 0, (SOCKET *)sock)) {
		printf("socketpair error\n");
		goto error3;
	}

	for (i = 0;i < ADDRESS_NUM; i++) {
		pattern_to_find[i] = (uint8_t*)malloc(sizeof(uint8_t) * 21);
		if (!add_prefix2(adr_to_find[i], pattern_to_find[i])) {
			goto error2;
		}
	}
	
	thread_count = (uint64_t *)malloc(sizeof(uint64_t) * THREADS_NUM);

	for (i = 0;i < THREADS_NUM;i++) {
		args[i] = (thread_data_t *)malloc(sizeof(thread_data_t));
		args[i]->thread_num = i;
		args[i]->sock = sock[1];
		args[i]->thread_count = &thread_count[i];
		args[i]->pattern_to_find = (uint8_t **)&pattern_to_find[0];
		hThreadArray[i] = CreateThread(NULL, 0, threadEngine, args[i], 0, NULL);
		if (hThreadArray[i] == NULL) {
			printf("main error: can't create thread\n");
			goto error1;
		}
	}

	manager_loop(THREADS_NUM, sock[0]);
	printf("Exit\n");

error1:
	WSACleanup();
	for (i = 0; i< THREADS_NUM;i++)
		free(args[i]);
	free(thread_count);
error2:
	for (i = 0; i< ADDRESS_NUM;i++)
		free(pattern_to_find[i]);
error3:
	return 1;
}

static void manager_loop(int threads, SOCKET sock)
{

	fd_set readset;
	struct timeval tv = { 1, 0 };
	char msg[256];
	uint8_t result[52];
	uint64_t prev = 0, last_result = 0, count, avg, count_avg[8];
	int i, ret, len, found = 0, count_index = 0, count_max = 0;

	FD_ZERO(&readset);

	while (1) {
		/* Wait up to 1 second for hashes to be reported */
		FD_SET(sock, &readset);
		if ((ret = select(sock + 1, &readset, NULL, NULL, &tv)) == -1) {
			perror("select");
			return;
		}

		if (ret) {
			/* Read the (PrivKey,PubKey) tuple from the socket */
			//if ((len = read(sock, result, 52)) != 52) {
			if ((len = recv(sock, result, 52, 0)) != 52) {
				/* Datagram read wasn't 52 bytes; ignore message */
				if (len != -1)
					continue;

				/* Something went very wrong if this happens; exit */
				perror("read");
				return;
			}

			/* Verify we received a valid (PrivKey,PubKey) tuple */
			if (!verify_key(result))
				continue;

			announce_result(++found, result);

			/* Reset hash count */
			for (i = 0, count = 0;i < threads;i++)
				count += thread_count[i];
			last_result = count;
			continue;
		}

		/* Reset the select() timer */
		tv.tv_sec = 1, tv.tv_usec = 0;

		/* Collect updated hash counts */
		for (i = 0, count = 0;i < threads;i++)
			count += thread_count[i];
		count_avg[count_index] = count - prev;
		if (++count_index > count_max)
			count_max = count_index;
		if (count_index == NELEM(count_avg))
			count_index = 0;
		prev = count;
		count -= last_result;

		/* Average the last 8 seconds */
		for (i = 0, avg = 0;i < count_max;i++)
			avg += count_avg[i];
		avg /= count_max;

		sprintf(msg, "[%llu Kkey/s][Total %llu]", (avg + 500) / 1000, count);

		/* Display match count */
		if (found) {
			sprintf(msg + strlen(msg), "[Found %d]", found);
		}

		printf("\r%-78.78s", msg);
		fflush(stdout);
	}
}

DWORD WINAPI threadEngine(void *args_)
{
	thread_data_t *args = (thread_data_t *)args_;
	static secp256k1_gej base[STEP];
	static secp256k1_ge rslt[STEP];
	secp256k1_context *sec_ctx;
	secp256k1_scalar scalar_key, scalar_one = { { 1 } }, scalar_step;
	secp256k1_gej temp;
	secp256k1_ge offset;
	int thread = args->thread_num;
	uint8_t **pattern_to_find;
	pattern_to_find = args->pattern_to_find;

	uint8_t sha_block[SHA256_DIGEST_LENGTH + 1], rmd_block[SHA256_DIGEST_LENGTH], result[52], *pubkey = result + 32;
	uint64_t *key = (uint64_t *)result;
	int i, k;//, fd, len;
	int j;

	/* Initialize the secp256k1 context */
	sec_ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN);
	secp256k1_scalar_set_int(&scalar_step, STEP);
	srand(time(NULL) + thread);
rekey:

	randScalar7Bytes(&scalar_key, 0x05 + thread);

	/* Create group elements for both the random private key and the value 1 */
	secp256k1_ecmult_gen(&sec_ctx->ecmult_gen_ctx, &base[STEP - 1], &scalar_key);
	secp256k1_ecmult_gen(&sec_ctx->ecmult_gen_ctx, &temp, &scalar_one);
	secp256k1_ge_set_gej_var(&offset, &temp);

	/* Main Loop */
	printf("\r");  // This magically makes the loop faster by a smidge
	while (1) {
		/* Add 1 in Jacobian coordinates and save the result; repeat STEP times */
		my_secp256k1_gej_add_ge_var(&base[0], &base[STEP - 1], &offset);
		for (k = 1;k < STEP;k++)
			my_secp256k1_gej_add_ge_var(&base[k], &base[k - 1], &offset);

		/* Convert all group elements from Jacobian to affine coordinates */
		my_secp256k1_ge_set_all_gej_var(rslt, base);

		for (k = 0;k < STEP;k++) {
			thread_count[thread]++;

			/* Extract the 33-byte compressed public key from the group element */
			sha_block[0] = (secp256k1_fe_is_odd(&rslt[k].y) ? 0x03 : 0x02);
			secp256k1_fe_get_b32(sha_block + 1, &rslt[k].x);

			/* Hash public key */
			SHA256(sha_block, sizeof(sha_block), rmd_block);
			RIPEMD160(rmd_block, sizeof(rmd_block), pubkey);

			for (i = 0;i < ADDRESS_NUM; i++) {
				if (0 == memcmp(pattern_to_find[i], pubkey, 15)) {
					secp256k1_scalar val1, val2;
					secp256k1_scalar_set_int(&val1, k + 1);
					if (secp256k1_scalar_add(&val2, &scalar_key, &val1))
						printf("\nOverflow \n");
					secp256k1_scalar_get_b32((uint8_t*)key, &val2);
					printf("\nPrivate key found ");
					for (j = 0;j < 32;j++) printf("%02x", result[j]);
					printf(" >>> ");
					for (;j < 52;j++) printf("%02x", result[j]);
					printf("\n");

					if (send(args->sock, result, 52, 0) != 52)
						return 1;
					//goto rekey;
				}
			}
		}

		/* Increment privkey by STEP */
		if (secp256k1_scalar_add(&scalar_key, &scalar_key, &scalar_step)) {
			printf("\nOverflow \n");
			goto rekey;
		}
	}
	return 1;
}
