/*
 * Copyright (c) 2014 GENOME RESEARCH LIMITED
 * Author(s): Robert Davies
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *  . Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimer.
 * 
 *  . Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimer in the documentation 
 * and/or other materials provided with the distribution.
 * 
 *  . Neither the name of the GENOME RESEARCH LIMITED, the WELLCOME TRUST
 * SANGER INSTITUTE nor the names of its contributors may be used to endorse or 
 * promote products derived from this software without specific prior written 
 * permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*
  gcc -Wall -g -I .. -I ../../Misc/ -I../../build_x86_64 `io_lib-config --cflags` -I/usr/include/tcl8.4  ../tg_utils.c test_tg_utils.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "tg_utils.h"

#define TEST_ROUND_TRIP(NAME, ENC_FUNC, DEC_FUNC, TYPE, PRND, PRNX)	\
    int test_ ## NAME(TYPE val, int *test_num) {			\
	unsigned char encoded[16] = { 0 };				\
	char msg[80];							\
	TYPE decoded;							\
	int len1 = ENC_FUNC(val, encoded);				\
	int len2 = DEC_FUNC(encoded, &decoded);				\
	int fail = !(len1 == len2 && decoded == val);			\
	snprintf(msg, sizeof(msg),					\
		 "%-11s val: %" PRND " (0x%" PRNX ")",			\
		 #NAME, val, val);					\
	printf("  Test %2d: %-58s .. %s\n", ++(*test_num), msg,		\
	       !fail ? "pass" : "fail");				\
	if (fail) {							\
	    printf("    => len1 = %d len2 = %d decoded = %" PRNX "\n",	\
		   len1, len2, decoded);				\
	}								\
	return fail ? 1 : 0;						\
    }

TEST_ROUND_TRIP(u7_uint32_t, int2u7, u72int, uint32_t, PRIu32, PRIx32);
TEST_ROUND_TRIP(u7_int32_t,  int2s7, s72int, int32_t,  PRId32, PRIx32);
TEST_ROUND_TRIP(u7_uint64_t, intw2u7, u72intw, uint64_t, PRIu64, PRIx64);
TEST_ROUND_TRIP(u7_int64_t,  intw2s7, s72intw, int64_t,  PRId64, PRIx64);

int test_tg_utils() {
    int test_num = 0;
    int failed = 0;
    printf("Testing tg_utils.c..\n");
    failed += test_u7_uint32_t(0, &test_num);
    failed += test_u7_uint32_t(1, &test_num);
    failed += test_u7_uint32_t(UINT32_MAX - 1, &test_num);
    failed += test_u7_uint32_t(UINT32_MAX, &test_num);
    failed += test_u7_int32_t(INT32_MIN, &test_num);
    failed += test_u7_int32_t(INT32_MIN + 1, &test_num);
    failed += test_u7_int32_t(-1, &test_num);
    failed += test_u7_int32_t(0, &test_num);
    failed += test_u7_int32_t(1, &test_num);
    failed += test_u7_int32_t(INT32_MAX -1, &test_num);
    failed += test_u7_int32_t(INT32_MAX, &test_num);
    failed += test_u7_uint64_t(0, &test_num);
    failed += test_u7_uint64_t(1, &test_num);
    failed += test_u7_uint64_t(UINT64_MAX - 1, &test_num);
    failed += test_u7_uint64_t(UINT64_MAX, &test_num);
    failed += test_u7_int64_t(INT64_MIN, &test_num);
    failed += test_u7_int64_t(INT64_MIN + 1, &test_num);
    failed += test_u7_int64_t(-1, &test_num);
    failed += test_u7_int64_t(0, &test_num);
    failed += test_u7_int64_t(1, &test_num);
    failed += test_u7_int64_t(INT64_MAX -1, &test_num);
    failed += test_u7_int64_t(INT64_MAX, &test_num);
    printf("%d tests, %d passed, %d failed\n",
	   test_num, test_num - failed, failed);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    return test_tg_utils();
}
