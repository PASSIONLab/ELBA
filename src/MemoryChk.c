#ifndef __MEMORY_CHK_H
#define __MEMORY_CHK_H
#include "../include/MemoryChk.h"

#include <inttypes.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

void *pad_memory_check_at_line(void *paddedPtr, char *msg, int PAD_ALLOC, const char *which_file, int which_line);

void *pad_memory_at_line(void *realPtr, int64_t paddedSize, int PAD_ALLOC, const char *which_file, int which_line) {
	if (PAD_ALLOC < sizeof(int64_t)*2) return realPtr;
	assert(paddedSize > 2 * PAD_ALLOC);
	assert(PAD_ALLOC >= sizeof(int64_t)*2);
	assert(PAD_ALLOC % 8 == 0);
	char *orig = (char*) realPtr;
	int64_t addrValue = (int64_t) orig;
	void *userMemStart = orig + PAD_ALLOC;
	void *userMemEnd = orig + (paddedSize - PAD_ALLOC);
	int64_t *beginPad = ((int64_t*) userMemStart);
	while ( --beginPad >= (int64_t*) orig ) {
		*beginPad = addrValue;
		if ( --beginPad >= (int64_t*) orig ) {
			*beginPad = paddedSize;
		}
	}
	assert((char*) beginPad < orig);
	int64_t *endPad = (int64_t*) userMemEnd;
	while (endPad < (int64_t*) (orig + paddedSize)) {
		*(endPad++) = addrValue;
		if (endPad < (int64_t*) (orig + paddedSize)) {
			*(endPad++) = paddedSize;
		}
	}
	assert((char*) endPad == orig + paddedSize);
	assert(userMemStart < userMemEnd);
	assert( pad_memory_check_at_line(userMemStart, NULL, PAD_ALLOC, which_file, which_line) == realPtr );
	return userMemStart;
}

void *pad_memory_real_ptr(void *paddedPtr, int PAD_ALLOC) {
	assert(paddedPtr);
	char *ptr = (char*) paddedPtr;
	if (ptr && PAD_ALLOC) {
		return ptr - PAD_ALLOC;
	} else {
		return ptr;
	}
}

int64_t pad_memory_real_size(void *paddedPtr, int PAD_ALLOC) {
	assert(paddedPtr);
	int64_t *ptr = (int64_t*) paddedPtr;
	if (PAD_ALLOC) {
		ptr--;
		ptr--;
		return *ptr;
	} else {
		return 0;
	}
}

void * pad_memory_check_prefix_at_line(void *realPtr, int PAD_ALLOC, int64_t *origAddr, int64_t *origSize, char *msg, const char * which_file, int which_line) {
	assert(realPtr);
	assert(origAddr);
	assert(origSize);
	char *paddedPtr = ((char*) realPtr) + PAD_ALLOC;
	int64_t *tmp = (int64_t*) paddedPtr;
	if (*origAddr == 0) *origAddr = *(--tmp);
	if (*origSize == 0) *origSize = *(--tmp);
	tmp = (int64_t*) paddedPtr;
	while ( (char*) (--tmp) >= (char*) realPtr ) {
		if (*tmp != *origAddr) {
			if (msg) sprintf(msg + strlen(msg), "%s:%d: pad_memory_check found corruption (%lld bytes into prefix pad)! expected pointer %p got %p\n", which_file, which_line, (long long int) ((char*) tmp - (char*) realPtr), (char*) *origAddr, (char*) *tmp);
			return NULL;
		}
		if ( (char*) (--tmp) >= (char*) realPtr ) {
			if (*tmp != *origSize) {
				if (msg) sprintf(msg + strlen(msg), "%s:%d: pad_memory_check found corruption (%lld bytes into prefix pad)! expected size %lld got %lld\n", which_file, which_line, (long long int) ((char*) tmp - (char*) realPtr), (long long int) *origSize, (long long int) *tmp);
				return NULL;
			}
		}
	}
	return realPtr;
}

void * pad_memory_check_suffix_at_line(void *realPtr, int PAD_ALLOC, int64_t *origAddr, int64_t *origSize, char *msg, const char * which_file, int which_line) {
	assert(realPtr);
	assert(origAddr);
	assert(origSize);
	assert(*origAddr);
	assert(*origSize > 0);
	char *paddedPtr = ((char*) realPtr) + (*origSize - PAD_ALLOC);
	int64_t *tmp = (int64_t*) paddedPtr;
	while ( (char*) tmp < paddedPtr + PAD_ALLOC ) {
		if (*tmp != *origAddr) {
			if (msg) sprintf(msg + strlen(msg), "%s:%d: pad_memory_check found corruption (%lld bytes into suffix pad)! expected pointer %p got %p\n", which_file, which_line, (long long int) ((char*) tmp - (char*) paddedPtr), (char*) *origAddr, (char*) *tmp);
fprintf(stderr, "suffix failed %p vs %p\n", (char*) *origAddr, (char*) *tmp);
			return NULL;
		}
		tmp++;
		if ( (char*) tmp < paddedPtr + PAD_ALLOC ) {
			if (*tmp != *origSize) {
				if (msg) sprintf(msg + strlen(msg), "%s:%d: pad_memory_check found corruption (%lld bytes in to suffix pad)! expected size %lld got %lld\n", which_file, which_line, (long long int) ((char*) tmp - (char*) paddedPtr), (long long int) *origSize, (long long int) *tmp);
				return NULL;
			}
		}
		tmp++;
	}
	return realPtr;
}

void *pad_memory_check_at_line(void *paddedPtr, char *msg, int PAD_ALLOC, const char *which_file, int which_line) {
	if (PAD_ALLOC < sizeof(int64_t)*2) return paddedPtr;
	char *ptr = (char*) paddedPtr;
	int64_t *orig, origAddr, *tmp;
	orig = (int64_t*) pad_memory_real_ptr(paddedPtr, PAD_ALLOC);
	int64_t origSize = pad_memory_real_size(paddedPtr, PAD_ALLOC);
	if (msg) msg[0] = '\0';
	if ( (((char*) orig) + PAD_ALLOC) != (char*) ptr ) {
		if (msg) sprintf(msg + strlen(msg), "%s:%d: pad_memory_check found corruption just before! expected %p got %p \n", which_file, which_line, ptr - PAD_ALLOC, (char*) orig);
		return NULL;
	}
	origAddr = (int64_t) orig;

	int64_t testAddr = origAddr, testSize = origSize;
	void *x = pad_memory_check_prefix_at_line(orig, PAD_ALLOC, &testAddr, &testSize, msg, which_file, which_line);

	if (x == NULL || testAddr != origAddr || testSize != origSize || (msg && strlen(msg))) {
		if (msg) sprintf(msg + strlen(msg), "pad_memory_check_prefix returned=%p testAddr=%lld testSize=%lld \n", x, (long long int) testAddr, (long long int) testSize);
		return NULL;
	}

	x = pad_memory_check_suffix_at_line(orig, PAD_ALLOC, &testAddr, &testSize, msg, which_file, which_line);
	if (x == NULL || testAddr != origAddr || testSize != origSize || (msg && strlen(msg))) {
		if (msg) sprintf(msg + strlen(msg), "pad_memory_check_suffix returned=%p testAddr=%lld testSize=%lld \n", x, (long long int) testAddr, (long long int) testSize);
		return NULL;
	}

	if (x == NULL || (int64_t*) x != orig || (msg && strlen(msg) > 0)) {
		if (msg) sprintf(msg + strlen(msg), "pad_memory_check returned=%p testAddr=%lld testSize=%lld \n", x, (long long int) testAddr, (long long int) testSize);
		return NULL;
	}

	return x;
}

#endif // MEMORY_CHK_H
