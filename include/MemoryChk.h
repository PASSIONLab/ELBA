#ifndef __MEMORY_CHK_H
#define __MEMORY_CHK_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <stdarg.h>

#ifndef __FILENAME__
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#ifndef PAD_ALLOC_BYTES
#ifdef DEBUG
#define PAD_ALLOC_BYTES 64
//#define PAD_ALLOC_BYTES 0
#else
#define PAD_ALLOC_BYTES 0
#endif
#endif


#if defined (__cplusplus)
extern "C" {
#endif

// returns a padded pointer to the middle of the memory block, pads with int64 values of orig
#define pad_memory(_orig, paddedSize) pad_memory_at_line(_orig, paddedSize, PAD_ALLOC_BYTES, __FILENAME__, __LINE__)
void *pad_memory_at_line(void *realPtr, int64_t paddedSize, int PAD_ALLOC, const char *which_file, int which_line);

// returns the real pointer when given a padded pointer
void *pad_memory_real_ptr(void *paddedPtr, int PAD_ALLOC);
// returns te real allocation size (requested+2*PAD_ALLOC) when given a padded pointer
int64_t pad_memory_real_size(void *paddedPtr, int PAD_ALLOC);

// returns the original allocated pointer, and verifies no memory was corrupted.
// Returns NULL on detection of corrup padding and populates msg (if not NULL)
void * pad_memory_check_prefix_at_line(void *realPtr, int PAD_ALLOC, int64_t *origAddr, int64_t *origSize, char *msg, const char * which_file, int which_line);
void * pad_memory_check_suffix_at_line(void *realPtr, int PAD_ALLOC, int64_t *origAddr, int64_t *origSize, char *msg, const char * which_file, int which_line);
#define pad_memory_check(_ptr, msg) pad_memory_check_at_line(_ptr, msg, PAD_ALLOC_BYTES, __FILENAME__, __LINE__)
void *pad_memory_check_at_line(void *paddedPtr, char *msg, int PAD_ALLOC, const char *which_file, int which_line);

#if defined (__cplusplus)
}
#endif

#endif // MEMORY_CHK_H
