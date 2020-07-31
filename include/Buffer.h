#ifndef _BUFFER_H_
#define _BUFFER_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#if defined (__cplusplus)
extern "C" {
#endif

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM (128*1024*1024) /* 128 MB */
//#define MAX_ALLTOALL_MEM (12*1024*1024) /* 12 MB */
#endif

typedef enum BUFFER_FLAGS {FIXED_BUFFER = 1, READONLY_BUFFER = 2, EXTERNAL_ALLOC_BUFFER = 4, EOF_BUFFER = 8} BUFFER_FLAGS_TYPE;
typedef struct {
        char *buf;
        size_t pos, len, size, flags;
} BufferBase;
typedef BufferBase *Buffer;

typedef struct _BufferList {
	Buffer buffer;
	struct _BufferList *next;
} _BufferList;
typedef struct {
	_BufferList *head;
	_BufferList *tail;
} BufferList;

Buffer initBuffer(size_t initSize); 
void attachBuffer(Buffer b, char *buf, size_t pos, size_t len, size_t size);
size_t growBuffer(Buffer b, size_t appendSize); 
size_t growBufferMax(Buffer b, size_t requestedSize);
void resetBuffer1(Buffer b);
void resetBuffer(Buffer b);
char * resetRawBuffer(Buffer b, size_t newSize);
void rewindBuffer(Buffer b);

#define freeBuffer(b) _freeBuffer(&(b), __FILENAME__, __LINE__);
void _freeBuffer(Buffer *b, const char *which_file, int which_line);
char * releaseBuffer(Buffer b);
int isValidBuffer(Buffer b);
int isValidBufferNullOK(Buffer b);

BufferList initBufferList();
void freeBufferList(BufferList bl);
Buffer extendBufferList(BufferList bl, size_t initSize);
Buffer getBuffer(BufferList bl);

size_t getPosBuffer(Buffer b);
size_t getLengthBuffer(Buffer b);
size_t getSizeBuffer(Buffer b);
size_t appendNullBuffer(Buffer b);
char * getStartBuffer(Buffer b);
char *getEndBuffer(Buffer b);
char *getCurBuffer(Buffer b);
void truncateBuffer(Buffer b, int pos);

size_t writeBuffer(Buffer b, const void *data, size_t size);
size_t readBuffer(Buffer b, void *data, size_t size);

size_t writeFileBuffer(Buffer b, FILE *f, size_t myOffset);
size_t readFileBuffer(Buffer b, FILE *f, size_t myOffset, size_t len);

size_t scanfBuffer(Buffer b, const char *fmt, ...);

size_t printfBuffer(Buffer b, const char *fmt, ...);
char *fgetsBuffer(Buffer b, size_t size, FILE *stream);
#ifndef NO_GZIP
#include <zlib.h>
char *gzgetsBuffer(Buffer b, size_t size, gzFile gz);
#else
#define gzgetsBuffer fgetsBuffer
#endif
void *memcpyBuffer(Buffer b, const void *src, size_t n);
void *memsetBuffer(Buffer b, int c, size_t n);
char *strcatBuffer(Buffer b, const char *src);
char *strcpyBuffer(Buffer b, const char *src);
char *strncpyBuffer(Buffer b, const char *src, size_t n);
size_t chompBuffer(Buffer b);

void setBufferForFile(Buffer b, FILE *f);

void swapBuffer(Buffer *a, Buffer *b);

#if defined (__cplusplus)
}
#endif

#endif