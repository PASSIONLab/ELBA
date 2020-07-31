#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "../include/Buffer.h"
#include "../include/MemoryChk.h"

// returns a new Buffer
Buffer initBuffer(size_t initSize) {
    Buffer b = (BufferBase *) calloc(1, sizeof(BufferBase) + 2 * PAD_ALLOC_BYTES);
    if (PAD_ALLOC_BYTES > 0) b = (BufferBase*)  pad_memory(b, sizeof(BufferBase) + 2 * PAD_ALLOC_BYTES);
    if (b == NULL) { fprintf(stderr, "Could not allocate new Buffer!\n"); fflush(stderr); assert(0); exit(1); }
    assert(b->buf == NULL);
    assert(b->size == 0);
    if (initSize == 0) return b; // 0 initSize just callocs memory
    b->buf = (char*) calloc(initSize + 2 * PAD_ALLOC_BYTES, sizeof(char));
    if (b->buf == NULL) { fprintf(stderr, "Could not allocate %lld bytes into Buffer!\n", (long long int) initSize); fflush(stderr); assert(0); exit(1); }
    if (PAD_ALLOC_BYTES > 0) b->buf = (char*) pad_memory(b->buf, initSize + 2 * PAD_ALLOC_BYTES);
    b->pos = b->len = 0;
    b->size = initSize;
    assert(b->len < b->size);
    b->buf[b->len] = '\0';
    assert(isValidBuffer(b));
    return b;
}

// refer to a raw Buffer
void attachBuffer(Buffer b, char *buf, size_t pos, size_t len, size_t size) {
    assert(b != NULL);
    if (b == NULL) { fprintf(stderr, "Could not attach existing Buffer to null BufferBase pointer!\n"); fflush(stderr); assert(0); exit(1); }
    if (b->buf != NULL || b->size > 0) { fprintf(stderr, "Can not attachBuffer to a non-empty Buffer!  call releaseBuffer first!\n"); fflush(stderr); assert(0); exit(1); }
    b->flags = FIXED_BUFFER | EXTERNAL_ALLOC_BUFFER;
    b->buf = buf;
    b->pos = pos;
    b->len = len;
    b->size = size;
    if (size > len)
        appendNullBuffer(b);
    assert(isValidBuffer(b));
}

// returns the remaining size of the Buffer
size_t growBuffer(Buffer b, size_t appendSize) {
    assert(isValidBufferNullOK(b));
    assert(b->size > 0);
    if (appendSize < 0) { fprintf(stderr, "growBuffer(): Can not grow negative length!\n"); fflush(stderr); assert(0); exit(1); }
    size_t requiredSize = b->len + appendSize + 1;
    if (requiredSize >= b->size) {
        while (requiredSize > b->size) {
            b->size *= 2;
        }
        if ( b->flags & FIXED_BUFFER ) { fprintf(stderr, "growBuffer(): Can not grow a fixed-sized Buffer!\n"); fflush(stderr); assert(0); exit(1); }
        if ( b->flags & EXTERNAL_ALLOC_BUFFER ) { fprintf(stderr, "growBuffer(): Can not grow an externally allocated Buffer!\n"); fflush(stderr); assert(0); exit(1); }
        assert(b->size >= requiredSize);
        if (PAD_ALLOC_BYTES > 0) {
            char msg[256];
            b->buf = pad_memory_check(b->buf, msg);
            if (b->buf == NULL) { fprintf(stderr, "growBuffer(): padding failed! %s\n", msg); fflush(stderr); assert(0); exit(1);  }
        }
        b->buf = (char*) realloc(b->buf, b->size + 2 * PAD_ALLOC_BYTES);
        if (b->buf == NULL)  { fprintf(stderr, "Could not reallocate %lld bytes into Buffer!", (long long int) b->size); fflush(stderr); assert(0); exit(1); }
        if (PAD_ALLOC_BYTES > 0) {
            b->buf = (char*) pad_memory(b->buf, b->size + PAD_ALLOC_BYTES * 2);
        }
        memset(b->buf + b->len, 0, b->size - b->len); // initialize new memory to 0
    }
    assert(b->len + appendSize < b->size);
    assert(b->buf);
    assert(b->size > 0);
    b->buf[b->size-1] = '\0';
    return b->size - b->len;
}

// grows the buffer to the maximum available size up to requestedSize
size_t growBufferMax(Buffer b, size_t requestedSize) {
    assert(isValidBufferNullOK(b));
    assert(b->size > 0);
    if ( b->flags & FIXED_BUFFER ) { fprintf(stderr, "growBufferMax(): Can not grow a fixed-sized Buffer!\n"); fflush(stderr); assert(0); exit(1); }
    if ( b->flags & EXTERNAL_ALLOC_BUFFER ) { fprintf(stderr, "growBufferMax(): Can not grow an externally allocated Buffer!\n"); fflush(stderr); assert(0); exit(1); }
    size_t requiredSize = b->len + requestedSize + 1;
    char *newBuf = NULL;
    if (PAD_ALLOC_BYTES) {
        char msg[256];
        b->buf = pad_memory_check(b->buf, msg);
        if ( b->buf == 0 ) { fprintf(stderr, "growBufferMax(): corrupted memory before realloc! %s\n", msg); fflush(stderr); assert(0); exit(1); }
    }
    while ((newBuf = realloc(b->buf, requiredSize + 2*PAD_ALLOC_BYTES)) == NULL) {
        requiredSize = 3 * requiredSize / 4;
        if (requiredSize < b->size) { fprintf(stderr, "Could not growBufferMax to %lld (beyond %lld)\n", (long long int) requestedSize, (long long int) b->size); fflush(stderr); assert(0); exit(1); }
    }
    b->buf = newBuf;
    b->size = requiredSize;
    if (PAD_ALLOC_BYTES > 0)
        b->buf = pad_memory(b->buf, requiredSize + 2*PAD_ALLOC_BYTES);
    memset(b->buf + b->len, 0, b->size - b->len); // initialize new memory to 0
    return b->size - b->len;
}

// sets the buffer len to 0, does not affect size
void resetBuffer1(Buffer b) {
    assert(b != NULL);
    assert(b->buf != NULL);
    assert(b->size > 0);
    b->pos = b->len = 0;
    b->flags &= ~EOF_BUFFER;
}
void resetBuffer(Buffer b) {
    resetBuffer1(b);
    assert(b->len == 0);
    assert(b->pos == 0);
    b->buf[b->len] = '\0';
}
char * resetRawBuffer(Buffer b, size_t newSize) {
    assert(isValidBuffer(b));
    growBuffer(b, newSize);
    resetBuffer1(b);
    assert(b->size > newSize);
    b->len = newSize;
    b->buf[b->len] = '\0';
    if (b->pos > b->len) b->pos = b->len;
    return b->buf;
}
void rewindBuffer(Buffer b) {
    assert(isValidBuffer(b));
    b->pos = 0;
    b->flags &= ~EOF_BUFFER;
}

// destroys a Buffer
void _freeBuffer(Buffer *_b, const char *which_file, int which_line) {
    assert(_b);
    Buffer b = *_b;
    if (b == NULL) return;
    char msg[256];
    if (b->buf != NULL) {
        if (PAD_ALLOC_BYTES > 0) {
            b->buf = pad_memory_check(b->buf, msg);
            if (b->buf == NULL || strlen(msg) > 0) { fprintf(stderr, "freeBuffer(%p) %s:%d: found corruption in Buffer data pointer! %s\n", b, which_file, which_line, msg); fflush(stderr); assert(0); exit(1); }
        }
        free(b->buf); 
    }
    if ( b->flags & EXTERNAL_ALLOC_BUFFER ) { fprintf(stderr, "freeBuffer(%p) %s:%d: Can not free an externally allocated Buffer!\n", b, which_file, which_line); fflush(stderr); assert(0); exit(1); }
    b->buf = NULL;
    b->pos = b->len = b->size = b->flags = 0;
    if (PAD_ALLOC_BYTES > 0) {
       b = pad_memory_check(b, msg);
       if (b == NULL) { fprintf(stderr, "freeBuffer(%p) %s:%d: found corruption in Buffer pointer! %s\n", *_b, which_file, which_line, msg); fflush(stderr); assert(0); exit(1); }
    }
    free(b);
    *_b = NULL;
}

// releases the buffer memory (free() now expected to be called by user not freeBuffer)
char * releaseBuffer(Buffer b) {
    char *buf = b->buf;
    b->buf = NULL;
    b->pos = b->len = b->size = b->flags = 0;
    return buf;
}

int isValidBuffer(Buffer b) {
     return (b != NULL) & (b->buf != NULL) & (b->size > 0) & (b->pos >= 0) & (b->len >= b->pos) & (b->size >= b->len);
}

int isValidBufferNullOK(Buffer b) {
     return (b != NULL) & (b->buf == NULL ? ((b->size == 0) & (b->pos == 0) & (b->len == 0)) : isValidBuffer(b) );
}


BufferList initBufferList(size_t initSize) {
	BufferList bl;
	bl.head = NULL;
	bl.tail = NULL;
	extendBufferList(bl, initSize);
	return bl;
}

void freeBufferList(BufferList bl) {
	struct _BufferList *ptr = bl.head;
	while (ptr != NULL) {
		struct _BufferList *copy = ptr;
		Buffer b = copy->buffer;
		freeBuffer(b);
		ptr = copy->next;
		copy->buffer = NULL;
		copy->next = NULL;
		char msg[256];
		copy = pad_memory_check(copy, msg);
		if (copy == NULL) { fprintf(stderr, "Error: padded memory is corrupted: %s\n", msg); fflush(stderr); assert(0); exit(1); }
		free(copy);
	}
	bl.head = NULL;
	bl.tail = NULL;
}
	
Buffer extendBufferList(BufferList bl, size_t initSize) {
	struct _BufferList *newnode = (struct _BufferList *) calloc(1, sizeof(struct _BufferList) + 2*PAD_ALLOC_BYTES);
	if (PAD_ALLOC_BYTES > 0) newnode = (struct _BufferList *) pad_memory(newnode, sizeof(struct _BufferList) + 2*PAD_ALLOC_BYTES);
	newnode->buffer = initBuffer(initSize);
	newnode->next = NULL;
	if (bl.tail) {
		bl.tail->next = newnode;
		bl.tail = newnode;
	} else {
		assert(!bl.head);
		bl.head = bl.tail = newnode;
	}
	return getBuffer(bl);
}

Buffer getBuffer(BufferList bl) {
	assert(bl.tail);
	assert(bl.tail->buffer);
	return bl.tail->buffer;
}

size_t getPosBuffer(Buffer b) {
    assert(isValidBufferNullOK(b));
    return b->pos;
}
size_t getLengthBuffer(Buffer b) {
    assert(isValidBufferNullOK(b));
    return b->len;
}
size_t getSizeBuffer(Buffer b) {
    return b->size;
}

size_t appendNullBuffer(Buffer b) {
    assert(isValidBuffer(b));
    if (b->buf[b->len] != '\0') {
        growBuffer(b, 1);
        b->buf[b->len] = '\0';
    }
    return b->len;
}

char * getStartBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf;
}
char *getEndBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf + b->len;
}
char *getCurBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf + b->pos;
}

void truncateBuffer(Buffer b, int pos) {
    assert(isValidBuffer(b));
    if (pos > b->len) {
        growBuffer(b, pos - b->len);
    }
    b->len = pos;
    b->buf[b->len] = '\0';
}

size_t writeBuffer(Buffer b, const void *data, size_t size) {
    assert(isValidBuffer(b));
    memcpyBuffer(b, data, size);
    return size;
}

size_t writeFileBuffer(Buffer b, FILE *f, size_t myOffset) {
    assert(isValidBuffer(b));
    size_t wrote, len = getLengthBuffer(b);
    if ( len == 0 ) return len;
    if ( fseek(f, myOffset, SEEK_SET)  != 0 ) { fprintf(stderr, "ERROR: Could not fseek to %lld! %s\n", (long long int) myOffset, strerror(errno)); return 0; }
    if ( (wrote = fwrite(getStartBuffer(b), 1, len, f)) != len ) { fprintf(stderr, "ERROR: Only wrote %lld of %lld bytes! %s\n", (long long int) wrote, (long long int) len, strerror(errno)); return 0; }
    return len;
}

size_t readFileBuffer(Buffer b, FILE *f, size_t myOffset, size_t len) {
    assert(isValidBuffer(b));
    size_t read;
    if ( fseek(f, myOffset, SEEK_SET)  != 0 ) { fprintf(stderr, "ERROR: Could not fseek to %lld! %s\n", (long long int) myOffset, strerror(errno)); return 0; }
    growBuffer(b, len + 1);
    if ( (read = fread(getEndBuffer(b), 1, len, f)) != len ) { fprintf(stderr, "ERROR: Only read %lld of %lld bytes! %s\n", (long long int) read, (long long int) len, strerror(errno)); return 0; }
    b->len += read;
    b->buf[b->len] = '\0';
    return read;
}

size_t readBuffer(Buffer b, void *data, size_t size) {
    assert(isValidBuffer(b));
    assert(b->pos <= b->len);
    size_t maxReadSize = b->len - b->pos;
    if (size > maxReadSize) {
        if (b->pos >= b->len) {
            if ( (b->flags & EOF_BUFFER) == EOF_BUFFER) {
                fprintf(stderr, "Warning: Attempt to read past buffer pos=%lld len=%lld size=%lld requestSize=%lld\n", (long long int) b->pos, (long long int) b->len, (long long int) b->size, (long long int) size);
            }
            b->flags |= EOF_BUFFER;
        }
        size = maxReadSize;
    }
    if (size) {
        char *x = memcpy(data, b->buf + b->pos, size);
        assert(x == (char*) data);
        b->pos += size;
    }
    return size;
}

size_t scanfBuffer(Buffer b, const char *fmt, ...)
{
    assert(isValidBuffer(b));

    va_list args;
    va_start(args, fmt);
    size_t count = vsscanf(getEndBuffer(b), fmt, args);
    va_end(args);
    return count;
}

// writes to a buffer, reallocating, if necessary
// returns the length of the write
size_t printfBuffer(Buffer b, const char *fmt, ...)
{
    assert(isValidBuffer(b));
    size_t requiredSize = -1;
    assert(b != NULL);
    {
        va_list args;

        va_start(args, fmt);
        requiredSize = vsnprintf(NULL, 0, fmt, args);
        va_end(args);
        if (requiredSize < 0) {
            fprintf(stderr, "Error!  printfBuffer can not process '%s' vsnprintf returned %lld\n", fmt, (long long int) requiredSize); fflush(stderr);
            assert(0);
            return 0;
        }

        assert(requiredSize > 0);
        growBuffer(b, requiredSize+1);
    }
    assert(b->size - b->len > 0);

    va_list args;
    va_start(args, fmt);
    size_t len = vsnprintf(getEndBuffer(b), b->size - b->len, fmt, args);
    va_end(args);

    assert(len == requiredSize);
    b->len += len;
    assert(b->len < b->size);
    return len;
}

// reads one line or end of file, regardless of the line length
char *fgetsBuffer(Buffer b, size_t size, FILE *stream)
{
    assert(isValidBuffer(b));
    size_t oldLen = b->len;
    assert(oldLen < b->size);
    do {
        growBuffer(b, size + 1);
        assert(b->size >= b->len + size + 1);
        char *pos = fgets(getEndBuffer(b), size, stream);
        if (pos) {
            assert(pos == getEndBuffer(b));
            size_t readlen = strlen(pos);
            assert(readlen < size);
            assert(readlen > 0);
            b->len += readlen;
            assert(b->len < b->size);
            assert(pos + readlen == getEndBuffer(b));
            assert(*(getEndBuffer(b)) == '\0');
            if (b->buf[b->len-1] == '\n') break;
        } else {
            // nothing read
            assert(b->len < b->size);
            b->buf[b->len] = '\0';
            break;
        }
    } while(!feof(stream));
    // return the pointer to the beginning of this line
    return (b->len > oldLen) ? (getStartBuffer(b) + oldLen) : NULL;
}

#ifndef NO_GZIP
#include <zlib.h>
char *gzgetsBuffer(Buffer b, size_t size, gzFile gz) {
    assert(isValidBuffer(b));
    size_t oldLen = b->len;
    assert(oldLen < b->size);
    do {
        growBuffer(b, size + 1);
        assert(b->size >= b->len + size + 1);
        char *pos = gzgets(gz, getEndBuffer(b), size);
        if (pos) {
            assert(pos == getEndBuffer(b));
            size_t readlen = strlen(pos);
            assert(readlen < size);
            assert(readlen > 0);
            b->len += readlen;
            assert(b->len < b->size);
            assert(pos + readlen == getEndBuffer(b));
            assert(*(getEndBuffer(b)) == '\0');
            if (b->buf[b->len-1] == '\n') break;
        } else {
            // nothing read
            assert(b->len < b->size);
            b->buf[b->len] = '\0';
            break;
        }
    } while(!gzeof(gz));
    // return the pointer to the beginning of this line
    return (b->len > oldLen) ? (getStartBuffer(b) + oldLen) : NULL;
}
#endif

void *memcpyBuffer(Buffer b, const void *src, size_t len)
{
    assert(isValidBuffer(b));
    growBuffer(b, len+1);
    assert(b->size >= b->len + len + 1);
    void *dst = getEndBuffer(b);
    assert(dst != NULL);
    assert(src != NULL);
    void *ret = memcpy(dst, src, len);
    assert(dst == ret);
    b->len += len;
    assert(b->size > b->len);
    return ret;
}

void *memsetBuffer(Buffer b, int c, size_t n)
{
    assert(isValidBuffer(b));
    growBuffer(b, n+1);
    assert(b->size >= b->len + n + 1);
    void *dst = getEndBuffer(b);
    assert(dst != NULL);
    void *ret = memset(dst, c, n);
    assert(dst == ret);
    b->len += n;
    assert(b->size > b->len);
    return ret;
}

char *strcatBuffer(Buffer b, const char *src)
{
    assert(isValidBuffer(b));
    assert(src != NULL);
    growBuffer(b, strlen(src) + 1);
    char *ret = strcat(getStartBuffer(b), src);
    b->len += strlen(src);
    return ret;
}

char *strcpyBuffer(Buffer b, const char *src)
{
    assert(isValidBuffer(b));
    assert(src != NULL);
    size_t len = strlen(src);
    return strncpyBuffer(b, src, len);
}

char *strncpyBuffer(Buffer b, const char *src, size_t len)
{
    assert(isValidBuffer(b));
    growBuffer(b, len + 1);
    char *ret = strcpy(getEndBuffer(b), src);
    b->len += len;
    assert(b->size > b->len);
    assert(b->buf[b->len] == '\0');
    return ret;
}

void setBufferForFile(Buffer b, FILE *f) {
    memset(b->buf, 0, b->size);
    setvbuf(f, (char*) b->buf, _IOFBF, b->size);
}

// removes a trailing new line, if it exists
// returns the possibly new length
size_t chompBuffer(Buffer b)
{
    assert(isValidBuffer(b));
    assert(b->buf[b->len] == '\0');
    if (b->buf[b->len - 1] == '\n') {
        b->buf[b->len - 1] = '\0';
        b->len--;
    }
    return b->len;
}

void swapBuffer(Buffer *a, Buffer *b) {
    Buffer tmp = *a;
    *a = *b;
    *b = tmp;
}
