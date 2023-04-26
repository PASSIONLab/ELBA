#ifndef HASH_FUNCS_H
#define HASH_FUNCS_H

#include <cstdint>
#include <cstddef>

void wanghash64(const void *key, void *hashval);
void wanghash64_inv(const void *hashval, void *key);

void murmurhash3_128(const void *key, uint32_t numbytes, void *out);
void murmurhash3_64(const void *key, uint32_t numbytes, void *out);
void murmurhash3_32(const void *key, uint32_t numbytes, void *out);

uint32_t murmurhash3(const void *key, size_t len, uint32_t seed);

#endif // HASH_FUNCS_H
