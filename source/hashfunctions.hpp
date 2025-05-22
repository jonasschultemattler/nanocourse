#include "farmhash.h"
// #include <murmurhash.h>

static inline uint64_t stdhash(const uint64_t x) {
    std::hash<uint64_t> hasher;
    return hasher(x);
}

static inline uint64_t wyhash(const uint64_t x) {
    __uint128_t result = x;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

static inline uint64_t farmhash(uint64_t x) {
    return util::Fingerprint(x);
}
