#include "mp.hpp"
#include "fq.hpp"

#include <chrono>
#include <cstdio>
#include "gtest/gtest.h"

static inline void set(uint64_t r[MP_N64], uint64_t l0, uint64_t l1, uint64_t l2, uint64_t l3) {
    r[0] = l0; r[1] = l1; r[2] = l2; r[3] = l3;
}

// 64x64 -> 128 without __uint128_t via 32-bit decomposition
static inline void mul64_no128(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
    const uint64_t a0 = (uint32_t)a;
    const uint64_t a1 = a >> MP_N;
    const uint64_t b0 = (uint32_t)b;
    const uint64_t b1 = b >> MP_N;

    const uint64_t p00 = a0 * b0;
    const uint64_t p01 = a0 * b1;
    const uint64_t p10 = a1 * b0;
    const uint64_t p11 = a1 * b1;

    uint64_t mid = (p00 >> MP_N) + (uint32_t)p01 + (uint32_t)p10;
    lo = (p00 & 0xFFFFFFFFULL) | (mid << MP_N);
    hi = p11 + (p01 >> MP_N) + (p10 >> MP_N) + (mid >> MP_N);
}

static inline uint64_t ref_addmul_no128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t k) {
    uint64_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        uint64_t hi, lo;
        mul64_no128(a[i], k, lo, hi);

        // t = r[i] + lo + carry
        uint64_t t = r[i] + lo;
        uint64_t c1 = (t < r[i]) ? 1 : 0;
        t += carry;
        uint64_t c2 = (t < carry) ? 1 : 0;

        r[i] = t;

        carry = hi + c1 + c2;
    }
    return carry;
}

#if defined(__SIZEOF_INT128__)
static inline uint64_t ref_addmul_128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t k) {
    uint64_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        unsigned __int128 prod = (unsigned __int128)a[i] * (unsigned __int128)k;
        uint64_t lo = (uint64_t)prod;
        uint64_t hi = (uint64_t)(prod >> 2*MP_N);

        unsigned __int128 sum = (unsigned __int128)r[i] + lo + carry;
        r[i] = (uint64_t)sum;

        carry = hi + (uint64_t)(sum >> 2*MP_N);
    }
    return carry;
}
#endif

static inline uint64_t ref_mul_128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        __uint128_t p = (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)p;
        carry = p >> 2*MP_N;
    }
    return (uint64_t)carry;
#else
    (void)r; (void)a; (void)b;
    return 0;
#endif
}

static inline uint64_t ref_mul_no128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t b) {
    uint64_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        uint64_t plo, phi;
        mul64_no128(a[i], b, plo, phi);

        // add carry to lo
        uint64_t lo2 = plo + carry;
        uint64_t c1 = (lo2 < plo) ? 1u : 0u;

        r[i] = lo2;
        carry = phi + c1;
    }
    return carry;
}

static inline uint64_t now_ns() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
}

TEST(mp_perf_1m, addmul_and_mul_timing) {
    constexpr size_t N = 1000000;

    mp_uint_t a{}, r{};
    volatile uint64_t sink = 0;

    auto bench_addmul = [&](auto &&fn, const char *name) {
        uint64_t t0 = now_ns();
        for (size_t i = 0; i < N; i++) {
            set(a,
                      0x0123456789ABCDEFULL + i,
                      0x0FEDCBA987654321ULL ^ (i * 23),
                      0xAAAAAAAAAAAAAAAAULL + (i * 29),
                      0x5555555555555555ULL ^ (i * 31));

            set(r,
                      0x9E3779B97F4A7C15ULL + i,
                      0xBF58476D1CE4E5B9ULL ^ (i << 1),
                      0x94D049BB133111EBULL + (i * 3),
                      0xD6E8FEB86659FD93ULL ^ (i * 7));

            uint64_t k = 0x9E3779B185EBCA87ULL + (i * 37);

            uint64_t c = fn(r, a, k);
            sink ^= r[0] + c;
        }
        uint64_t t1 = now_ns();
        std::printf("[perf] %s: %.3f ms (N=%zu)\n", name, (double)(t1 - t0) / 1e6, N);
    };

    auto bench_mul = [&](auto &&fn, const char *name) {
        uint64_t t0 = now_ns();
        for (size_t i = 0; i < N; i++) {
            set(a,
                      0x0123456789ABCDEFULL + i,
                      0x0FEDCBA987654321ULL ^ (i * 23),
                      0xAAAAAAAAAAAAAAAAULL + (i * 29),
                      0x5555555555555555ULL ^ (i * 31));
            uint64_t k = 0x9E3779B185EBCA87ULL + (i * 37);

            uint64_t c = fn(r, a, k);
            sink ^= r[1] + c;
        }
        uint64_t t1 = now_ns();
        double ms = (double)(t1 - t0) / 1e6;
        std::printf("[perf] %s: %.3f ms (N=%zu)\n", name, ms, N);
    };

    // --- mp ---
    bench_addmul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return mp_addmul(out, x, k);
    }, "mp_addmul       ");

    // --- reference without __uint128_t ---
    bench_addmul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return ref_addmul_no128(out, x, k);
    }, "ref_addmul_no128");

#if defined(__SIZEOF_INT128__)
    bench_addmul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return ref_addmul_128(out, x, k);
    }, "ref_addmul_128  ");
#else
    std::printf("[perf] __uint128_t not available on this toolchain\n");
#endif

    bench_mul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return mp_mul(out, x, k);
    }, "mp_mul          ");

    bench_mul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return ref_mul_no128(out, x, k);
    }, "ref_mul_no128   ");

#if defined(__SIZEOF_INT128__)
    bench_mul([](mp_uint_t out, const mp_uint_t x, uint64_t k) -> uint64_t {
        return ref_mul_128(out, x, k);
    }, "ref_mul_128     ");
#else
    std::printf("[perf] __uint128_t not available on this toolchain\n");
#endif

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(fq_perf_inv, rawfq_inv_1m) {
    RawFq Fq;

    constexpr size_t N = 1'000'000;

    RawFq::Element a, inv;

    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) {
        Fq.inv(inv, a);
    }

    {
        RawFq::Element prod;
        RawFq::Element one = Fq.one();
        // 1^-1 == 1
        Fq.set(a, 1);
        Fq.inv(inv, a);
        Fq.mul(prod, a, inv);
        ASSERT_EQ(Fq.toString(prod, 16), Fq.toString(one, 16));

        // 2^-1 * 2 == 1
        Fq.set(a, 2);
        Fq.inv(inv, a);
        Fq.mul(prod, a, inv);
        ASSERT_EQ(Fq.toString(prod, 16), Fq.toString(one, 16));
    }

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
        // a = i+1
        Fq.set(a, (int)((i % 1000000) + 1));
        Fq.inv(inv, a);
        sink ^= inv.v[0];
    }
    uint64_t t1 = now_ns();

    double ms = (double)(t1 - t0) / 1e6;
    double ns_per = (double)(t1 - t0) / (double)N;

    std::printf("[perf] RawFq::inv: %.3f ms (N=%zu)  => %.2f ns/op\n", ms, N, ns_per);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

