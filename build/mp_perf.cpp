#include "mp.hpp"
#include "fq.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
#include "gtest/gtest.h"
#include <cstring>

static inline void set(uint64_t r[MP_N64], uint64_t l0, uint64_t l1, uint64_t l2, uint64_t l3) {
    r[0] = l0; r[1] = l1; r[2] = l2; r[3] = l3;
}

static inline void mp_add_mod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    uint64_t carry = mp_add(r, a, b);

    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
    mp_uint_t aa, bb, q, rem;

    if (mp_cmp(a, mod) >= 0) {
        mp_div(q, rem, a, mod);
        mp_copy(aa, rem);
    } else {
        mp_copy(aa, a);
    }

    if (mp_cmp(b, mod) >= 0) {
        mp_div(q, rem, b, mod);
        mp_copy(bb, rem);
    } else {
        mp_copy(bb, b);
    }

    const uint64_t car = mp_add(r, aa, bb);

    if (car || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
        if (mp_cmp(r, mod) >= 0) mp_sub(r, r, mod);
    }
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

static inline void ref_mul_full_no128(uint64_t out2n[2 * MP_N64], const uint64_t a[MP_N64], const uint64_t b[MP_N64]) {
    for (int i = 0; i < 2 * MP_N64; i++) out2n[i] = 0;

    for (int i = 0; i < MP_N64; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < MP_N64; j++) {
            uint64_t lo, hi;
            mul64_no128(a[i], b[j], lo, hi);

            // out2n[i+j] += lo + carry
            uint64_t t = out2n[i + j] + lo;
            uint64_t c1 = (t < out2n[i + j]) ? 1u : 0u;
            t += carry;
            uint64_t c2 = (t < carry) ? 1u : 0u;

            out2n[i + j] = t;

            // carry = hi + c1 + c2  (как в твоих ref_*; возможный оверфлоу >64 игнорим)
            carry = hi + c1 + c2;
        }

        // propagate carry
        int k = i + MP_N64;
        while (carry && k < 2 * MP_N64) {
            uint64_t t = out2n[k] + carry;
            uint64_t c = (t < out2n[k]) ? 1u : 0u;
            out2n[k] = t;
            carry = c;
            k++;
        }
    }
}

#if defined(__SIZEOF_INT128__)
static inline void ref_mul_full_128(uint64_t out2n[2 * MP_N64], const uint64_t a[MP_N64], const uint64_t b[MP_N64]) {
    for (int i = 0; i < 2 * MP_N64; i++) out2n[i] = 0;

    for (int i = 0; i < MP_N64; i++) {
        __uint128_t carry = 0;
        for (int j = 0; j < MP_N64; j++) {
            __uint128_t cur = (__uint128_t)a[i] * (__uint128_t)b[j]
                            + (__uint128_t)out2n[i + j]
                            + carry;
            out2n[i + j] = (uint64_t)cur;
            carry = cur >> 64;
        }

        int k = i + MP_N64;
        while (carry && k < 2 * MP_N64) {
            __uint128_t cur = (__uint128_t)out2n[k] + carry;
            out2n[k] = (uint64_t)cur;
            carry = cur >> 64;
            k++;
        }
    }
}
#endif

static inline void ref_mod_2n_n_no128(uint64_t r[MP_N64], const uint64_t num2n[2 * MP_N64], const uint64_t den[MP_N64]) {
    uint64_t rem[MP_N64 + 1];
    std::memset(rem, 0, sizeof(rem));

    for (int bit = 511; bit >= 0; --bit) {
        // rem <<= 1
        uint64_t carry = 0;
        for (int i = 0; i < MP_N64 + 1; i++) {
            uint64_t nc = rem[i] >> 63;
            rem[i] = (rem[i] << 1) | carry;
            carry = nc;
        }

        // rem |= bit(num2n)
        const uint64_t w = num2n[(unsigned)bit >> 6];
        const uint64_t b = (w >> ((unsigned)bit & 63u)) & 1ULL;
        rem[0] |= b;

        // if rem >= den then rem -= den
        bool ge = false;
        if (rem[MP_N64] != 0) {
            ge = true;
        } else {
            for (int i = MP_N64 - 1; i >= 0; --i) {
                if (rem[i] > den[i]) { ge = true; break; }
                if (rem[i] < den[i]) { ge = false; break; }
                if (i == 0) { ge = true; } // equal
            }
        }

        if (ge) {
            uint64_t br = mp_sub(rem, rem, den);
            rem[MP_N64] -= br;
        }
    }

    mp_copy(r, rem);
}

static inline void ref_mulmod_no128(uint64_t r[MP_N64], const uint64_t a[MP_N64], const uint64_t b[MP_N64], const uint64_t mod[MP_N64]) {
    if (mp_is_zero(mod)) { mp_set(r, 0); return; }

    uint64_t res[MP_N64]; mp_set(res, 0);

    uint64_t cur[MP_N64];
    if (mp_cmp(a, mod) >= 0) {
        uint64_t q[MP_N64], rem[MP_N64];
        mp_div(q, rem, a, mod);
        mp_copy(cur, rem);
    } else {
        mp_copy(cur, a);
    }

    for (int bit = 0; bit < 256; ++bit) {
        if (mp_tstbit(b, (size_t)bit)) {
            uint64_t tmp[MP_N64];
            mp_add_mod(tmp, res, cur, mod);
            mp_copy(res, tmp);
        }

        if (bit != 255) {
            uint64_t tmp2[MP_N64];
            mp_add_mod(tmp2, cur, cur, mod);
            mp_copy(cur, tmp2);
        }
    }

    mp_copy(r, res);
}

static inline void ref_pow_mod_no128(uint64_t r[MP_N64], const uint64_t base[MP_N64], const uint64_t exp[MP_N64], const uint64_t mod[MP_N64]) {
    mp_uint_t one; mp_set(one, 1u);
    if (mp_cmp(mod, one) == 0) { mp_set(r, 0); return; }

    mp_uint_t bcur;
    if (mp_cmp(base, mod) >= 0) {
        mp_uint_t q, rem;
        mp_div(q, rem, base, mod);
        mp_copy(bcur, rem);
    } else {
        mp_copy(bcur, base);
    }

    int topBit = -1;
    for (int limb = MP_N64 - 1; limb >= 0 && topBit < 0; --limb) {
        uint64_t w = exp[limb];
        if (!w) continue;
        for (int bit = 63; bit >= 0; --bit) {
            if ((w >> bit) & 1u) {
                topBit = limb * 64 + bit;
                break;
            }
        }
    }

    if (topBit < 0) {
        if (mp_cmp(one, mod) >= 0) {
            mp_uint_t q, rem;
            mp_div(q, rem, one, mod);
            mp_copy(r, rem);
        } else {
            mp_copy(r, one);
        }
        return;
    }

    mp_uint_t acc;
    mp_copy(acc, bcur);

    for (int i = topBit - 1; i >= 0; --i) {
        mp_uint_t sq;
        ref_mulmod_no128(sq, acc, acc, mod);
        mp_copy(acc, sq);

        const int limb = (i >> 6);
        const int bit  = (i & 63);
        if ((exp[limb] >> bit) & 1u) {
            mp_uint_t tmp;
            ref_mulmod_no128(tmp, acc, bcur, mod);
            mp_copy(acc, tmp);
        }
    }

    mp_copy(r, acc);
}

static inline unsigned clz64_perf(uint64_t x) {
    return x ? (unsigned)__builtin_clzll(x) : 64u;
}

// u[0..n] -= qhat * v[0..n-1]; return borrow (0/1)
// u must have length n+1
static inline uint64_t mul_sub_knuth_perf(uint64_t *u, const uint64_t *v, int n, uint64_t qhat) {
#if !defined(__SIZEOF_INT128__)
    (void)u; (void)v; (void)n; (void)qhat;
    return 0;
#else
    __uint128_t carry = 0;
    uint64_t borrow = 0;

    for (int i = 0; i < n; i++) {
        __uint128_t p = (__uint128_t)qhat * (__uint128_t)v[i] + carry;
        uint64_t plo = (uint64_t)p;
        uint64_t phi = (uint64_t)(p >> 64);

        uint64_t ui = u[i];
        uint64_t t = ui - plo;
        uint64_t b1 = (ui < plo) ? 1u : 0u;

        uint64_t t2 = t - borrow;
        uint64_t b2 = (t < borrow) ? 1u : 0u;

        u[i] = t2;

        borrow = b1 | b2;
        carry = phi;
    }

    // subtract carry and borrow from u[n]
    {
        uint64_t ui = u[n];
        uint64_t clo = (uint64_t)carry;

        uint64_t t = ui - clo;
        uint64_t b1 = (ui < clo) ? 1u : 0u;

        uint64_t t2 = t - borrow;
        uint64_t b2 = (t < borrow) ? 1u : 0u;

        u[n] = t2;
        return (b1 | b2);
    }
#endif
}

static inline void add_back_knuth_perf(uint64_t *u, const uint64_t *v, int n) {
#if !defined(__SIZEOF_INT128__)
    (void)u; (void)v; (void)n;
#else
    __uint128_t carry = 0;
    for (int i = 0; i < n; i++) {
        __uint128_t s = (__uint128_t)u[i] + (__uint128_t)v[i] + carry;
        u[i] = (uint64_t)s;
        carry = s >> 64;
    }
    u[n] += (uint64_t)carry;
#endif
}

// --------- 128-optimized versions (what you call “mp_*”) ---------
static inline void mp_mul_full(uint64_t *out2n, const uint64_t *a, const uint64_t *b) {
#if defined(__SIZEOF_INT128__)
    for (int i = 0; i < 2 * MP_N64; i++) out2n[i] = 0;

    for (int i = 0; i < MP_N64; i++) {
        __uint128_t carry = 0;
        for (int j = 0; j < MP_N64; j++) {
            __uint128_t cur = (__uint128_t)a[i] * (__uint128_t)b[j]
                            + (__uint128_t)out2n[i + j]
                            + carry;
            out2n[i + j] = (uint64_t)cur;
            carry = cur >> 64;
        }

        int k = i + MP_N64;
        while (carry && k < 2 * MP_N64) {
            __uint128_t cur = (__uint128_t)out2n[k] + carry;
            out2n[k] = (uint64_t)cur;
            carry = cur >> 64;
            k++;
        }
    }
#else
    // if no int128, fall back to no128 full mul
    ref_mul_full_no128(out2n, a, b);
#endif
}

static inline void mp_mod_2n_n(uint64_t *r, const uint64_t *num2n, const uint64_t *den) {
#if defined(__SIZEOF_INT128__)
    // Knuth D (n=4) variant, requires den[3]!=0 to be worth it
    if (den[MP_N64 - 1] != 0) {
        uint64_t vnorm[MP_N64];
        uint64_t unorm[2 * MP_N64 + 1];

        unsigned s = clz64_perf(den[MP_N64 - 1]);

        if (s == 0) {
            for (int i = 0; i < MP_N64; i++) vnorm[i] = den[i];
            for (int i = 0; i < 2 * MP_N64; i++) unorm[i] = num2n[i];
            unorm[2 * MP_N64] = 0;
        } else {
            uint64_t carry = 0;
            for (int i = 0; i < MP_N64; i++) {
                uint64_t x = den[i];
                vnorm[i] = (x << s) | carry;
                carry = x >> (64 - s);
            }

            carry = 0;
            for (int i = 0; i < 2 * MP_N64; i++) {
                uint64_t x = num2n[i];
                unorm[i] = (x << s) | carry;
                carry = x >> (64 - s);
            }
            unorm[2 * MP_N64] = carry;
        }

        const uint64_t v1 = vnorm[MP_N64 - 1];
        const uint64_t v2 = vnorm[MP_N64 - 2];

        for (int j = MP_N64; j >= 0; --j) {
            __uint128_t uj2 = ((__uint128_t)unorm[j + MP_N64] << 64)
                            | (__uint128_t)unorm[j + MP_N64 - 1];

            uint64_t qhat = (uint64_t)(uj2 / v1);
            uint64_t rhat = (uint64_t)(uj2 % v1);

            for (;;) {
                __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
                __uint128_t right = ((__uint128_t)rhat << 64) | (__uint128_t)unorm[j + MP_N64 - 2];
                if (left <= right) break;
                qhat--;
                rhat += v1;
                if (rhat < v1) break;
            }

            uint64_t borrow_out = mul_sub_knuth_perf(&unorm[j], vnorm, MP_N64, qhat);
            if (borrow_out) add_back_knuth_perf(&unorm[j], vnorm, MP_N64);
        }

        if (s == 0) {
            for (int i = 0; i < MP_N64; i++) r[i] = unorm[i];
        } else {
            uint64_t carry = 0;
            for (int i = MP_N64 - 1; i >= 0; --i) {
                uint64_t x = unorm[i];
                r[i] = (x >> s) | carry;
                carry = x << (64 - s);
            }
        }

        if (mp_cmp(r, den) >= 0) mp_sub(r, r, den);
        return;
    }
#endif

    // fallback bit-by-bit remainder
    ref_mod_2n_n_no128(r, num2n, den);
}

static inline void mp_mulmod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
#if defined(__SIZEOF_INT128__)
    uint64_t prod[2 * MP_N64];
    mp_mul_full(prod, a, b);
    mp_mod_2n_n(r, prod, mod);
#else
    ref_mulmod_no128(r, a, b, mod);
#endif
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

static inline uint64_t splitmix64_next(uint64_t &state) {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

TEST(fq_perf_inv, rawfq_inv_1m_random) {
    RawFq Fq;
    constexpr size_t N = 1'000'000;

    std::vector<RawFq::Element> inputs(N);

    uint64_t rng = 0x123456789ABCDEF0ULL;

    // prepare inputs
    for (size_t i = 0; i < N; i++) {
        uint64_t t[MP_N64];
        t[0] = splitmix64_next(rng);
        t[1] = splitmix64_next(rng);
        t[2] = splitmix64_next(rng);
        t[3] = splitmix64_next(rng);

        while (mp_cmp(t, Fq_q.longVal) >= 0) {
            mp_sub(t, t, Fq_q.longVal);
        }

        // raw -> mont
        mp_copy(inputs[i].v, t);
        Fq_rawMMul(inputs[i].v, inputs[i].v, Fq_R2.longVal);
    }

    // warmup
    RawFq::Element a, inv;
    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) Fq.inv(inv, a);

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
        Fq.inv(inv, inputs[i]);
        sink ^= inv.v[0];
    }
    uint64_t t1 = now_ns();

    std::printf("[perf] RawFq::inv random: %.3f ms (N=%zu) => %.2f ns/op\n",
                (double)(t1 - t0) / 1e6, N, (double)(t1 - t0) / (double)N);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(fq_perf_inv, rawfq_inv_1m_random_nonzero) {
    RawFq Fq;
    constexpr size_t N = 1'000'000;

    std::vector<RawFq::Element> inputs(N);

    uint64_t rng = 0x0FEDCBA987654321ULL;

    // prepare inputs (nonzero)
    for (size_t i = 0; i < N; i++) {
        uint64_t t[MP_N64];
        do {
            t[0] = splitmix64_next(rng);
            t[1] = splitmix64_next(rng);
            t[2] = splitmix64_next(rng);
            t[3] = splitmix64_next(rng);

            while (mp_cmp(t, Fq_q.longVal) >= 0) {
                mp_sub(t, t, Fq_q.longVal);
            }
        } while (mp_is_zero(t));

        mp_copy(inputs[i].v, t);
        Fq_rawMMul(inputs[i].v, inputs[i].v, Fq_R2.longVal);
    }

    // warmup
    RawFq::Element a, inv;
    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) Fq.inv(inv, a);

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
        Fq.inv(inv, inputs[i]);
        sink ^= inv.v[0];
    }
    uint64_t t1 = now_ns();

    std::printf("[perf] RawFq::inv random nonzero: %.3f ms (N=%zu) => %.2f ns/op\n",
                (double)(t1 - t0) / 1e6, N, (double)(t1 - t0) / (double)N);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(mp_perf_more, mulfull_mod_mulmod_powmod_timing) {
    constexpr size_t N_MULFULL = 300000;
    constexpr size_t N_MOD     = 100000;
    constexpr size_t N_MULMOD  = 80000;
    constexpr size_t N_POW     = 5000;

    mp_uint_t a{}, b{}, r{};
    uint64_t prod[2 * MP_N64];

    volatile uint64_t sink = 0;

    auto prep_ab = [&](size_t i) {
        set(a,
            0x0123456789ABCDEFULL + i,
            0x0FEDCBA987654321ULL ^ (i * 23),
            0xAAAAAAAAAAAAAAAAULL + (i * 29),
            0x5555555555555555ULL ^ (i * 31));

        set(b,
            0x9E3779B97F4A7C15ULL + i,
            0xBF58476D1CE4E5B9ULL ^ (i << 1),
            0x94D049BB133111EBULL + (i * 3),
            0xD6E8FEB86659FD93ULL ^ (i * 7));

        // редуцируем в поле, чтобы не улетать в дивизию постоянно
        while (mp_cmp(a, Fq_q.longVal) >= 0) mp_sub(a, a, Fq_q.longVal);
        while (mp_cmp(b, Fq_q.longVal) >= 0) mp_sub(b, b, Fq_q.longVal);
    };

    auto bench_void = [&](auto &&fn, const char *name, size_t N) {
        uint64_t t0 = now_ns();
        for (size_t i = 0; i < N; i++) {
            prep_ab(i);
            fn(i);
        }
        uint64_t t1 = now_ns();
        std::printf("[perf] %s: %.3f ms (N=%zu)\n", name, (double)(t1 - t0) / 1e6, N);
    };

    // --- mul_full ---
#if defined(__SIZEOF_INT128__)
    bench_void([&](size_t i) {
        (void)i;
        mp_mul_full(prod, a, b);
        sink ^= prod[0] ^ prod[7];
    }, "mp_mul_full            ", N_MULFULL);
#else
    std::printf("[perf] mp_mul_full: __uint128_t not available on this toolchain\n");
#endif

    bench_void([&](size_t i) {
        (void)i;
        ref_mul_full_no128(prod, a, b);
        sink ^= prod[1] ^ prod[6];
    }, "ref_mul_full_no128     ", N_MULFULL);

#if defined(__SIZEOF_INT128__)
    bench_void([&](size_t i) {
        (void)i;
        ref_mul_full_128(prod, a, b);
        sink ^= prod[2] ^ prod[5];
    }, "ref_mul_full_128       ", N_MULFULL);
#endif

    // --- mod_2n_n ---
#if defined(__SIZEOF_INT128__)
    bench_void([&](size_t i) {
        (void)i;
        ref_mul_full_128(prod, a, b);
        mp_mod_2n_n(r, prod, Fq_q.longVal);
        sink ^= r[0];
    }, "mp_mod_2n_n            ", N_MOD);
#else
    std::printf("[perf] mp_mod_2n_n: __uint128_t not available on this toolchain\n");
#endif

    bench_void([&](size_t i) {
        (void)i;
        ref_mul_full_no128(prod, a, b);
        ref_mod_2n_n_no128(r, prod, Fq_q.longVal);
        sink ^= r[1];
    }, "ref_mod_2n_n_no128     ", N_MOD);

    // --- mulmod ---
#if defined(__SIZEOF_INT128__)
    bench_void([&](size_t i) {
        (void)i;
        mp_mulmod(r, a, b, Fq_q.longVal);
        sink ^= r[2];
    }, "mp_mulmod              ", N_MULMOD);
#else
    std::printf("[perf] mp_mulmod: __uint128_t not available on this toolchain\n");
#endif

    bench_void([&](size_t i) {
        (void)i;
        ref_mulmod_no128(r, a, b, Fq_q.longVal);
        sink ^= r[3];
    }, "ref_mulmod_no128       ", N_MULMOD);

    mp_uint_t e{};
    mp_set(e, 0);
    e[0] = 0xD6E8FEB86659FD93ULL; // 64-bit exp

    bench_void([&](size_t i) {
        (void)i;
        mp_pow_mod(r, a, e, Fq_q.longVal);
        sink ^= r[0];
    }, "mp_pow_mod             ", N_POW);

    bench_void([&](size_t i) {
        (void)i;
        ref_pow_mod_no128(r, a, e, Fq_q.longVal);
        sink ^= r[1];
    }, "ref_pow_mod_no128      ", N_POW);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}