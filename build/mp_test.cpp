#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <array>
#include <climits>
#include "gtest/gtest.h"
#include "mp.hpp"

namespace {

typedef uint64_t FqRawElement[MP_N64];

typedef struct __attribute__((__packed__)) {
    int32_t shortVal;
    uint32_t type;
    FqRawElement longVal;
} FqElement;

FqElement Fq_q = {0, 0x80000000, {0x3c208c16d87cfd47,0x97816a916871ca8d,0xb85045b68181585d,0x30644e72e131a029}};

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

// Reference addmul for 256-bit: r += a*b, return carry-out (overflow beyond 256 bits)
static inline uint64_t ref_addmul_no128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t b) {
    uint64_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        uint64_t plo, phi;
        mul64_no128(a[i], b, plo, phi);

        uint64_t t = r[i] + plo;
        uint64_t c1 = (t < r[i]) ? 1u : 0u;
        uint64_t out = t + carry;
        uint64_t c2 = (out < t) ? 1u : 0u;
        r[i] = out;

        uint64_t c = phi;
        c += c1;
        c += c2;
        carry = c;
    }
    return carry;
}

static inline uint64_t ref_addmul_128(uint64_t r[MP_N64], const uint64_t a[MP_N64], uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)r[i] + (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
#else
    (void)r; (void)a; (void)b;
    return 0;
#endif
}

// Reference addmul for 256-bit using __uint128_t
static inline uint64_t ref_add_128(uint64_t r[MP_N64], const uint64_t a[MP_N64], const uint64_t b[MP_N64]) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        __uint128_t s = (__uint128_t)a[i] + (__uint128_t)b[i] + carry;
        r[i] = (uint64_t)s;
        carry = s >> 2*MP_N;
    }
    return (uint64_t)carry;
#else
    (void)r; (void)a; (void)b;
    return 0;
#endif
}

// Reference addmul for 256-bit WITHOUT using __uint128_t
static inline uint64_t ref_add_no128(uint64_t r[MP_N64], const uint64_t a[MP_N64], const uint64_t b[MP_N64]) {
    uint64_t carry = 0;
    for (size_t i = 0; i < MP_N64; i++) {
        uint64_t x = a[i];
        uint64_t y = b[i];
        uint64_t s1 = x + y;
        uint64_t c1 = (s1 < x) ? 1u : 0u;
        uint64_t s2 = s1 + carry;
        uint64_t c2 = (s2 < s1) ? 1u : 0u;
        r[i] = s2;
        carry = (c1 | c2);
    }
    return carry;
}

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

// Reduce x (0<=x<2^256) modulo mod using mp_div remainder
static inline void reduce_mod(uint64_t r[MP_N64], const uint64_t x[MP_N64], const uint64_t mod[MP_N64]) {
    uint64_t q[MP_N64], rem[MP_N64];
    mp_div(q, rem, x, mod);
    mp_copy(r, rem);
}

} // namespace


TEST(mp_set_u64, mp_set) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_TRUE(mp_is_zero(a));

    mp_set(a, 123456789ULL);
    EXPECT_EQ(a[0], 123456789ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0ULL);
}

TEST(mp_copy, self_noop) {
    mp_uint_t a{};
    set(a, 1, 2, 3, 4);
    mp_copy(a, a); // must be no-op
    EXPECT_EQ(a[0], 1ULL);
    EXPECT_EQ(a[1], 2ULL);
    EXPECT_EQ(a[2], 3ULL);
    EXPECT_EQ(a[3], 4ULL);
}

TEST(mp_copy, non_overlapping) {
    mp_uint_t a{}, b{};
    set(a, 1, 2, 3, 4);
    mp_copy(b, a);
    EXPECT_EQ(mp_cmp(a, b), 0);
}

TEST(mp_copy, overlapping_forward_memmove_semantics) {
    // r overlaps a and starts inside source range: requires backward copy
    // Layout: buf[0..4], copy 4 limbs from &buf[0] to &buf[1]
    uint64_t buf[MP_N64 + 1] = {11, 22, 33, 44, 55};

    uint64_t *a = &buf[0];
    uint64_t *r = &buf[1];

    mp_copy(r, a);

    // Expected like memmove: r[0..3] = old a[0..3]
    EXPECT_EQ(buf[0], 11ULL);
    EXPECT_EQ(buf[1], 11ULL);
    EXPECT_EQ(buf[2], 22ULL);
    EXPECT_EQ(buf[3], 33ULL);
    EXPECT_EQ(buf[4], 44ULL);
}

TEST(mp_copy, overlapping_backward_memmove_semantics) {
    // r overlaps a but starts before source: forward copy is OK
    uint64_t buf[MP_N64 + 1] = {11, 22, 33, 44, 55};

    uint64_t *a = &buf[1];
    uint64_t *r = &buf[0];

    mp_copy(r, a);

    EXPECT_EQ(buf[0], 22ULL);
    EXPECT_EQ(buf[1], 33ULL);
    EXPECT_EQ(buf[2], 44ULL);
    EXPECT_EQ(buf[3], 55ULL);
    EXPECT_EQ(buf[4], 55ULL); // unchanged tail
}

TEST(mp_cmp, mp_cmp) {
    mp_uint_t a{}, b{};
    mp_set(a, 0);
    mp_set(b, 0);
    EXPECT_EQ(mp_cmp(a, b), 0);

    set(a, 1, 0, 0, 0);
    set(b, 2, 0, 0, 0);
    EXPECT_LT(mp_cmp(a, b), 0);
    EXPECT_GT(mp_cmp(b, a), 0);

    // different high limb
    set(a, 0, 0, 0, 1);
    set(b, 0, 0, 0, 2);
    EXPECT_LT(mp_cmp(a, b), 0);
}

TEST(mp_is_zero, mp_is_zero) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_TRUE(mp_is_zero(a));
    a[2] = 1;
    EXPECT_FALSE(mp_is_zero(a));
}

TEST(mp_add_256, mp_add) {
    const std::array<std::array<uint64_t, MP_N64>, 6> cases_a = {{
        {{0,0,0,0}},
        {{1,0,0,0}},
        {{0xFFFFFFFFFFFFFFFFULL,0,0,0}},
        {{0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0,0}},
        {{0,0,0,0x8000000000000000ULL}},
        {{0x0123456789ABCDEFULL,0x0,0xFFFFFFFFFFFFFFFFULL,0x7ULL}},
    }};

    for (size_t i = 0; i < cases_a.size(); i++) {
        const std::array<std::array<uint64_t, MP_N64>, 6> cases_b = {{
            {{0,0,0,0}},
            {{2,0,0,0}},
            {{1,0,0,0}},
            {{1,0,0,0}},
            {{0,0,0,0x8000000000000000ULL}},
            {{0x1111111111111111ULL,0x2222222222222222ULL,0,0}},
        }};
        mp_uint_t r{}, ref1{}, ref2{};
        uint64_t c = mp_add(r, cases_a[i].data(), cases_b[i].data());
        uint64_t c1 = ref_add_no128(ref1, cases_a[i].data(), cases_b[i].data());
        EXPECT_EQ(c, c1);
        EXPECT_EQ(mp_cmp(r, ref1), 0);

#if defined(__SIZEOF_INT128__)
        uint64_t c2 = ref_add_128(ref2, cases_a[i].data(), cases_b[i].data());
        EXPECT_EQ(c, c2);
        EXPECT_EQ(mp_cmp(r, ref2), 0);
#endif
    }

    // boundary: max + 1 => 0 carry 1
    mp_uint_t max{}, one{}, out{};
    set(max, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
    mp_set(one, 1);
    EXPECT_EQ(mp_add(out, max, one), 1ULL);
    EXPECT_TRUE(mp_is_zero(out));
}

TEST(mp_add_varlimbs, mp_add) {
    {
        uint64_t a1[1] = {0xFFFFFFFFFFFFFFFFULL};
        uint64_t b1[1] = {1};
        uint64_t r1[1] = {0};
        EXPECT_EQ(mp_add(r1, a1, 1, b1, 1), 1ULL);
        EXPECT_EQ(r1[0], 0ULL);
    }
    {
        uint64_t a2[2] = {0xFFFFFFFFFFFFFFFFULL, 0x0ULL};
        uint64_t b1[1] = {1};
        uint64_t r2[2] = {0,0};
        EXPECT_EQ(mp_add(r2, a2, 2, b1, 1), 0ULL);
        EXPECT_EQ(r2[0], 0ULL);
        EXPECT_EQ(r2[1], 1ULL);
    }
    {
        uint64_t a3[3] = {1,2,3};
        uint64_t b4[4] = {4,5,6,7};
        uint64_t r4[4] = {0,0,0,0};
        EXPECT_EQ(mp_add(r4, a3, 3, b4, 4), 0ULL);
        EXPECT_EQ(r4[0], 5ULL);
        EXPECT_EQ(r4[1], 7ULL);
        EXPECT_EQ(r4[2], 9ULL);
        EXPECT_EQ(r4[3], 7ULL);
    }
    // carry ripple across all provided limbs
    {
        uint64_t a4[4] = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
        uint64_t b1[1] = {1};
        uint64_t r4[4] = {0,0,0,0};
        EXPECT_EQ(mp_add(r4, a4, 4, b1, 1), 1ULL);
        EXPECT_EQ(r4[0], 0ULL);
        EXPECT_EQ(r4[1], 0ULL);
        EXPECT_EQ(r4[2], 0ULL);
        EXPECT_EQ(r4[3], 0ULL);
    }
    // an=0 should behave like adding 0 + b
    {
        uint64_t a0[1] = {0}; // not used
        uint64_t b4[4] = {1,2,3,4};
        uint64_t r4[4] = {0xAA,0xAA,0xAA,0xAA};
        EXPECT_EQ(mp_add(r4, a0, 0, b4, 4), 0ULL);
        EXPECT_EQ(r4[0], 1ULL);
        EXPECT_EQ(r4[1], 2ULL);
        EXPECT_EQ(r4[2], 3ULL);
        EXPECT_EQ(r4[3], 4ULL);
    }
}

TEST(mp_sub_256, mp_sub) {
    mp_uint_t a{}, b{}, r{};

    // 0 - 1 => borrow 1, result = 2^256-1
    mp_set(a, 0);
    mp_set(b, 1);
    EXPECT_EQ(mp_sub(r, a, b), 1ULL);
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);

    // max - max = 0
    set(a, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
    mp_copy(b, a);
    EXPECT_EQ(mp_sub(r, a, b), 0ULL);
    EXPECT_TRUE(mp_is_zero(r));

    // borrow ripple across all limbs: [0,0,0,1] - 1 = [~0,~0,~0,0]
    set(a, 0,0,0,1);
    mp_set(b, 1);
    EXPECT_EQ(mp_sub(r, a, b), 0ULL);
    EXPECT_EQ(r[0], ~0ULL);
    EXPECT_EQ(r[1], ~0ULL);
    EXPECT_EQ(r[2], ~0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_add_u64, mp_add) {
    mp_uint_t a{}, r{};
    set(a, 0xFFFFFFFFFFFFFFFFULL, 0, 0, 0);
    EXPECT_EQ(mp_add(r, a, 1ULL), 0ULL);
    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);

    // overflow
    set(a, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(mp_add(r, a, 1ULL), 1ULL);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_sub_u64, mp_sub) {
    mp_uint_t a{}, r{};
    mp_set(a, 0);
    EXPECT_EQ(mp_sub(r, a, 1ULL), 1ULL);
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);

    set(a, 0, 1, 0, 0);
    EXPECT_EQ(mp_sub(r, a, 1ULL), 0ULL);
    EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[1], 0ULL);
}

TEST(mp_mul_u64, mp_mul) {
    const std::array<std::array<uint64_t, MP_N64>, 5> cases_a = {{
        {{0,0,0,0}},
        {{1,0,0,0}},
        {{0xFFFFFFFFFFFFFFFFULL,0,0,0}},
        {{0,0xFFFFFFFFFFFFFFFFULL,0,0}},
        {{0x0123456789ABCDEFULL,0x0ULL,0xFFFFFFFFFFFFFFFFULL,0x7ULL}},
    }};

    for (size_t i = 0; i < cases_a.size(); i++) {
        const std::array<uint64_t, 5> cases_b = {{0ULL, 2ULL, 3ULL, 5ULL, 0xFEDCBA9876543211ULL}};
        mp_uint_t r{}, ref_no128{}, ref_128{};
        uint64_t c = mp_mul(r, cases_a[i].data(), cases_b[i]);

        uint64_t c_no128 = ref_mul_no128(ref_no128, cases_a[i].data(), cases_b[i]);
        EXPECT_EQ(c, c_no128);
        EXPECT_EQ(mp_cmp(r, ref_no128), 0);

#if defined(__SIZEOF_INT128__)
        uint64_t c_128 = ref_mul_128(ref_128, cases_a[i].data(), cases_b[i]);
        EXPECT_EQ(c, c_128);
        EXPECT_EQ(mp_cmp(r, ref_128), 0);
#endif
    }
    // boundary: (2^256-1) * 2 = (2^256-2) with carry=1
    {
        mp_uint_t a{}, r{}, ref{};
        set(a, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
        uint64_t c = mp_mul(r, a, 2);
        mp_copy(ref, a);
        mp_shl(ref, ref, 1); // mod 2^256
        EXPECT_EQ(c, 1ULL);
        EXPECT_EQ(mp_cmp(r, ref), 0);
    }
}

TEST(mp_addmul_u64, mp_addmul) {
    // r += a*b
    mp_uint_t a{}, r{}, ref{};
    set(a, 3, 0, 0, 0);
    set(r, 5, 0, 0, 0);

    // expected: 5 + 3*7 = 26
    mp_copy(ref, r);
    mp_uint_t prod{};
    mp_mul(prod, a, 7);
    mp_add(ref, ref, prod);

    mp_addmul(r, a, 7);
    EXPECT_EQ(mp_cmp(r, ref), 0);
    // also check carry-out and overflow case
    {
        mp_uint_t aa{}, rr{}, rr_ref{};
        set(aa, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
        set(rr, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
        mp_copy(rr_ref, rr);

        uint64_t c_ref = ref_addmul_no128(rr_ref, aa, 2);
        uint64_t c = mp_addmul(rr, aa, 2);

        EXPECT_EQ(c, c_ref);
        EXPECT_EQ(mp_cmp(rr, rr_ref), 0);

#if defined(__SIZEOF_INT128__)
        mp_uint_t rr_ref2{};
        set(rr_ref2, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
        uint64_t c_ref2 = ref_addmul_128(rr_ref2, aa, 2);
        EXPECT_EQ(c, c_ref2);
        EXPECT_EQ(mp_cmp(rr, rr_ref2), 0);
#endif
    }
}

TEST(mp_and, mp_and) {
    mp_uint_t a{}, b{}, r{};
    set(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_and(r, a, b);
    EXPECT_EQ(r[0], a[0] & b[0]);
    EXPECT_EQ(r[1], a[1] & b[1]);
    EXPECT_EQ(r[2], a[2] & b[2]);
    EXPECT_EQ(r[3], a[3] & b[3]);
}

TEST(mp_or, mp_or) {
    mp_uint_t a{}, b{}, r{};
    set(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_or(r, a, b);
    EXPECT_EQ(r[0], a[0] | b[0]);
    EXPECT_EQ(r[1], a[1] | b[1]);
    EXPECT_EQ(r[2], a[2] | b[2]);
    EXPECT_EQ(r[3], a[3] | b[3]);
}

TEST(mp_xor, mp_xor) {
    mp_uint_t a{}, b{}, r{};
    set(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_xor(r, a, b);
    EXPECT_EQ(r[0], a[0] ^ b[0]);
    EXPECT_EQ(r[1], a[1] ^ b[1]);
    EXPECT_EQ(r[2], a[2] ^ b[2]);
    EXPECT_EQ(r[3], a[3] ^ b[3]);
}

TEST(mp_not, mp_not) {
    mp_uint_t a{}, r{};
    set(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    mp_not(r, a);
    EXPECT_EQ(r[0], ~a[0]);
    EXPECT_EQ(r[1], ~a[1]);
    EXPECT_EQ(r[2], ~a[2]);
    EXPECT_EQ(r[3], ~a[3]);
}

TEST(mp_tstbit, mp_tstbit) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_FALSE(mp_tstbit(a, 0));
    a[0] = 1;
    EXPECT_TRUE(mp_tstbit(a, 0));
    EXPECT_FALSE(mp_tstbit(a, 1));

    mp_set(a, 0);
    a[1] = 1ULL << 63; // bit 127
    EXPECT_TRUE(mp_tstbit(a, 127));
    EXPECT_FALSE(mp_tstbit(a, 126));

    mp_set(a, 0);
    a[3] = 1ULL << 63; // bit 255
    EXPECT_TRUE(mp_tstbit(a, 255));
    EXPECT_FALSE(mp_tstbit(a, 254));

    // boundary bits: 63 and 64
    mp_set(a, 0);
    a[0] = 1ULL << 63;
    EXPECT_TRUE(mp_tstbit(a, 63));
    EXPECT_FALSE(mp_tstbit(a, 62));

    mp_set(a, 0);
    a[1] = 1ULL; // bit 64
    EXPECT_TRUE(mp_tstbit(a, 64));
    EXPECT_FALSE(mp_tstbit(a, 65));

    // out of range
    EXPECT_FALSE(mp_tstbit(a, 256));
}

// -------------------- mp_shl --------------------

TEST(mp_shl, identity_shift0) {
    mp_uint_t a{}, r{};
    set(a, 0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL);
    mp_shl(r, a, 0);
    EXPECT_EQ(mp_cmp(r, a), 0);
}

TEST(mp_shl, zero_stays_zero) {
    mp_uint_t a{}, r{};
    mp_set(a, 0);

    mp_shl(r, a, 1);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shl(r, a, 255);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shl(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shl(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shl, basic_small_shift) {
    mp_uint_t a{}, r{};
    set(a, 1, 0, 0, 0);
    mp_shl(r, a, 1);

    EXPECT_EQ(r[0], 2ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_63) {
    mp_uint_t a{}, r{};
    // a = 3 => 3<<63 has bit63 and bit64 set:
    // limb0 = 0x8000.., limb1 = 1
    set(a, 3ULL, 0ULL, 0ULL, 0ULL);
    mp_shl(r, a, 63);

    EXPECT_EQ(r[0], 0x8000000000000000ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_64_bitShift0_path) {
    mp_uint_t a{}, r{};
    // a = 1 => 1<<64 moves into limb1 (bitShift==0 path)
    set(a, 1ULL, 0ULL, 0ULL, 0ULL);
    mp_shl(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_65) {
    mp_uint_t a{}, r{};
    // a = 1 => 1<<65 => limb1 bit1 => 2
    set(a, 1ULL, 0ULL, 0ULL, 0ULL);
    mp_shl(r, a, 65);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 2ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, shift_across_limbs_example) {
    mp_uint_t a{}, r{};
    // a = 1 in limb1 => <<64 moves into limb2
    set(a, 0ULL, 1ULL, 0ULL, 0ULL);
    mp_shl(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 1ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, shift_255_bit0_to_top_bit) {
    mp_uint_t a{}, r{};
    // a = 1 => 1<<255 => top bit (bit63 of limb3)
    set(a, 1ULL, 0ULL, 0ULL, 0ULL);
    mp_shl(r, a, 255);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0x8000000000000000ULL);
}

TEST(mp_shl, shift_ge_256_is_zero) {
    mp_uint_t a{}, r{};
    set(a,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL);

    mp_shl(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shl(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shl, all_ones_signature_shift1) {
    mp_uint_t a{}, r{};
    // (2^256-1)<<1 mod 2^256 = 2^256-2
    // => limbs: [..FE, ..FF, ..FF, ..FF]
    set(a,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL);

    mp_shl(r, a, 1);

    EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFEULL);
    EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[3], 0xFFFFFFFFFFFFFFFFULL);
}

// -------------------- mp_shr --------------------

TEST(mp_shr, identity_shift0) {
    mp_uint_t a{}, r{};
    set(a, 0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL);
    mp_shr(r, a, 0);
    EXPECT_EQ(mp_cmp(r, a), 0);
}

TEST(mp_shr, zero_stays_zero) {
    mp_uint_t a{}, r{};
    mp_set(a, 0);

    mp_shr(r, a, 1);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shr(r, a, 255);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shr(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shr(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shr, basic_small_shift) {
    mp_uint_t a{}, r{};
    set(a, 2, 0, 0, 0);
    mp_shr(r, a, 1);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_63) {
    mp_uint_t a{}, r{};
    // a = 2^63 + 2^64  => a >> 63 = 3
    set(a, 0x8000000000000000ULL, 0x1ULL, 0x0ULL, 0x0ULL);
    mp_shr(r, a, 63);

    EXPECT_EQ(r[0], 3ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_64_bitShift0_path) {
    mp_uint_t a{}, r{};
    // a = 2^64 => >>64 = 1 (bitShift==0 path)
    set(a, 0x0ULL, 0x1ULL, 0x0ULL, 0x0ULL);
    mp_shr(r, a, 64);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_65) {
    mp_uint_t a{}, r{};
    // a = 2^64 => >>65 = 0
    set(a, 0x0ULL, 0x1ULL, 0x0ULL, 0x0ULL);
    mp_shr(r, a, 65);

    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shr, shift_across_limbs_example) {
    mp_uint_t a{}, r{};
    set(a, 0, 0, 1, 0);
    mp_shr(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, shift_255_top_bit_to_bit0) {
    mp_uint_t a{}, r{};
    // a = 2^255 => a >> 255 = 1
    set(a, 0x0ULL, 0x0ULL, 0x0ULL, 0x8000000000000000ULL);
    mp_shr(r, a, 255);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, shift_ge_256_is_zero) {
    mp_uint_t a{}, r{};
    set(a,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL);

    mp_shr(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shr(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shr, all_ones_signature_shift1) {
    mp_uint_t a{}, r{};
    // (2^256-1)>>1 = 2^255-1
    set(a,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
              0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL);

    mp_shr(r, a, 1);

    EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[3], 0x7FFFFFFFFFFFFFFFULL);
}

TEST(mp_get_int32, mp_get_int32) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_EQ(mp_get_int32(a), 0);

    mp_set(a, 123);
    EXPECT_EQ(mp_get_int32(a), 123);

    set(a, 0xFFFFFFFFULL, 1, 0, 0);
    // mp_get_int32 returns low 32-bit signed value of a[0]
    EXPECT_EQ(mp_get_int32(a), -1);

    set(a, 0x80000000ULL, 0, 0, 0);
    EXPECT_EQ(mp_get_int32(a), (int32_t)0x80000000u);
}

TEST(mp_fits_int32, mp_fits_int32) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX + 1ULL);
    EXPECT_FALSE(mp_fits_int32(a));

    set(a, 1, 1, 0, 0);
    EXPECT_FALSE(mp_fits_int32(a));
}

TEST(mp_set_str, valid_cases) {

    mp_uint_t a{};

    ASSERT_TRUE(mp_set(a, "0", 10));
    EXPECT_TRUE(mp_is_zero(a));

    ASSERT_TRUE(mp_set(a, "1", 16));
    EXPECT_EQ(a[0], 1ULL);

    // leading spaces and leading '+'
    ASSERT_TRUE(mp_set(a, "   +ff", 16));
    EXPECT_EQ(a[0], 255ULL);

    // base=0 => defaults to 10
    ASSERT_TRUE(mp_set(a, "42", 0));
    EXPECT_EQ(a[0], 42ULL);

    // max 256-bit value in hex
    ASSERT_TRUE(mp_set(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16));
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(a[i], 0xFFFFFFFFFFFFFFFFULL);

    // base 10 boundary: 2^64-1
    ASSERT_TRUE(mp_set(a, "18446744073709551615", 10));
    EXPECT_EQ(a[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(a[1], 0ULL);
}

TEST(mp_set_str, invalid_cases) {
#if !defined(__SIZEOF_INT128__)
    GTEST_SKIP() << "mp_set(str) requires __int128 in this build";
#endif

    mp_uint_t a{};

    // invalid base
    EXPECT_FALSE(mp_set(a, "1", 1));
    EXPECT_FALSE(mp_set(a, "1", 17));

    // negative not allowed
    EXPECT_FALSE(mp_set(a, "-1", 10));

    // invalid digit for base
    EXPECT_FALSE(mp_set(a, "2", 2));
    EXPECT_FALSE(mp_set(a, "g", 16));

    // empty / only spaces
    EXPECT_FALSE(mp_set(a, "", 10));
    EXPECT_FALSE(mp_set(a, "   ", 10));
    EXPECT_FALSE(mp_set(a, "+", 10));

    // trailing garbage
    EXPECT_FALSE(mp_set(a, "123x", 10));
    EXPECT_FALSE(mp_set(a, "ff zz", 16));

    // overflow: 2^256 is 1 followed by 64 hex zeros => must fail
    EXPECT_FALSE(mp_set(a, "10000000000000000000000000000000000000000000000000000000000000000", 16));
}

TEST(mp_get_str, basic_cases) {

    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_EQ(mp_get_str(a, 10), std::string("0"));

    mp_set(a, 255);
    EXPECT_EQ(mp_get_str(a, 16), std::string("ff"));
    EXPECT_EQ(mp_get_str(a, 2), std::string("11111111"));

    // max value in hex should be 64 f's
    set(a, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
    EXPECT_EQ(mp_get_str(a, 16), std::string(64, 'f'));

    // normalization: leading zeros in input -> output without leading zeros
    mp_uint_t b{};
    ASSERT_TRUE(mp_set(b, "00000100", 16));
    EXPECT_EQ(mp_get_str(b, 16), std::string("100"));
}

TEST(mp_get_str, invalid_base_returns_empty) {

    mp_uint_t a{};
    mp_set(a, 123);

    EXPECT_EQ(mp_get_str(a, 1), std::string());
    EXPECT_EQ(mp_get_str(a, 17), std::string());
}

TEST(mp_set_mod_i64, mp_set_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t r{}, tmp{}, q{}, rem{};

    mp_set_mod(r, 0, mod);
    EXPECT_TRUE(mp_is_zero(r));

    // positive: r = 5
    mp_set_mod(r, 5, mod);
    mp_set(tmp, 5);
    reduce_mod(tmp, tmp, mod);
    EXPECT_EQ(mp_cmp(r, tmp), 0);

    // negative: r = mod - 5
    mp_set_mod(r, -5, mod);
    mp_set(tmp, 5);
    reduce_mod(rem, tmp, mod);
    mp_sub(tmp, mod, rem);
    EXPECT_EQ(mp_cmp(r, tmp), 0);

    // large positive: (mod + 1) mod mod = 1
    // Construct x = mod + 1 (mod is 256-bit) and pass as int64? cannot.
    // Instead test int64 boundary: INT64_MIN/INT64_MAX.
    mp_set_mod(r, INT64_MAX, mod);
    mp_set(tmp, INT64_MAX);
    reduce_mod(tmp, tmp, mod);
    EXPECT_EQ(mp_cmp(r, tmp), 0);

    (void)q; (void)rem;
}

TEST(mp_set_mod_i64, small_mod_extremes) {
    mp_uint_t mod{}, r{};
    mp_set(mod, 17);

    auto expect_i64 = [&](int64_t x) {
        int64_t e = x % 17;
        if (e < 0) e += 17;
        mp_uint_t got{};
        mp_set_mod(got, x, mod);
        EXPECT_EQ(got[0], (uint64_t)e);
        EXPECT_EQ(got[1], 0ULL);
        EXPECT_EQ(got[2], 0ULL);
        EXPECT_EQ(got[3], 0ULL);
    };

    expect_i64(0);
    expect_i64(1);
    expect_i64(-1);
    expect_i64(INT64_MAX);
    expect_i64(INT64_MIN);

    // mod = 1 => always 0
    mp_uint_t mod1{};
    mp_set(mod1, 1);
    mp_set_mod(r, INT64_MIN, mod1);
    EXPECT_TRUE(mp_is_zero(r));
    mp_set_mod(r, INT64_MAX, mod1);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_set_mod_str, mp_set_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t r{}, ref{};

    ASSERT_TRUE(mp_set_mod(r, "0", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    // hex value bigger than mod: parse then reduce via mp_div
    ASSERT_TRUE(mp_set_mod(r,
        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16, mod));

    mp_uint_t max{};
    set(max, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
    reduce_mod(ref, max, mod);
    EXPECT_EQ(mp_cmp(r, ref), 0);

    // negative
    ASSERT_TRUE(mp_set_mod(r, "-1", 10, mod));
    mp_uint_t one{};
    mp_set(one, 1);
    mp_sub(ref, mod, one);
    EXPECT_EQ(mp_cmp(r, ref), 0);

    // "-0" must be 0
    ASSERT_TRUE(mp_set_mod(r, "-0", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    // base=0 defaults to 10
    ASSERT_TRUE(mp_set_mod(r, "15", 0, mod));
    mp_set(ref, 15);
    reduce_mod(ref, ref, mod);
    EXPECT_EQ(mp_cmp(r, ref), 0);
}

TEST(mp_set_mod_str, digit_must_be_lt_base) {
    mp_uint_t mod{}, r{};
    mp_set(mod, 17);

    // digit '2' is invalid in base 2
    EXPECT_FALSE(mp_set_mod(r, "2", 2, mod));
    EXPECT_FALSE(mp_set_mod(r, "10 1", 2, mod)); // trailing garbage after spaces

    // invalid base
    EXPECT_FALSE(mp_set_mod(r, "1", 1, mod));
    EXPECT_FALSE(mp_set_mod(r, "1", 17, mod));

    // empty / spaces only
    EXPECT_FALSE(mp_set_mod(r, "", 10, mod));
    EXPECT_FALSE(mp_set_mod(r, "   ", 10, mod));
    EXPECT_FALSE(mp_set_mod(r, "+", 10, mod));

    // exact multiple of mod: "34" mod 17 = 0
    ASSERT_TRUE(mp_set_mod(r, "34", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    // negative: "-1" mod 17 = 16
    ASSERT_TRUE(mp_set_mod(r, "-1", 10, mod));
    EXPECT_EQ(r[0], 16ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_export_be, mp_export_be) {
    mp_uint_t a{};
    set(a, 0x1122334455667788ULL, 0x99AABBCCDDEEFF00ULL, 0x0102030405060708ULL, 0x1112131415161718ULL);

    uint8_t out[MP_N] = {0};
    mp_export_be(out, a);

    // Big-endian: first byte is MSB of limb3
    EXPECT_EQ(out[0], (uint8_t)0x11);
    EXPECT_EQ(out[1], (uint8_t)0x12);
    EXPECT_EQ(out[2], (uint8_t)0x13);
    EXPECT_EQ(out[3], (uint8_t)0x14);
    EXPECT_EQ(out[31], (uint8_t)0x88);
}

TEST(mp_import_be, mp_import_be) {
    uint8_t in[MP_N] = {0};
    in[31] = 0x01;

    mp_uint_t a{};
    mp_import_be(a, in);
    EXPECT_EQ(a[0], 1ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0ULL);
}

TEST(mp_import_be, msb_sets_bit255) {
    uint8_t in[MP_N] = {0};
    in[0] = 0x80; // MSB of whole 256-bit number
    mp_uint_t a{};
    mp_import_be(a, in);
    EXPECT_EQ(a[0], 0ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0x8000000000000000ULL);
}

TEST(mp_export_import_be, roundtrip) {
    mp_uint_t a{}, b{};
    set(a, 0x1122334455667788ULL, 0x99AABBCCDDEEFF00ULL, 0x0102030405060708ULL, 0x1112131415161718ULL);

    uint8_t buf[MP_N] = {0};
    mp_export_be(buf, a);
    mp_import_be(b, buf);

    EXPECT_EQ(mp_cmp(a, b), 0);
}

TEST(mp_export_import_be, all_ones_roundtrip) {
    mp_uint_t a{}, b{};
    set(a, ~0ULL, ~0ULL, ~0ULL, ~0ULL);

    uint8_t buf[MP_N] = {0};
    mp_export_be(buf, a);
    mp_import_be(b, buf);

    EXPECT_EQ(mp_cmp(a, b), 0);
}

TEST(mp_div, num_lt_den_q0_rnum) {
    mp_uint_t den{}, num{}, q{}, r{};
    set(den, 0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0);
    mp_set(num, 123);

    mp_div(q, r, num, den);

    mp_uint_t qexp{}, rexp{};
    mp_set(qexp, 0);
    mp_set(rexp, 123);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, num_eq_den_q1_r0) {
    mp_uint_t den{}, num{}, q{}, r{}, one{}, zero{};
    set(den, 0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0);
    mp_copy(num, den);

    mp_div(q, r, num, den);

    mp_set(one, 1);
    mp_set(zero, 0);
    EXPECT_EQ(mp_cmp(q, one), 0);
    EXPECT_EQ(mp_cmp(r, zero), 0);
}

TEST(mp_div, den_1word_path) {
    mp_uint_t den{}, num{}, q{}, r{}, qexp{}, rexp{};

    mp_set(den, 7);
    mp_set(num, 1000);

    mp_div(q, r, num, den);

    mp_set(qexp, 142); // 1000/7
    mp_set(rexp, 6);   // 1000%7
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, knuth_path_s_nonzero) {
    mp_uint_t den{}, qv{}, rv{}, num{}, tmp{}, q{}, r{};

    // top limb NOT having MSB set => clz != 0 (s != 0)
    set(den, 0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0);

    mp_set(qv, 17);
    mp_mul(tmp, den, 17);
    mp_set(rv, 5);
    mp_add(num, tmp, rv);

    mp_div(q, r, num, den);

    mp_uint_t qexp{}, rexp{};
    mp_set(qexp, 17);
    mp_set(rexp, 5);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, knuth_path_s_zero) {
    mp_uint_t den{}, num{}, tmp{}, q{}, r{};

    // top limb HAS MSB set => clz == 0 (s == 0)
    set(den, 0x0123456789ABCDEFULL, 0x8000000000000000ULL, 0, 0);

    mp_mul(tmp, den, 9);
    mp_uint_t seven{};
    mp_set(seven, 7);
    mp_add(num, tmp, seven);

    mp_div(q, r, num, den);

    mp_uint_t qexp{}, rexp{};
    mp_set(qexp, 9);
    mp_set(rexp, 7);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, den_one_returns_num) {
    mp_uint_t den{}, num{}, q{}, r{};
    mp_set(den, 1);
    set(num, 0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL);

    mp_div(q, r, num, den);

    EXPECT_EQ(mp_cmp(q, num), 0);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_div, top_bit_quotient_case) {
    // num = 2^255, den = 2 => q = 2^254, r = 0
    mp_uint_t num{}, den{}, q{}, r{}, qexp{}, z{};
    set(num, 0,0,0, 0x8000000000000000ULL);
    mp_set(den, 2);
    mp_div(q, r, num, den);
    mp_set(z, 0);
    set(qexp, 0,0,0, 0x4000000000000000ULL);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, z), 0);
}

TEST(mp_pow_mod, exp_zero_returns_one_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t base{}, exp{}, out{}, ref{};
    mp_set(base, 5);
    mp_set(exp, 0);

    mp_pow_mod(out, base, exp, mod);

    mp_set(ref, 1);
    reduce_mod(ref, ref, mod);
    EXPECT_EQ(mp_cmp(out, ref), 0);
}

TEST(mp_pow_mod, exp_one_returns_base_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t base{}, exp{}, out{}, ref{};
    mp_set(base, 5);
    mp_set(exp, 1);

    mp_pow_mod(out, base, exp, mod);

    mp_set(ref, 5);
    reduce_mod(ref, ref, mod);
    EXPECT_EQ(mp_cmp(out, ref), 0);
}

TEST(mp_pow_mod, base_zero_exp_zero_defined_as_one_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t base{}, exp{}, out{}, ref{};
    mp_set(base, 0);
    mp_set(exp, 0);

    mp_pow_mod(out, base, exp, mod);

    mp_set(ref, 1);
    reduce_mod(ref, ref, mod);
    EXPECT_EQ(mp_cmp(out, ref), 0);
}

TEST(mp_pow_mod, mod_one_returns_zero) {
    mp_uint_t base{}, exp{}, mod1{}, out{};
    mp_set(base, 123);
    mp_set(exp, 456);
    mp_set(mod1, 1);
    mp_pow_mod(out, base, exp, mod1);
    EXPECT_TRUE(mp_is_zero(out));
}

TEST(mp_pow_mod, base_reduction_path_base_ge_mod) {
    // Use small mod so we can reason easily
    mp_uint_t mod{}, base{}, exp{}, out{}, ref{};
    mp_set(mod, 17);

    // base = 17+5 => should reduce to 5
    mp_set(base, 22);
    mp_set(exp, 3);

    mp_pow_mod(out, base, exp, mod);

    // 5^3 mod 17 = 125 mod 17 = 6
    mp_set(ref, 6);
    EXPECT_EQ(mp_cmp(out, ref), 0);
}

TEST(mp_inv_mod, zero_not_invertible) {
    const uint64_t *mod = Fq_q.longVal;
    mp_uint_t a{}, inv{};
    mp_set(a, 0);
    EXPECT_FALSE(mp_inv_mod(inv, a, mod));
}

TEST(mp_inv_mod, one_inverts_to_one) {
    const uint64_t *mod = Fq_q.longVal;
    mp_uint_t a{}, inv{}, one{};
    mp_set(a, 1);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(inv, one), 0);
}

TEST(mp_inv_mod, check_a_times_inv_is_one_small_mod) {
#if !defined(__SIZEOF_INT128__)
    GTEST_SKIP() << "mp_inv_mod/mp_div/mod reduction tests require __int128 in this build";
#endif
    mp_uint_t mod{}, a{}, inv{}, one{}, prod{}, rem{};
    mp_set(mod, 17);
    mp_set(a, 3);

    ASSERT_TRUE(mp_inv_mod(inv, a, mod));

    mp_mul(prod, inv, 3);
    reduce_mod(rem, prod, mod);

    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(rem, one), 0);
}

TEST(mp_inv_mod, inv_of_minus_one_is_minus_one_small_mod) {
    mp_uint_t mod{}, a{}, inv{}, one{}, prod{}, rem{};
    mp_set(mod, 17);
    mp_set(a, 16); // -1 mod 17
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    // (-1)*(-1)=1 mod 17
    mp_mul(prod, inv, 16);
    reduce_mod(rem, prod, mod);
    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(rem, one), 0);
}
