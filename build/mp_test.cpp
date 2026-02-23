#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <array>
#include "gtest/gtest.h"
#include "mp.hpp"
#include "fq.hpp"

#include <chrono>
#include <cstdio>

namespace {

static inline void set_limbs(uint64_t r[MP_N64], uint64_t l0, uint64_t l1, uint64_t l2, uint64_t l3) {
    r[0] = l0; r[1] = l1; r[2] = l2; r[3] = l3;
}

static inline void expect_eq(const uint64_t a[MP_N64], const uint64_t b[MP_N64]) {
    EXPECT_EQ(mp_cmp(a, b), 0);
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
        mul64_no128(a[i], k, hi, lo);

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

TEST(mp_copy, mp_copy) {
    mp_uint_t a{}, b{};
    set_limbs(a, 1, 2, 3, 4);
    mp_copy(b, a);
    expect_eq(a, b);
}

TEST(mp_cmp, mp_cmp) {
    mp_uint_t a{}, b{};
    mp_set(a, 0);
    mp_set(b, 0);
    EXPECT_EQ(mp_cmp(a, b), 0);

    set_limbs(a, 1, 0, 0, 0);
    set_limbs(b, 2, 0, 0, 0);
    EXPECT_LT(mp_cmp(a, b), 0);
    EXPECT_GT(mp_cmp(b, a), 0);

    // different high limb
    set_limbs(a, 0, 0, 0, 1);
    set_limbs(b, 0, 0, 0, 2);
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
    const std::array<std::array<uint64_t, MP_N64>, 6> cases_b = {{
        {{0,0,0,0}},
        {{2,0,0,0}},
        {{1,0,0,0}},
        {{1,0,0,0}},
        {{0,0,0,0x8000000000000000ULL}},
        {{0x1111111111111111ULL,0x2222222222222222ULL,0,0}},
    }};

    for (size_t i = 0; i < cases_a.size(); i++) {
        mp_uint_t r{}, ref1{}, ref2{};
        uint64_t c = mp_add(r, cases_a[i].data(), cases_b[i].data());
        uint64_t c1 = ref_add_no128(ref1, cases_a[i].data(), cases_b[i].data());
        EXPECT_EQ(c, c1);
        expect_eq(r, ref1);

#if defined(__SIZEOF_INT128__)
        uint64_t c2 = ref_add_128(ref2, cases_a[i].data(), cases_b[i].data());
        EXPECT_EQ(c, c2);
        expect_eq(r, ref2);
#endif
    }

    // boundary: max + 1 => 0 carry 1
    mp_uint_t max{}, one{}, out{};
    set_limbs(max, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
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
}

TEST(mp_sub_256, mp_sub) {
    mp_uint_t a{}, b{}, r{};

    // 0 - 1 => borrow 1, result = 2^256-1
    mp_set(a, 0);
    mp_set(b, 1);
    EXPECT_EQ(mp_sub(r, a, b), 1ULL);
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);

    // max - max = 0
    set_limbs(a, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
    mp_copy(b, a);
    EXPECT_EQ(mp_sub(r, a, b), 0ULL);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_add_u64, mp_add) {
    mp_uint_t a{}, r{};
    set_limbs(a, 0xFFFFFFFFFFFFFFFFULL, 0, 0, 0);
    EXPECT_EQ(mp_add(r, a, 1ULL), 0ULL);
    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);

    // overflow
    set_limbs(a, 0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(mp_add(r, a, 1ULL), 1ULL);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_sub_u64, mp_sub) {
    mp_uint_t a{}, r{};
    mp_set(a, 0);
    EXPECT_EQ(mp_sub(r, a, 1ULL), 1ULL);
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);

    set_limbs(a, 0, 1, 0, 0);
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
    const std::array<uint64_t, 5> cases_b = {{0ULL, 2ULL, 3ULL, 5ULL, 0xFEDCBA9876543211ULL}};

    for (size_t i = 0; i < cases_a.size(); i++) {
        mp_uint_t r{}, ref_no128{}, ref_128{};
        uint64_t c = mp_mul(r, cases_a[i].data(), cases_b[i]);

        uint64_t c_no128 = ref_mul_no128(ref_no128, cases_a[i].data(), cases_b[i]);
        EXPECT_EQ(c, c_no128);
        expect_eq(r, ref_no128);

#if defined(__SIZEOF_INT128__)
        uint64_t c_128 = ref_mul_128(ref_128, cases_a[i].data(), cases_b[i]);
        EXPECT_EQ(c, c_128);
        expect_eq(r, ref_128);
#endif
    }
}

TEST(mp_addmul_u64, mp_addmul) {
    // r += a*b
    mp_uint_t a{}, r{}, ref{};
    set_limbs(a, 3, 0, 0, 0);
    set_limbs(r, 5, 0, 0, 0);

    // expected: 5 + 3*7 = 26
    mp_copy(ref, r);
    mp_uint_t prod{};
    mp_mul(prod, a, 7);
    mp_add(ref, ref, prod);

    mp_addmul(r, a, 7);
    expect_eq(r, ref);
}

TEST(mp_and, mp_and) {
    mp_uint_t a{}, b{}, r{};
    set_limbs(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set_limbs(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_and(r, a, b);
    EXPECT_EQ(r[0], a[0] & b[0]);
    EXPECT_EQ(r[1], a[1] & b[1]);
    EXPECT_EQ(r[2], a[2] & b[2]);
    EXPECT_EQ(r[3], a[3] & b[3]);
}

TEST(mp_or, mp_or) {
    mp_uint_t a{}, b{}, r{};
    set_limbs(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set_limbs(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_or(r, a, b);
    EXPECT_EQ(r[0], a[0] | b[0]);
    EXPECT_EQ(r[1], a[1] | b[1]);
    EXPECT_EQ(r[2], a[2] | b[2]);
    EXPECT_EQ(r[3], a[3] | b[3]);
}

TEST(mp_xor, mp_xor) {
    mp_uint_t a{}, b{}, r{};
    set_limbs(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
    set_limbs(b, 0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL);
    mp_xor(r, a, b);
    EXPECT_EQ(r[0], a[0] ^ b[0]);
    EXPECT_EQ(r[1], a[1] ^ b[1]);
    EXPECT_EQ(r[2], a[2] ^ b[2]);
    EXPECT_EQ(r[3], a[3] ^ b[3]);
}

TEST(mp_not, mp_not) {
    mp_uint_t a{}, r{};
    set_limbs(a, 0xF0F0ULL, 0xAAAAULL, 0, ~0ULL);
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

    // out of range
    EXPECT_FALSE(mp_tstbit(a, 256));
}

TEST(mp_shl, mp_shl) {
    mp_uint_t a{}, r{};
    set_limbs(a, 1, 0, 0, 0);
    mp_shl(r, a, 1);
    EXPECT_EQ(r[0], 2ULL);

    // shift across limbs
    set_limbs(a, 0, 1, 0, 0);
    mp_shl(r, a, 64);
    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 1ULL);
}

TEST(mp_shr, mp_shr) {
    mp_uint_t a{}, r{};
    set_limbs(a, 2, 0, 0, 0);
    mp_shr(r, a, 1);
    EXPECT_EQ(r[0], 1ULL);

    // shift across limbs
    set_limbs(a, 0, 0, 1, 0);
    mp_shr(r, a, 64);
    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
}

TEST(mp_get_int32, mp_get_int32) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_EQ(mp_get_int32(a), 0);

    mp_set(a, 123);
    EXPECT_EQ(mp_get_int32(a), 123);

    set_limbs(a, 0xFFFFFFFFULL, 1, 0, 0);
    (void)mp_get_int32(a);
}

TEST(mp_fits_int32, mp_fits_int32) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX + 1ULL);
    EXPECT_FALSE(mp_fits_int32(a));

    set_limbs(a, 1, 1, 0, 0);
    EXPECT_FALSE(mp_fits_int32(a));
}

TEST(mp_set_str, mp_set) {
    mp_uint_t a{};
    ASSERT_TRUE(mp_set(a, "0", 10));
    EXPECT_TRUE(mp_is_zero(a));

    ASSERT_TRUE(mp_set(a, "1", 16));
    EXPECT_EQ(a[0], 1ULL);

    // max 256-bit value in hex
    ASSERT_TRUE(mp_set(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16));
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(a[i], 0xFFFFFFFFFFFFFFFFULL);

    // base 10 boundary: 2^64-1
    ASSERT_TRUE(mp_set(a, "18446744073709551615", 10));
    EXPECT_EQ(a[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(a[1], 0ULL);
}

TEST(mp_get_str, mp_get_str) {
    mp_uint_t a{};
    mp_set(a, 0);
    EXPECT_EQ(mp_get_str(a, 10), std::string("0"));

    mp_set(a, 255);
    EXPECT_EQ(mp_get_str(a, 16), std::string("ff"));
    EXPECT_EQ(mp_get_str(a, 2), std::string("11111111"));
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
    expect_eq(r, tmp);

    // negative: r = mod - 5
    mp_set_mod(r, -5, mod);
    mp_set(tmp, 5);
    reduce_mod(rem, tmp, mod);
    mp_sub(tmp, mod, rem);
    expect_eq(r, tmp);

    // large positive: (mod + 1) mod mod = 1
    // Construct x = mod + 1 (mod is 256-bit) and pass as int64? cannot.
    // Instead test int64 boundary: INT64_MIN/INT64_MAX.
    mp_set_mod(r, INT64_MAX, mod);
    mp_set(tmp, INT64_MAX);
    reduce_mod(tmp, tmp, mod);
    expect_eq(r, tmp);

    (void)q; (void)rem;
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
    set_limbs(max, ~0ULL, ~0ULL, ~0ULL, ~0ULL);
    reduce_mod(ref, max, mod);
    expect_eq(r, ref);

    // negative
    ASSERT_TRUE(mp_set_mod(r, "-1", 10, mod));
    mp_uint_t one{};
    mp_set(one, 1);
    mp_sub(ref, mod, one);
    expect_eq(r, ref);
}

TEST(mp_export_be, mp_export_be) {
    mp_uint_t a{};
    set_limbs(a, 0x1122334455667788ULL, 0x99AABBCCDDEEFF00ULL, 0x0102030405060708ULL, 0x1112131415161718ULL);

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

TEST(mp_div, mp_div) {
#if !defined(__SIZEOF_INT128__)
    GTEST_SKIP() << "mp_div requires __int128 in this build";
#endif

    // Construct num = den*q + r with q small (uint64), r < den.
    mp_uint_t den{}, qv{}, rv{}, num{}, tmp{}, q{}, r{};

    // den is 128-bit-ish
    set_limbs(den, 0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0);

    // q = 17
    mp_set(qv, 17);

    // tmp = den * 17
    mp_mul(tmp, den, 17);

    // r = 5
    mp_set(rv, 5);

    // num = tmp + 5
    mp_add(num, tmp, rv);

    mp_div(q, r, num, den);

    // expect quotient=17 remainder=5
    mp_uint_t qexp{}, rexp{};
    mp_set(qexp, 17);
    mp_set(rexp, 5);
    expect_eq(q, qexp);
    expect_eq(r, rexp);

    // boundary: num < den => q=0 r=num
    mp_set(num, 123);
    mp_div(q, r, num, den);
    mp_set(qexp, 0);
    mp_set(rexp, 123);
    expect_eq(q, qexp);
    expect_eq(r, rexp);
}

TEST(mp_pow_mod, mp_pow_mod) {
    const uint64_t *mod = Fq_q.longVal;

    mp_uint_t base{}, exp{}, out{}, ref{};

    // base=5 exp=0 => 1
    mp_set(base, 5);
    mp_set(exp, 0);
    mp_pow_mod(out, base, exp, mod);
    mp_set(ref, 1);
    reduce_mod(ref, ref, mod);
    expect_eq(out, ref);

    // base=5 exp=1 => 5
    mp_set(exp, 1);
    mp_pow_mod(out, base, exp, mod);
    mp_set(ref, 5);
    reduce_mod(ref, ref, mod);
    expect_eq(out, ref);

    // base=0 exp=0 => 1 (as implemented in pow loop)
    mp_set(base, 0);
    mp_set(exp, 0);
    mp_pow_mod(out, base, exp, mod);
    mp_set(ref, 1);
    reduce_mod(ref, ref, mod);
    expect_eq(out, ref);
}

TEST(mp_inv_mod, mp_inv_mod) {
    const uint64_t *mod = Fq_q.longVal;

    // a=0 not invertible
    mp_uint_t a{}, inv{};
    mp_set(a, 0);
    EXPECT_FALSE(mp_inv_mod(inv, a, mod));

    // a=1 invertible, inv=1
    mp_set(a, 1);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    mp_uint_t one{};
    mp_set(one, 1);
    expect_eq(inv, one);

    // a=2: check (a*inv) % mod == 1 using mp_mul + mp_div
    mp_set(a, 2);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));

    mp_uint_t prod{}, rem{};
    mp_mul(prod, inv, 2);
    reduce_mod(rem, prod, mod);
    expect_eq(rem, one);
}

static inline uint64_t now_ns() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
}

TEST(mp_perf_1m, addmul_and_mul_timing) {
    constexpr size_t N = 1000000;

    mp_uint_t a{}, b{}, r{};
    volatile uint64_t sink = 0;

    auto bench_addmul = [&](auto &&fn, const char *name) {
        uint64_t t0 = now_ns();
        for (size_t i = 0; i < N; i++) {
            set_limbs(a,
                      0x0123456789ABCDEFULL + i,
                      0x0FEDCBA987654321ULL ^ (i * 23),
                      0xAAAAAAAAAAAAAAAAULL + (i * 29),
                      0x5555555555555555ULL ^ (i * 31));

            set_limbs(r,
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
            set_limbs(a,
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

    RawFq::Element a, inv, prod;

    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) {
        Fq.inv(inv, a);
    }

    {
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