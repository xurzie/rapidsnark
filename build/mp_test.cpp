#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <array>
#include <climits>
#include "gtest/gtest.h"
#include "mp.hpp"


static constexpr mp_uint_t Fq_q = {0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL, 0xb85045b68181585dULL, 0x30644e72e131a029ULL};

TEST(mp_set_u64, mp_set) {
    mp_uint_t a;
    mp_set(a, 0);
    EXPECT_TRUE(mp_is_zero(a));

    mp_set(a, 123456789ULL);
    EXPECT_EQ(a[0], 123456789ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0ULL);
}

TEST(mp_copy, self_noop) {
    mp_uint_t a{1,2,3,4};
    mp_copy(a, a); // must be no-op
    EXPECT_EQ(a[0], 1ULL);
    EXPECT_EQ(a[1], 2ULL);
    EXPECT_EQ(a[2], 3ULL);
    EXPECT_EQ(a[3], 4ULL);
}

TEST(mp_copy, non_overlapping) {
    mp_uint_t a{1,2,3,4}, b;
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
    mp_uint_t a, b;
    mp_set(a, 0);
    mp_set(b, 0);
    EXPECT_EQ(mp_cmp(a, b), 0);

    {
        mp_uint_t aa = {1,0,0,0};
        mp_uint_t bb = {2,0,0,0};
        EXPECT_LT(mp_cmp(aa, bb), 0);
        EXPECT_GT(mp_cmp(bb, aa), 0);
    }

    // different high limb
    {
        mp_uint_t aa = {0,0,0,1};
        mp_uint_t bb = {0,0,0,2};
        EXPECT_LT(mp_cmp(aa, bb), 0);
    }
}

TEST(mp_is_zero, mp_is_zero) {
    mp_uint_t a;
    mp_set(a, 0);
    EXPECT_TRUE(mp_is_zero(a));
    a[2] = 1;
    EXPECT_FALSE(mp_is_zero(a));
}

TEST(mp_add_256, mp_add) {
    // Case 0: 0 + 0
    {
        mp_uint_t a = {0,0,0,0};
        mp_uint_t b = {0,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 0ULL);
        EXPECT_TRUE(mp_is_zero(r));
    }

    // Case 1: 1 + 2 = 3
    {
        mp_uint_t a = {1,0,0,0};
        mp_uint_t b = {2,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 0ULL);
        EXPECT_EQ(r[0], 3ULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // Case 2: (2^64-1) + 1 => 0 carry into limb1
    {
        mp_uint_t a = {0xFFFFFFFFFFFFFFFFULL,0,0,0};
        mp_uint_t b = {1,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 0ULL);
        EXPECT_EQ(r[0], 0ULL);
        EXPECT_EQ(r[1], 1ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // Case 3: (2^128-1) + 1 => limb0/1 zero, carry into limb2
    {
        mp_uint_t a = {0xFFFFFFFFFFFFFFFFULL,0xFFFFFFFFFFFFFFFFULL,0,0};
        mp_uint_t b = {1,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 0ULL);
        EXPECT_EQ(r[0], 0ULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 1ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // Case 4: top half addition without carry out
    {
        mp_uint_t a = {0,0,0,0x8000000000000000ULL};
        mp_uint_t b = {0,0,0,0x8000000000000000ULL};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 1ULL); // overflow beyond 2^256
        EXPECT_TRUE(mp_is_zero(r));
    }

    // Case 5: mixed limbs
    {
        mp_uint_t a = {0x0123456789ABCDEFULL,0x0ULL,0xFFFFFFFFFFFFFFFFULL,0x7ULL};
        mp_uint_t b = {0x1111111111111111ULL,0x2222222222222222ULL,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, b), 0ULL);
        EXPECT_EQ(r[0], 0x123456789ABCDF00ULL);
        EXPECT_EQ(r[1], 0x2222222222222222ULL);
        EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[3], 0x7ULL);
    }

    // boundary: max + 1 => 0 carry 1
    {
        mp_uint_t max = {~0ULL,~0ULL,~0ULL,~0ULL};
        mp_uint_t one, out;
        mp_set(one, 1);
        EXPECT_EQ(mp_add(out, max, one), 1ULL);
        EXPECT_TRUE(mp_is_zero(out));
    }
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
    mp_uint_t a, b, r;

    // 0 - 1 => borrow 1, result = 2^256-1
    mp_set(a, 0);
    mp_set(b, 1);
    EXPECT_EQ(mp_sub(r, a, b), 1ULL);
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);

    // max - max = 0
    {
        mp_uint_t aa = {~0ULL,~0ULL,~0ULL,~0ULL};
        mp_uint_t bb; mp_copy(bb, aa);
        EXPECT_EQ(mp_sub(r, aa, bb), 0ULL);
        EXPECT_TRUE(mp_is_zero(r));
    }

    // borrow ripple across all limbs: [0,0,0,1] - 1 = [~0,~0,~0,0]
    {
        mp_uint_t aa = {0,0,0,1};
        mp_uint_t bb; mp_set(bb, 1);
        EXPECT_EQ(mp_sub(r, aa, bb), 0ULL);
        EXPECT_EQ(r[0], ~0ULL);
        EXPECT_EQ(r[1], ~0ULL);
        EXPECT_EQ(r[2], ~0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }
}

TEST(mp_add_u64, mp_add) {
    {
        mp_uint_t a = {~0ULL, 0, 0, 0};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, 1ULL), 0ULL);
        EXPECT_EQ(r[0], 0ULL);
        EXPECT_EQ(r[1], 1ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // overflow
    {
        mp_uint_t a = {~0ULL,~0ULL,~0ULL,~0ULL};
        mp_uint_t r;
        EXPECT_EQ(mp_add(r, a, 1ULL), 1ULL);
        EXPECT_TRUE(mp_is_zero(r));
    }
}

TEST(mp_sub_u64, mp_sub) {
    {
        mp_uint_t a, r;
        mp_set(a, 0);
        EXPECT_EQ(mp_sub(r, a, 1ULL), 1ULL);
        for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(r[i], 0xFFFFFFFFFFFFFFFFULL);
    }

    {
        mp_uint_t a = {0,1,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_sub(r, a, 1ULL), 0ULL);
        EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }
}

TEST(mp_mul_u64, mp_mul) {
    // 0 * x
    {
        mp_uint_t a = {0,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_mul(r, a, 0xDEADBEEF), 0ULL);
        EXPECT_TRUE(mp_is_zero(r));
    }

    // 1 * 2
    {
        mp_uint_t a = {1,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_mul(r, a, 2), 0ULL);
        EXPECT_EQ(r[0], 2ULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // (2^64-1) * 3 => low = 0xFFFFFFFFFFFFFFFD, carry into limb1 = 2
    {
        mp_uint_t a = {0xFFFFFFFFFFFFFFFFULL,0,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_mul(r, a, 3), 0ULL);
        EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFDULL);
        EXPECT_EQ(r[1], 0x2ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // limb1 all ones * 5 (acts like shifting product into limb1..)
    {
        mp_uint_t a = {0,0xFFFFFFFFFFFFFFFFULL,0,0};
        mp_uint_t r;
        EXPECT_EQ(mp_mul(r, a, 5), 0ULL);
        EXPECT_EQ(r[0], 0ULL);
        EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFBULL);
        EXPECT_EQ(r[2], 0x4ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // boundary: (2^256-1) * 2 = (2^256-2) with carry=1
    {
        mp_uint_t a = {~0ULL,~0ULL,~0ULL,~0ULL};
        mp_uint_t r;
        EXPECT_EQ(mp_mul(r, a, 2), 1ULL);
        EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFEULL);
        EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[3], 0xFFFFFFFFFFFFFFFFULL);
    }
}
/*
TEST(mp_addmul_u64, mp_addmul) {
    // r += a*b : 5 + 3*7 = 26
    {
        mp_uint_t a = {3,0,0,0};
        mp_uint_t r = {5,0,0,0};

        const uint64_t c = mp_addmul(r, a, 7);
        EXPECT_EQ(c, 0ULL);
        EXPECT_EQ(r[0], 26ULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // overflow-heavy: r=max, a=max, b=2
    // a*2 (mod 2^256) = 2^256-2, plus r=max => 2^256-3, carry=2
    {
        mp_uint_t a = {~0ULL,~0ULL,~0ULL,~0ULL};
        mp_uint_t r = {~0ULL,~0ULL,~0ULL,~0ULL};

        const uint64_t c = mp_addmul(r, a, 2);
        EXPECT_EQ(c, 2ULL);
        EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFDULL);
        EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
        EXPECT_EQ(r[3], 0xFFFFFFFFFFFFFFFFULL);
    }
}
*/
TEST(mp_and, mp_and) {
    mp_uint_t a = {0xF0F0ULL, 0xAAAAULL, 0, ~0ULL};
    mp_uint_t b = {0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL};
    mp_uint_t r;

    mp_and(r, a, b);
    EXPECT_EQ(r[0], a[0] & b[0]);
    EXPECT_EQ(r[1], a[1] & b[1]);
    EXPECT_EQ(r[2], a[2] & b[2]);
    EXPECT_EQ(r[3], a[3] & b[3]);
}

TEST(mp_or, mp_or) {
    mp_uint_t a = {0xF0F0ULL, 0xAAAAULL, 0, ~0ULL};
    mp_uint_t b = {0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL};
    mp_uint_t r;

    mp_or(r, a, b);
    EXPECT_EQ(r[0], a[0] | b[0]);
    EXPECT_EQ(r[1], a[1] | b[1]);
    EXPECT_EQ(r[2], a[2] | b[2]);
    EXPECT_EQ(r[3], a[3] | b[3]);
}

TEST(mp_xor, mp_xor) {
    mp_uint_t a = {0xF0F0ULL, 0xAAAAULL, 0, ~0ULL};
    mp_uint_t b = {0x0FF0ULL, 0x0F0FULL, ~0ULL, 0x1234ULL};
    mp_uint_t r;

    mp_xor(r, a, b);
    EXPECT_EQ(r[0], a[0] ^ b[0]);
    EXPECT_EQ(r[1], a[1] ^ b[1]);
    EXPECT_EQ(r[2], a[2] ^ b[2]);
    EXPECT_EQ(r[3], a[3] ^ b[3]);
}

TEST(mp_not, mp_not) {
    mp_uint_t a = {0xF0F0ULL, 0xAAAAULL, 0, ~0ULL};
    mp_uint_t r;

    mp_not(r, a);
    EXPECT_EQ(r[0], ~a[0]);
    EXPECT_EQ(r[1], ~a[1]);
    EXPECT_EQ(r[2], ~a[2]);
    EXPECT_EQ(r[3], ~a[3]);
}

TEST(mp_tstbit, mp_tstbit) {
    mp_uint_t a;
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
    mp_uint_t a = {0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL};
    mp_uint_t r;
    mp_shl(r, a, 0);
    EXPECT_EQ(mp_cmp(r, a), 0);
}

TEST(mp_shl, zero_stays_zero) {
    mp_uint_t a, r;
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
    mp_uint_t a = {1,0,0,0};
    mp_uint_t r;
    mp_shl(r, a, 1);

    EXPECT_EQ(r[0], 2ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_63) {
    mp_uint_t a = {3ULL, 0ULL, 0ULL, 0ULL};
    mp_uint_t r;
    mp_shl(r, a, 63);

    EXPECT_EQ(r[0], 0x8000000000000000ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_64_bitShift0_path) {
    mp_uint_t a = {1ULL, 0ULL, 0ULL, 0ULL};
    mp_uint_t r;
    mp_shl(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, limb_boundary_65) {
    mp_uint_t a = {1ULL, 0ULL, 0ULL, 0ULL};
    mp_uint_t r;
    mp_shl(r, a, 65);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 2ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, shift_across_limbs_example) {
    mp_uint_t a = {0ULL, 1ULL, 0ULL, 0ULL};
    mp_uint_t r;
    mp_shl(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 1ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shl, shift_255_bit0_to_top_bit) {
    mp_uint_t a = {1ULL, 0ULL, 0ULL, 0ULL};
    mp_uint_t r;
    mp_shl(r, a, 255);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0x8000000000000000ULL);
}

TEST(mp_shl, shift_ge_256_is_zero) {
    mp_uint_t a = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
    mp_uint_t r;

    mp_shl(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shl(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shl, all_ones_signature_shift1) {
    mp_uint_t a = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
    mp_uint_t r;
    mp_shl(r, a, 1);

    EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFEULL);
    EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[3], 0xFFFFFFFFFFFFFFFFULL);
}

// -------------------- mp_shr --------------------

TEST(mp_shr, identity_shift0) {
    mp_uint_t a = {0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL};
    mp_uint_t r;
    mp_shr(r, a, 0);
    EXPECT_EQ(mp_cmp(r, a), 0);
}

TEST(mp_shr, zero_stays_zero) {
    mp_uint_t a, r;
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
    mp_uint_t a = {2,0,0,0};
    mp_uint_t r;
    mp_shr(r, a, 1);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_63) {
    mp_uint_t a = {0x8000000000000000ULL, 0x1ULL, 0x0ULL, 0x0ULL};
    mp_uint_t r;
    mp_shr(r, a, 63);

    EXPECT_EQ(r[0], 3ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_64_bitShift0_path) {
    mp_uint_t a = {0x0ULL, 0x1ULL, 0x0ULL, 0x0ULL};
    mp_uint_t r;
    mp_shr(r, a, 64);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, limb_boundary_65) {
    mp_uint_t a = {0x0ULL, 0x1ULL, 0x0ULL, 0x0ULL};
    mp_uint_t r;
    mp_shr(r, a, 65);

    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shr, shift_across_limbs_example) {
    mp_uint_t a = {0, 0, 1, 0};
    mp_uint_t r;
    mp_shr(r, a, 64);

    EXPECT_EQ(r[0], 0ULL);
    EXPECT_EQ(r[1], 1ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, shift_255_top_bit_to_bit0) {
    mp_uint_t a = {0x0ULL, 0x0ULL, 0x0ULL, 0x8000000000000000ULL};
    mp_uint_t r;
    mp_shr(r, a, 255);

    EXPECT_EQ(r[0], 1ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_shr, shift_ge_256_is_zero) {
    mp_uint_t a = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
    mp_uint_t r;

    mp_shr(r, a, 256);
    EXPECT_TRUE(mp_is_zero(r));

    mp_shr(r, a, 257);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_shr, all_ones_signature_shift1) {
    mp_uint_t a = {~0ULL, ~0ULL, ~0ULL, ~0ULL};
    mp_uint_t r;
    mp_shr(r, a, 1);

    EXPECT_EQ(r[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[1], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[2], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(r[3], 0x7FFFFFFFFFFFFFFFULL);
}

TEST(mp_get_int32, mp_get_int32) {
    mp_uint_t a;
    mp_set(a, 0);
    EXPECT_EQ(mp_get_int32(a), 0);

    mp_set(a, 123);
    EXPECT_EQ(mp_get_int32(a), 123);

    {
        mp_uint_t aa = {0xFFFFFFFFULL, 1, 0, 0};
        EXPECT_EQ(mp_get_int32(aa), -1);
    }

    {
        mp_uint_t aa = {0x80000000ULL, 0, 0, 0};
        EXPECT_EQ(mp_get_int32(aa), (int32_t)0x80000000u);
    }
}

TEST(mp_fits_int32, mp_fits_int32) {
    mp_uint_t a;
    mp_set(a, 0);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX);
    EXPECT_TRUE(mp_fits_int32(a));

    mp_set(a, (uint64_t)INT32_MAX + 1ULL);
    EXPECT_FALSE(mp_fits_int32(a));

    {
        mp_uint_t aa = {1,1,0,0};
        EXPECT_FALSE(mp_fits_int32(aa));
    }
}

TEST(mp_set_str, valid_cases) {
    mp_uint_t a;

    ASSERT_TRUE(mp_set(a, "0", 10));
    EXPECT_TRUE(mp_is_zero(a));

    ASSERT_TRUE(mp_set(a, "1", 16));
    EXPECT_EQ(a[0], 1ULL);

    ASSERT_TRUE(mp_set(a, "   +ff", 16));
    EXPECT_EQ(a[0], 255ULL);

    ASSERT_TRUE(mp_set(a, "42", 0));
    EXPECT_EQ(a[0], 42ULL);

    ASSERT_TRUE(mp_set(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16));
    for (size_t i = 0; i < MP_N64; i++) EXPECT_EQ(a[i], 0xFFFFFFFFFFFFFFFFULL);

    ASSERT_TRUE(mp_set(a, "18446744073709551615", 10));
    EXPECT_EQ(a[0], 0xFFFFFFFFFFFFFFFFULL);
    EXPECT_EQ(a[1], 0ULL);
}

TEST(mp_set_str, invalid_cases) {
    mp_uint_t a;

    EXPECT_FALSE(mp_set(a, "1", 1));
    EXPECT_FALSE(mp_set(a, "1", 17));

    EXPECT_FALSE(mp_set(a, "-1", 10));

    EXPECT_FALSE(mp_set(a, "2", 2));
    EXPECT_FALSE(mp_set(a, "g", 16));

    EXPECT_FALSE(mp_set(a, "", 10));
    EXPECT_FALSE(mp_set(a, "   ", 10));
    EXPECT_FALSE(mp_set(a, "+", 10));

    EXPECT_FALSE(mp_set(a, "123x", 10));
    EXPECT_FALSE(mp_set(a, "ff zz", 16));

    // overflow: 2^256 is 1 followed by 64 hex zeros => must fail
    EXPECT_FALSE(mp_set(a, "10000000000000000000000000000000000000000000000000000000000000000", 16));
}

TEST(mp_get_str, basic_cases) {
    mp_uint_t a;
    mp_set(a, 0);
    EXPECT_EQ(mp_get_str(a, 10), std::string("0"));

    mp_set(a, 255);
    EXPECT_EQ(mp_get_str(a, 16), std::string("ff"));
    EXPECT_EQ(mp_get_str(a, 2), std::string("11111111"));

    {
        mp_uint_t aa = {~0ULL,~0ULL,~0ULL,~0ULL};
        EXPECT_EQ(mp_get_str(aa, 16), std::string(64, 'f'));
    }

    mp_uint_t b;
    ASSERT_TRUE(mp_set(b, "00000100", 16));
    EXPECT_EQ(mp_get_str(b, 16), std::string("100"));
}

TEST(mp_get_str, invalid_base_returns_empty) {
    mp_uint_t a;
    mp_set(a, 123);

    EXPECT_EQ(mp_get_str(a, 1), std::string());
    EXPECT_EQ(mp_get_str(a, 17), std::string());
}

TEST(mp_set_mod_i64, mp_set_mod_Fq_q_small_constants) {
    const uint64_t *mod = Fq_q;

    // 0
    {
        mp_uint_t r;
        mp_set_mod(r, 0, mod);
        EXPECT_TRUE(mp_is_zero(r));
    }

    // +5
    {
        mp_uint_t r;
        mp_set_mod(r, 5, mod);
        EXPECT_EQ(r[0], 5ULL);
        EXPECT_EQ(r[1], 0ULL);
        EXPECT_EQ(r[2], 0ULL);
        EXPECT_EQ(r[3], 0ULL);
    }

    // -1 => mod-1
    {
        mp_uint_t r;
        mp_uint_t q_minus_1 = {
            0x3c208c16d87cfd46ULL,
            0x97816a916871ca8dULL,
            0xb85045b68181585dULL,
            0x30644e72e131a029ULL
        };
        mp_set_mod(r, -1, mod);
        EXPECT_EQ(mp_cmp(r, q_minus_1), 0);
    }

    // -5 => mod-5
    {
        mp_uint_t q_minus_5 = {
            0x3c208c16d87cfd42ULL,
            0x97816a916871ca8dULL,
            0xb85045b68181585dULL,
            0x30644e72e131a029ULL
        };
        mp_uint_t r;
        mp_set_mod(r, -5, mod);
        EXPECT_EQ(mp_cmp(r, q_minus_5), 0);
    }
}

TEST(mp_set_mod_i64, small_mod_extremes) {
    mp_uint_t mod;
    mp_set(mod, 17);

    auto expect_i64 = [&](int64_t x, uint64_t expected) {
        mp_uint_t got;
        mp_set_mod(got, x, mod);
        EXPECT_EQ(got[0], expected);
        EXPECT_EQ(got[1], 0ULL);
        EXPECT_EQ(got[2], 0ULL);
        EXPECT_EQ(got[3], 0ULL);
    };

    expect_i64(0, 0);
    expect_i64(1, 1);
    expect_i64(-1, 16);
    expect_i64(17, 0);
    expect_i64(-17, 0);
    expect_i64(18, 1);
    expect_i64(-18, 16);

    // mod = 1 => always 0
    mp_uint_t mod1, r;
    mp_set(mod1, 1);
    mp_set_mod(r, INT64_MIN, mod1);
    EXPECT_TRUE(mp_is_zero(r));
    mp_set_mod(r, INT64_MAX, mod1);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_set_mod_str, mp_set_mod_Fq_q_small_constants) {
    const uint64_t *mod = Fq_q;

    mp_uint_t r;
    mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    ASSERT_TRUE(mp_set_mod(r, "0", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    ASSERT_TRUE(mp_set_mod(r, "5", 10, mod));
    EXPECT_EQ(r[0], 5ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);

    ASSERT_TRUE(mp_set_mod(r, "-1", 10, mod));
    EXPECT_EQ(mp_cmp(r, q_minus_1), 0);

    ASSERT_TRUE(mp_set_mod(r, "-0", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    // base=0 defaults to 10
    ASSERT_TRUE(mp_set_mod(r, "15", 0, mod));
    EXPECT_EQ(r[0], 15ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_set_mod_str, digit_must_be_lt_base_small_mod) {
    mp_uint_t mod, r;
    mp_set(mod, 17);

    EXPECT_FALSE(mp_set_mod(r, "2", 2, mod));
    EXPECT_FALSE(mp_set_mod(r, "10 1", 2, mod));

    EXPECT_FALSE(mp_set_mod(r, "1", 1, mod));
    EXPECT_FALSE(mp_set_mod(r, "1", 17, mod));

    EXPECT_FALSE(mp_set_mod(r, "", 10, mod));
    EXPECT_FALSE(mp_set_mod(r, "   ", 10, mod));
    EXPECT_FALSE(mp_set_mod(r, "+", 10, mod));

    // "34" mod 17 = 0
    ASSERT_TRUE(mp_set_mod(r, "34", 10, mod));
    EXPECT_TRUE(mp_is_zero(r));

    // "-1" mod 17 = 16
    ASSERT_TRUE(mp_set_mod(r, "-1", 10, mod));
    EXPECT_EQ(r[0], 16ULL);
    EXPECT_EQ(r[1], 0ULL);
    EXPECT_EQ(r[2], 0ULL);
    EXPECT_EQ(r[3], 0ULL);
}

TEST(mp_export_be, mp_export_be) {
    mp_uint_t a = {
        0x1122334455667788ULL,
        0x99AABBCCDDEEFF00ULL,
        0x0102030405060708ULL,
        0x1112131415161718ULL
    };

    uint8_t out[MP_N] = {0};
    mp_export_be(out, a);

    EXPECT_EQ(out[0], (uint8_t)0x11);
    EXPECT_EQ(out[1], (uint8_t)0x12);
    EXPECT_EQ(out[2], (uint8_t)0x13);
    EXPECT_EQ(out[3], (uint8_t)0x14);
    EXPECT_EQ(out[31], (uint8_t)0x88);
}

TEST(mp_import_be, mp_import_be) {
    uint8_t in[MP_N] = {0};
    in[31] = 0x01;

    mp_uint_t a;
    mp_import_be(a, in);
    EXPECT_EQ(a[0], 1ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0ULL);
}

TEST(mp_import_be, msb_sets_bit255) {
    uint8_t in[MP_N] = {0};
    in[0] = 0x80;
    mp_uint_t a;
    mp_import_be(a, in);
    EXPECT_EQ(a[0], 0ULL);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0x8000000000000000ULL);
}

TEST(mp_export_import_be, roundtrip) {
    mp_uint_t a = {
        0x1122334455667788ULL,
        0x99AABBCCDDEEFF00ULL,
        0x0102030405060708ULL,
        0x1112131415161718ULL
    };
    mp_uint_t b;

    uint8_t buf[MP_N] = {0};
    mp_export_be(buf, a);
    mp_import_be(b, buf);

    EXPECT_EQ(mp_cmp(a, b), 0);
}

TEST(mp_export_import_be, all_ones_roundtrip) {
    mp_uint_t a = {~0ULL,~0ULL,~0ULL,~0ULL};
    mp_uint_t b;

    uint8_t buf[MP_N] = {0};
    mp_export_be(buf, a);
    mp_import_be(b, buf);

    EXPECT_EQ(mp_cmp(a, b), 0);
}

TEST(mp_div, num_lt_den_q0_rnum) {
    mp_uint_t den = {0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0};
    mp_uint_t num, q, r;
    mp_set(num, 123);

    mp_div(q, r, num, den);

    mp_uint_t qexp, rexp;
    mp_set(qexp, 0);
    mp_set(rexp, 123);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, num_eq_den_q1_r0) {
    mp_uint_t den = {0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0};
    mp_uint_t num, q, r, one, zero;
    mp_copy(num, den);

    mp_div(q, r, num, den);

    mp_set(one, 1);
    mp_set(zero, 0);
    EXPECT_EQ(mp_cmp(q, one), 0);
    EXPECT_EQ(mp_cmp(r, zero), 0);
}

TEST(mp_div, den_1word_path) {
    mp_uint_t den, num, q, r, qexp, rexp;

    mp_set(den, 7);
    mp_set(num, 1000);

    mp_div(q, r, num, den);

    mp_set(qexp, 142);
    mp_set(rexp, 6);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, knuth_path_s_nonzero) {
    mp_uint_t den = {0xFFFFFFFFFFFFFFF1ULL, 0x123456789ABCDEF0ULL, 0, 0};
    mp_uint_t qv, rv, num, tmp, q, r;

    mp_set(qv, 17);
    mp_mul(tmp, den, 17);
    mp_set(rv, 5);
    mp_add(num, tmp, rv);

    mp_div(q, r, num, den);

    mp_uint_t qexp, rexp;
    mp_set(qexp, 17);
    mp_set(rexp, 5);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, knuth_path_s_zero) {
    mp_uint_t den = {0x0123456789ABCDEFULL, 0x8000000000000000ULL, 0, 0};
    mp_uint_t num, tmp, q, r;

    mp_mul(tmp, den, 9);
    mp_uint_t seven;
    mp_set(seven, 7);
    mp_add(num, tmp, seven);

    mp_div(q, r, num, den);

    mp_uint_t qexp, rexp;
    mp_set(qexp, 9);
    mp_set(rexp, 7);
    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, rexp), 0);
}

TEST(mp_div, den_one_returns_num) {
    mp_uint_t den, num, q, r;
    mp_set(den, 1);
    mp_uint_t n = {0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x8000000000000001ULL};
    mp_copy(num, n);

    mp_div(q, r, num, den);

    EXPECT_EQ(mp_cmp(q, num), 0);
    EXPECT_TRUE(mp_is_zero(r));
}

TEST(mp_div, top_bit_quotient_case) {
    mp_uint_t num = {0,0,0,0x8000000000000000ULL};
    mp_uint_t den, q, r, qexp, z;
    mp_set(den, 2);

    mp_div(q, r, num, den);

    mp_set(z, 0);
    mp_uint_t qq = {0,0,0,0x4000000000000000ULL};
    mp_copy(qexp, qq);

    EXPECT_EQ(mp_cmp(q, qexp), 0);
    EXPECT_EQ(mp_cmp(r, z), 0);
}

TEST(mp_pow_mod, exp_zero_returns_one_mod) {
    const uint64_t *mod = Fq_q;

    mp_uint_t base, exp, out;
    mp_set(base, 5);
    mp_set(exp, 0);

    mp_pow_mod(out, base, exp, mod);

    mp_uint_t one;
    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(out, one), 0);
}

TEST(mp_pow_mod, exp_one_returns_base_mod_small) {
    const uint64_t *mod = Fq_q;

    mp_uint_t base, exp, out;
    mp_set(base, 5);
    mp_set(exp, 1);

    mp_pow_mod(out, base, exp, mod);

    mp_uint_t five;
    mp_set(five, 5);
    EXPECT_EQ(mp_cmp(out, five), 0);
}

TEST(mp_pow_mod, base_zero_exp_zero_defined_as_one_mod) {
    const uint64_t *mod = Fq_q;

    mp_uint_t base, exp, out, one;
    mp_set(base, 0);
    mp_set(exp, 0);

    mp_pow_mod(out, base, exp, mod);

    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(out, one), 0);
}

TEST(mp_pow_mod, mod_one_returns_zero) {
    mp_uint_t base, exp, mod1, out;
    mp_set(base, 123);
    mp_set(exp, 456);
    mp_set(mod1, 1);

    mp_pow_mod(out, base, exp, mod1);
    EXPECT_TRUE(mp_is_zero(out));
}

TEST(mp_pow_mod, base_reduction_path_base_ge_mod_small_mod) {
    mp_uint_t mod, base, exp, out, ref;
    mp_set(mod, 17);

    mp_set(base, 22); // 17+5
    mp_set(exp, 3);

    mp_pow_mod(out, base, exp, mod);

    // 5^3 mod 17 = 6
    mp_set(ref, 6);
    EXPECT_EQ(mp_cmp(out, ref), 0);
}

TEST(mp_pow_mod, big_numbers_bn254_q_minus_1) {
    // Use BN254 scalar field order as mod, and base = mod-1.
    // (mod-1)^2 mod mod = 1, and (mod-1)^3 mod mod = mod-1.
    const uint64_t *mod = Fq_q;

    mp_uint_t base, exp, out, one;
    mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    // base = mod - 1 (constant)
    mp_copy(base, q_minus_1);

    // exponent = 2
    mp_set(exp, 2);
    mp_pow_mod(out, base, exp, mod);
    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(out, one), 0);

    // exponent = 3 => result = mod - 1
    mp_set(exp, 3);
    mp_pow_mod(out, base, exp, mod);
    EXPECT_EQ(mp_cmp(out, q_minus_1), 0);
}

TEST(mp_pow_mod, big_numbers_q_minus_5_pow2_pow3) {
    const uint64_t *mod = Fq_q;

    // q - 5
    const mp_uint_t q_minus_5 = {
        0x3c208c16d87cfd42ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    // q - 125
    const mp_uint_t q_minus_125 = {
        0x3c208c16d87cfccaULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    mp_uint_t base, exp, out, expect;

    mp_copy(base, q_minus_5);

    // (q-5)^2 mod q = 25
    mp_set(exp, 2);
    mp_pow_mod(out, base, exp, mod);
    mp_set(expect, 25);
    EXPECT_EQ(mp_cmp(out, expect), 0);

    // (q-5)^3 mod q = q-125
    mp_set(exp, 3);
    mp_pow_mod(out, base, exp, mod);
    EXPECT_EQ(mp_cmp(out, q_minus_125), 0);
}

TEST(mp_pow_mod, random_big_base_and_big_odd_exp_minus_one_result) {
    const uint64_t *mod = Fq_q;

    // q - 1  (то есть -1 mod q)
    const mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    // большой "рандомный" (но заданный константой) НЕЧЁТНЫЙ показатель
    const mp_uint_t exp_big_odd = {
        0xDEADBEEFCAFEBABFULL,  // LSB=1 => odd
        0x0123456789ABCDEFULL,
        0x0ULL,
        0x0ULL
    };

    mp_uint_t out;
    mp_pow_mod(out, q_minus_1, exp_big_odd, mod);

    // (-1)^odd == -1
    EXPECT_EQ(mp_cmp(out, q_minus_1), 0);
}

TEST(mp_inv_mod, zero_not_invertible) {
    const uint64_t *mod = Fq_q;
    mp_uint_t a, inv;
    mp_set(a, 0);
    EXPECT_FALSE(mp_inv_mod(inv, a, mod));
}

TEST(mp_inv_mod, one_inverts_to_one) {
    const uint64_t *mod = Fq_q;
    mp_uint_t a, inv, one;
    mp_set(a, 1);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    mp_set(one, 1);
    EXPECT_EQ(mp_cmp(inv, one), 0);
}

TEST(mp_inv_mod, small_mod_known_values) {
    // mod = 17:
    // inv(3) = 6 because 3*6 = 18 = 1 (mod 17)
    // inv(16) = 16 because (-1)^-1 = -1
    mp_uint_t mod, a, inv, ref;
    mp_set(mod, 17);

    mp_set(a, 3);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    mp_set(ref, 6);
    EXPECT_EQ(mp_cmp(inv, ref), 0);

    mp_set(a, 16);
    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    mp_set(ref, 16);
    EXPECT_EQ(mp_cmp(inv, ref), 0);
}

TEST(mp_inv_mod, big_numbers_bn254_q_minus_1) {
    // Inverse in BN254 scalar field:
    // inv(mod-1) == mod-1 (since (-1)^-1 = -1).
    const uint64_t *mod = Fq_q;
    mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    mp_uint_t a, inv;
    mp_copy(a, q_minus_1);

    ASSERT_TRUE(mp_inv_mod(inv, a, mod));
    EXPECT_EQ(mp_cmp(inv, q_minus_1), 0);
}

// ============================================================
// DIAGNOSTIC TESTS FOR MODULAR / ALIASING LAYER
// ============================================================

static inline uint64_t add_carry(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__clang__) || defined(__GNUC__)
    uint64_t t;
    unsigned c1 = __builtin_add_overflow(a, b, &t);
    unsigned c2 = __builtin_add_overflow(t, c, &t);
    *out = t;
    return (c1 | c2);
#else
    uint64_t t0 = a + b;
    uint64_t carry1 = (t0 < a);
    uint64_t t1 = t0 + c;
    uint64_t carry2 = (t1 < t0);
    *out = t1;
    return carry1 | carry2;
#endif
}

static inline uint64_t sub_borrow(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__clang__) || defined(__GNUC__)
    uint64_t t;
    unsigned b1 = __builtin_sub_overflow(a, b, &t);
    unsigned b2 = __builtin_sub_overflow(t, c, &t);
    *out = t;
    return (b1 | b2);
#else
    uint64_t t0 = a - b;
    uint64_t borrow1 = (a < b);
    uint64_t t1 = t0 - c;
    uint64_t borrow2 = (t0 < c);
    *out = t1;
    return borrow1 | borrow2;
#endif
}

static inline void mp_add_mod(uint64_t *r,
                              const uint64_t *a,
                              const uint64_t *b,
                              const uint64_t *mod)
{
    uint64_t carry = mp_add(r, a, b);

    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
}

static inline void mp_sub_mod(uint64_t *x, const uint64_t *y, const uint64_t *mod) {
    // x = x - y (mod mod), assumes 0 <= x,y < mod and mod odd
    uint64_t t0, t1, t2, t3;
    uint64_t br = 0;

    br = sub_borrow(&t0, x[0], y[0], br);
    br = sub_borrow(&t1, x[1], y[1], br);
    br = sub_borrow(&t2, x[2], y[2], br);
    br = sub_borrow(&t3, x[3], y[3], br);

    // if borrow, add mod back
    if (br) {
        uint64_t c = 0;
        c = add_carry(&t0, t0, mod[0], c);
        c = add_carry(&t1, t1, mod[1], c);
        c = add_carry(&t2, t2, mod[2], c);
        c = add_carry(&t3, t3, mod[3], c);
        (void)c;
    }

    x[0] = t0; x[1] = t1; x[2] = t2; x[3] = t3;
}

static inline void mp_div2_mod(uint64_t *x, const uint64_t *mod) {
    // if x is odd: x += mod
    if (x[0] & 1u) {
        uint64_t c = 0;
        c = add_carry(&x[0], x[0], mod[0], c);
        c = add_carry(&x[1], x[1], mod[1], c);
        c = add_carry(&x[2], x[2], mod[2], c);
        c = add_carry(&x[3], x[3], mod[3], c);
        (void)c;
    }
    // x >>= 1
    uint64_t b3 = x[3];
    uint64_t b2 = x[2];
    uint64_t b1 = x[1];
    uint64_t b0 = x[0];
    x[0] = (b0 >> 1) | (b1 << (2*MP_N - 1));
    x[1] = (b1 >> 1) | (b2 << (2*MP_N - 1));
    x[2] = (b2 >> 1) | (b3 << (2*MP_N - 1));
    x[3] = (b3 >> 1);
}

static void expect_mp_eq_u64(const mp_uint_t a, uint64_t x) {
    EXPECT_EQ(a[0], x);
    EXPECT_EQ(a[1], 0ULL);
    EXPECT_EQ(a[2], 0ULL);
    EXPECT_EQ(a[3], 0ULL);
}

TEST(mp_add_aliasing, r_eq_a) {
    mp_uint_t a = {5, 7, 11, 13};
    mp_uint_t b = {17, 19, 23, 29};

    mp_uint_t ref;
    mp_add(ref, a, b);

    mp_add(a, a, b);
    EXPECT_EQ(mp_cmp(a, ref), 0);
}

TEST(mp_add_aliasing, r_eq_b) {
    mp_uint_t a = {5, 7, 11, 13};
    mp_uint_t b = {17, 19, 23, 29};

    mp_uint_t ref;
    mp_add(ref, a, b);

    mp_add(b, a, b);
    EXPECT_EQ(mp_cmp(b, ref), 0);
}

TEST(mp_sub_aliasing, r_eq_a) {
    mp_uint_t a = {100, 7, 11, 13};
    mp_uint_t b = {17, 19, 23, 29};

    mp_uint_t ref;
    mp_sub(ref, a, b);

    mp_sub(a, a, b);
    EXPECT_EQ(mp_cmp(a, ref), 0);
}

TEST(mp_sub_aliasing, r_eq_b) {
    mp_uint_t a = {100, 7, 11, 13};
    mp_uint_t b = {17, 19, 23, 29};

    mp_uint_t ref;
    mp_sub(ref, a, b);

    // here result is different semantic: b := a-b
    mp_sub(b, a, b);
    EXPECT_EQ(mp_cmp(b, ref), 0);
}

TEST(mp_shl_aliasing, inplace) {
    mp_uint_t a = {0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x7ULL};
    mp_uint_t ref;
    mp_shl(ref, a, 17);
    mp_shl(a, a, 17);
    EXPECT_EQ(mp_cmp(a, ref), 0);
}

TEST(mp_shr_aliasing, inplace) {
    mp_uint_t a = {0x0123456789ABCDEFULL, 0x0ULL, 0xFFFFFFFFFFFFFFFFULL, 0x7ULL};
    mp_uint_t ref;
    mp_shr(ref, a, 17);
    mp_shr(a, a, 17);
    EXPECT_EQ(mp_cmp(a, ref), 0);
}

TEST(mp_add_mod, exhaustive_small_prime_17_normalized_inputs) {
    mp_uint_t mod;
    mp_set(mod, 17);

    for (uint64_t x = 0; x < 17; x++) {
        for (uint64_t y = 0; y < 17; y++) {
            mp_uint_t a, b, r;
            mp_set(a, x);
            mp_set(b, y);

            mp_add_mod(r, a, b, mod);

            uint64_t exp = (x + y) % 17;
            expect_mp_eq_u64(r, exp);
        }
    }
}

TEST(mp_sub_mod_inplace, exhaustive_small_prime_17_normalized_inputs) {
    mp_uint_t mod;
    mp_set(mod, 17);

    for (uint64_t x0 = 0; x0 < 17; x0++) {
        for (uint64_t y0 = 0; y0 < 17; y0++) {
            mp_uint_t x, y;
            mp_set(x, x0);
            mp_set(y, y0);

            mp_sub_mod(x, y, mod);

            uint64_t exp = (x0 + 17 - y0) % 17;
            EXPECT_EQ(x[0], exp) << "x0=" << x0 << " y0=" << y0;
            EXPECT_EQ(x[1], 0ULL);
            EXPECT_EQ(x[2], 0ULL);
            EXPECT_EQ(x[3], 0ULL);
        }
    }
}

TEST(mp_sub_mod_inplace, small_prime_17_non_normalized_x_should_fail_or_be_unsupported) {
    // Этот тест нужен только чтобы увидеть, не используешь ли ты mp_sub_mod
    // там, где x >= mod. По контракту у тебя assumes 0 <= x,y < mod.
    mp_uint_t mod;
    mp_set(mod, 17);

    mp_uint_t x, y;
    mp_set(x, 22); // 22 % 17 = 5, но функция НЕ обязана это нормализовать
    mp_set(y, 3);

    mp_sub_mod(x, y, mod);

    // Ничего не ASSERT-им как correctness.
    // Просто логика: если здесь вдруг получится 2, значит функция случайно
    // работает и на ненормализованных x; если нет — контракт подтверждается.
    SUCCEED();
}

TEST(mp_div2_mod, exhaustive_small_prime_17) {
    mp_uint_t mod;
    mp_set(mod, 17);

    for (uint64_t x0 = 0; x0 < 17; x0++) {
        mp_uint_t x;
        mp_set(x, x0);

        mp_div2_mod(x, mod);

        uint64_t exp = (x0 & 1ULL) ? ((x0 + 17) >> 1) : (x0 >> 1);
        expect_mp_eq_u64(x, exp);
    }
}

TEST(mp_div2_mod, zero) {
    mp_uint_t mod, x;
    mp_set(mod, 17);
    mp_set(x, 0);

    mp_div2_mod(x, mod);
    expect_mp_eq_u64(x, 0);
}

TEST(mp_div2_mod, bn254_edges) {
    const uint64_t *mod = Fq_q;

    const mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    mp_uint_t x, one;
    mp_set(one, 1);

    // 1 / 2 mod q = (q+1)/2
    mp_copy(x, one);
    mp_div2_mod(x, mod);

    mp_uint_t expect = {
        0x9e10460b6c3e7ea4ULL,
        0xcbc0b548b438e546ULL,
        0xdc2822db40c0ac2eULL,
        0x183227397098d014ULL
    };
    EXPECT_EQ(mp_cmp(x, expect), 0);

    // (q-1)/2 mod q = q/2 rounded down
    mp_copy(x, q_minus_1);
    mp_div2_mod(x, mod);

    mp_uint_t expect2 = {
        0x9e10460b6c3e7ea3ULL,
        0xcbc0b548b438e546ULL,
        0xdc2822db40c0ac2eULL,
        0x183227397098d014ULL
    };
    EXPECT_EQ(mp_cmp(x, expect2), 0);
}

TEST(mp_inv_mod, exhaustive_small_prime_17_table) {
    mp_uint_t mod;
    mp_set(mod, 17);

    const uint64_t inv17[17] = {
        0, 1, 9, 6, 13, 7, 3, 5, 15,
        2, 12, 14, 10, 4, 11, 8, 16
    };

    for (uint64_t x = 1; x < 17; x++) {
        mp_uint_t a, inv;
        mp_set(a, x);

        ASSERT_TRUE(mp_inv_mod(inv, a, mod)) << "x=" << x;
        expect_mp_eq_u64(inv, inv17[x]);
    }
}

TEST(mp_pow_mod, exhaustive_small_prime_17_against_naive) {
    mp_uint_t mod;
    mp_set(mod, 17);

    auto naive_pow_mod = [](uint64_t a, uint64_t e, uint64_t m) {
        uint64_t r = 1 % m;
        a %= m;
        while (e--) r = (r * a) % m;
        return r;
    };

    for (uint64_t a64 = 0; a64 < 17; a64++) {
        for (uint64_t e64 = 0; e64 <= 20; e64++) {
            mp_uint_t a, e, out;
            mp_set(a, a64);
            mp_set(e, e64);

            mp_pow_mod(out, a, e, mod);

            uint64_t exp = naive_pow_mod(a64, e64, 17);
            expect_mp_eq_u64(out, exp);
        }
    }
}

TEST(mp_add_mod, bn254_edges) {
    const uint64_t *mod = Fq_q;

    const mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    const mp_uint_t q_minus_2 = {
        0x3c208c16d87cfd45ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    mp_uint_t one, two, r;
    mp_set(one, 1);
    mp_set(two, 2);

    // (q-1) + 1 = 0 mod q
    mp_add_mod(r, q_minus_1, one, mod);
    EXPECT_TRUE(mp_is_zero(r));

    // (q-1) + 2 = 1 mod q
    mp_add_mod(r, q_minus_1, two, mod);
    expect_mp_eq_u64(r, 1);

    // (q-1) + (q-1) = q-2 mod q
    mp_add_mod(r, q_minus_1, q_minus_1, mod);
    EXPECT_EQ(mp_cmp(r, q_minus_2), 0);
}

TEST(mp_sub_mod_inplace, bn254_edges) {
    const uint64_t *mod = Fq_q;

    const mp_uint_t q_minus_1 = {
        0x3c208c16d87cfd46ULL,
        0x97816a916871ca8dULL,
        0xb85045b68181585dULL,
        0x30644e72e131a029ULL
    };

    mp_uint_t x, one, two, zero;
    mp_set(one, 1);
    mp_set(two, 2);
    mp_set(zero, 0);

    // 1 - 2 = q-1 mod q
    mp_copy(x, one);
    mp_sub_mod(x, two, mod);
    EXPECT_EQ(mp_cmp(x, q_minus_1), 0);

    // 1 - 1 = 0 mod q
    mp_copy(x, one);
    mp_sub_mod(x, one, mod);
    EXPECT_TRUE(mp_is_zero(x));

    // 0 - 1 = q-1 mod q
    mp_copy(x, zero);
    mp_sub_mod(x, one, mod);
    EXPECT_EQ(mp_cmp(x, q_minus_1), 0);
}

TEST(mp_add_mod, aliasing_r_eq_a_and_r_eq_b) {
    mp_uint_t mod;
    mp_set(mod, 17);

    // r == a
    {
        mp_uint_t a, b, ref;
        mp_set(a, 5);
        mp_set(b, 9);
        mp_set(ref, (5 + 9) % 17);

        mp_add_mod(a, a, b, mod);
        EXPECT_EQ(mp_cmp(a, ref), 0);
    }

    // r == b
    {
        mp_uint_t a, b, ref;
        mp_set(a, 5);
        mp_set(b, 9);
        mp_set(ref, (5 + 9) % 17);

        mp_add_mod(b, a, b, mod);
        EXPECT_EQ(mp_cmp(b, ref), 0);
    }
}

uint64_t mp_addmul(uint64_t *r, size_t rn, const uint64_t *a, size_t an, uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    size_t i = 0;

    for (; i < an && i < rn; ++i) {
        __uint128_t t = (__uint128_t)r[i] + (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }

    for (; i < rn; ++i) {
        __uint128_t t = (__uint128_t)r[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }

    return (uint64_t)carry;
#else
    uint64_t carry = 0;
    size_t i = 0;

    for (; i < an && i < rn; ++i) {
        uint64_t lo, hi;
        mp_mul64(a[i], b, lo, hi);

        uint64_t t1 = r[i] + lo;
        uint64_t c1 = (t1 < r[i]) ? 1ULL : 0ULL;

        uint64_t t2 = t1 + carry;
        uint64_t c2 = (t2 < t1) ? 1ULL : 0ULL;

        r[i] = t2;

        uint64_t new_carry = hi;
        new_carry += c1;
        new_carry += c2;
        carry = new_carry;
    }

    for (; i < rn; ++i) {
        uint64_t out = r[i] + carry;
        carry = (out < r[i]) ? 1ULL : 0ULL;
        r[i] = out;
    }

    return carry;
#endif
}

TEST(mp_addmul_varlen, carry_out_of_low_part_only)
{
    uint64_t r[5] = {
        0xffffffffffffffffULL,
        0xffffffffffffffffULL,
        0xffffffffffffffffULL,
        0xffffffffffffffffULL,
        7ULL
    };
    uint64_t a[4] = {
        0xffffffffffffffffULL,
        0xffffffffffffffffULL,
        0xffffffffffffffffULL,
        0xffffffffffffffffULL
    };

    uint64_t carry = mp_addmul(r, 4, a, 4, 2ULL);

    EXPECT_EQ(r[0], 0xfffffffffffffffdULL);
    EXPECT_EQ(r[1], 0xffffffffffffffffULL);
    EXPECT_EQ(r[2], 0xffffffffffffffffULL);
    EXPECT_EQ(r[3], 0xffffffffffffffffULL);
    EXPECT_EQ(r[4], 7ULL);
    EXPECT_EQ(carry, 2ULL);
}