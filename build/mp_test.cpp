
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <random>

#include "gtest/gtest.h"
#include "mp.hpp"

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

static cpp_int to_cpp_int(const uint64_t* a, size_t n) {
    cpp_int v = 0;
    for (size_t i = 0; i < n; i++) {
        v += (cpp_int(a[i]) << (64 * i));
    }
    return v;
}

static void from_cpp_int(uint64_t* r, size_t n, cpp_int v) {
    // clamp to n limbs (little-endian limbs)
    cpp_int mask = (cpp_int(1) << (64 * n)) - 1;
    v &= mask;
    for (size_t i = 0; i < n; i++) {
        // Boost's cpp_int uses expression templates; conversion to built-in types is explicit.
        // Extract each 64-bit limb via a concrete cpp_int and convert explicitly.
        cpp_int limb = (v >> (64 * i)) & cpp_int(0xFFFFFFFFFFFFFFFFULL);
        r[i] = limb.convert_to<uint64_t>();
    }
}

static void set_u64(uint64_t* r, uint64_t x) {
    mp_set(r, x);
}

static void set_i64_mod(uint64_t* r, int64_t x, const uint64_t* mod) {
    mp_set_mod(r, x, mod);
}

static std::vector<uint8_t> export_be(const uint64_t* a) {
    std::vector<uint8_t> out(MP_N64 * 8);
    mp_export_be(out.data(), a);
    return out;
}

static void import_be(uint64_t* r, const std::vector<uint8_t>& in) {
    ASSERT_EQ(in.size(), (size_t)MP_N64 * 8);
    mp_import_be(r, in.data());
}

static void expect_eq_mp(const uint64_t* a, const uint64_t* b) {
    EXPECT_EQ(mp_cmp(a, b), 0);
}

static cpp_int mod_pow(cpp_int base, cpp_int exp, const cpp_int& mod) {
    base %= mod;
    cpp_int res = 1 % mod;
    while (exp > 0) {
        if ((exp & 1) != 0) res = (res * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return res;
}

TEST(mp, set_copy_cmp_isZero) {
    mp_uint_t a{}, b{}, z{};
    set_u64(a, 0);
    EXPECT_TRUE(mp_is_zero(a));

    set_u64(a, 123456789ULL);
    EXPECT_FALSE(mp_is_zero(a));

    mp_copy(b, a);
    expect_eq_mp(a, b);

    set_u64(z, 0);
    EXPECT_GT(mp_cmp(a, z), 0);
    EXPECT_LT(mp_cmp(z, a), 0);
}

TEST(mp, shl_shr_roundtrip) {
    std::mt19937_64 rng(42);
    for (int t = 0; t < 200; t++) {
        mp_uint_t a{}, b{}, c{};
        for (int i = 0; i < MP_N64; i++) a[i] = rng();

        for (uint64_t k : {0ULL, 1ULL, 2ULL, 7ULL, 31ULL, 63ULL, 64ULL, 65ULL, 127ULL, 191ULL, 255ULL}) {
            mp_shl(b, a, k);
            mp_shr(c, b, k);

            // shifting left then right is not perfect if bits shifted out; compare against reference
            cpp_int av = to_cpp_int(a, MP_N64);
            cpp_int ref = (av << k) & ((cpp_int(1) << (64 * MP_N64)) - 1);
            ref >>= k;

            mp_uint_t ref_mp{};
            from_cpp_int(ref_mp, MP_N64, ref);
            expect_eq_mp(c, ref_mp);
        }
    }
}

TEST(mp, add_sub_carry_borrow) {
    mp_uint_t a{}, b{}, r{}, s{};
    // Max value to force carry
    for (int i = 0; i < MP_N64; i++) a[i] = 0xFFFFFFFFFFFFFFFFULL;
    set_u64(b, 1);

    uint64_t carry = mp_add(r, a, b);
    EXPECT_EQ(carry, 1ULL);
    EXPECT_TRUE(mp_is_zero(r));

    uint64_t borrow = mp_sub(s, r, 1ULL);
    // 0 - 1 underflows: borrow=1, result=2^(256)-1 (for MP_N64=4)
    EXPECT_EQ(borrow, 1ULL);
    expect_eq_mp(s, a);
}

TEST(mp, add_sub_bigint_fuzz) {
    std::mt19937_64 rng(7);
    for (int t = 0; t < 500; t++) {
        mp_uint_t a{}, b{}, r{}, s{};
        for (int i = 0; i < MP_N64; i++) { a[i] = rng(); b[i] = rng(); }

        cpp_int av = to_cpp_int(a, MP_N64);
        cpp_int bv = to_cpp_int(b, MP_N64);

        uint64_t carry = mp_add(r, a, b);
        cpp_int sum = av + bv;
        cpp_int mask = (cpp_int(1) << (64 * MP_N64)) - 1;
        cpp_int sum_lo = sum & mask;
        uint64_t ref_carry = ((sum >> (64 * MP_N64)) & cpp_int(0xFFFFFFFFFFFFFFFFULL)).convert_to<uint64_t>();

        mp_uint_t ref_sum{};
        from_cpp_int(ref_sum, MP_N64, sum_lo);

        expect_eq_mp(r, ref_sum);
        EXPECT_EQ(carry, ref_carry);

        uint64_t borrow = mp_sub(s, r, b);
        // since s = (a+b)-b mod 2^N, should equal a and borrow should equal carry? Not necessarily due to wrap,
        // but if we do exact arithmetic: (sum_lo - bv) underflows iff bv > sum_lo.
        cpp_int diff = sum_lo - bv;
        uint64_t ref_borrow = (diff < 0) ? 1ULL : 0ULL;
        if (diff < 0) diff += (cpp_int(1) << (64 * MP_N64));

        mp_uint_t ref_diff{};
        from_cpp_int(ref_diff, MP_N64, diff);
        expect_eq_mp(s, ref_diff);
        EXPECT_EQ(borrow, ref_borrow);
    }
}

TEST(mp, tstbit) {
    mp_uint_t a{};
    set_u64(a, 0);
    EXPECT_EQ(mp_tstbit(a, 0), 0);

    // set bit 0 and bit 127 and bit 255
    set_u64(a, 1);
    EXPECT_EQ(mp_tstbit(a, 0), 1);
    EXPECT_EQ(mp_tstbit(a, 1), 0);

    mp_uint_t b{};
    set_u64(b, 0);
    b[1] = (1ULL << 63); // bit 127
    EXPECT_EQ(mp_tstbit(b, 127), 1);
    EXPECT_EQ(mp_tstbit(b, 126), 0);

    mp_uint_t c{};
    set_u64(c, 0);
    c[3] = (1ULL << 63); // bit 255
    EXPECT_EQ(mp_tstbit(c, 255), 1);
    EXPECT_EQ(mp_tstbit(c, 254), 0);
}

TEST(mp, export_import_be_roundtrip) {
    std::mt19937_64 rng(123);
    for (int t = 0; t < 200; t++) {
        mp_uint_t a{}, b{};
        for (int i = 0; i < MP_N64; i++) a[i] = rng();
        auto bytes = export_be(a);
        import_be(b, bytes);
        expect_eq_mp(a, b);
    }

    // all-zero should serialize to all-zero
    mp_uint_t z{}, z2{};
    set_u64(z, 0);
    auto zb = export_be(z);
    for (auto c : zb) EXPECT_EQ(c, 0);
    import_be(z2, zb);
    expect_eq_mp(z, z2);
}

TEST(mp, set_str_and_back) {
    mp_uint_t a{};
    ASSERT_TRUE(mp_set(a, "0", 10));
    EXPECT_TRUE(mp_is_zero(a));
    EXPECT_EQ(mp_set_str(a, 10), "0");

    ASSERT_TRUE(mp_set(a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16));
    // should fill all limbs with 0xFF...
    for (int i = 0; i < MP_N64; i++) EXPECT_EQ(a[i], 0xFFFFFFFFFFFFFFFFULL);

    // roundtrip for random values in multiple bases
    std::mt19937_64 rng(9);
    for (int t = 0; t < 200; t++) {
        mp_uint_t x{}, y{};
        for (int i = 0; i < MP_N64; i++) x[i] = rng();

        for (uint32_t base : {2u, 10u, 16u}) {
            std::string s = mp_set_str(x, base);
            ASSERT_TRUE(mp_set(y, s.c_str(), base));
            expect_eq_mp(x, y);
        }
    }
}

TEST(mp, fits_get_int32) {
    mp_uint_t a{};
    set_u64(a, 0);
    EXPECT_TRUE(mp_fits_int32(a));
    EXPECT_EQ(mp_get_int32(a), 0);

    set_u64(a, 123);
    EXPECT_TRUE(mp_fits_int32(a));
    EXPECT_EQ(mp_get_int32(a), 123);

    // INT_MAX fits, INT_MAX+1 doesn't
    set_u64(a, (uint64_t)INT32_MAX);
    EXPECT_TRUE(mp_fits_int32(a));
    set_u64(a, (uint64_t)INT32_MAX + 1ULL);
    EXPECT_FALSE(mp_fits_int32(a));

    // any higher limb set -> false
    set_u64(a, 0);
    a[1] = 1;
    EXPECT_FALSE(mp_fits_int32(a));
}

TEST(mp, add_mod_and_set_mod_str) {
    // small prime modulus 97
    mp_uint_t mod{};
    set_u64(mod, 97);

    std::mt19937_64 rng(11);
    for (int t = 0; t < 300; t++) {
        mp_uint_t a{}, b{}, r{};
        set_u64(a, rng() % 5000);
        set_u64(b, rng() % 5000);

        cpp_int av = to_cpp_int(a, MP_N64);
        cpp_int bv = to_cpp_int(b, MP_N64);
        cpp_int mv = to_cpp_int(mod, MP_N64);

        mp_add_mod(r, a, b, mod);
        cpp_int ref_add = (av + bv) % mv;

        mp_uint_t ref_add_mp{};
        from_cpp_int(ref_add_mp, MP_N64, ref_add);
        expect_eq_mp(r, ref_add_mp);
    }

    mp_uint_t a{};
    ASSERT_TRUE(mp_set_mod(a, "-5", 10, mod));
    // -5 mod 97 = 92
    mp_uint_t exp{};
    set_u64(exp, 92);
    expect_eq_mp(a, exp);
}

// separate test for mp_mul_mod to avoid unused to avoid unused
TEST(mp, mul_mod_correctness) {
    mp_uint_t mod{};
    set_u64(mod, 97);
    std::mt19937_64 rng(12);

    for (int t = 0; t < 300; t++) {
        mp_uint_t a{}, r{};
        uint32_t b = (uint32_t)(rng() % 1000);
        set_u64(a, rng() % 50000);

        cpp_int av = to_cpp_int(a, MP_N64);
        cpp_int mv = to_cpp_int(mod, MP_N64);

        mp_mul_mod(r, a, b, mod);
        cpp_int ref = (av * b) % mv;

        mp_uint_t ref_mp{};
        from_cpp_int(ref_mp, MP_N64, ref);
        expect_eq_mp(r, ref_mp);
    }
}

TEST(mp, div_qr_big) {
    std::mt19937_64 rng(202);
    for (int t = 0; t < 200; t++) {
        mp_uint_t num{}, den{}, q{}, r{};
        // ensure den != 0
        for (int i = 0; i < MP_N64; i++) { num[i] = rng(); den[i] = rng(); }
        if (mp_is_zero(den)) den[0] = 1;

        cpp_int nv = to_cpp_int(num, MP_N64);
        cpp_int dv = to_cpp_int(den, MP_N64);
        if (dv == 0) dv = 1;

        bool ok = mp_div(q, r, num, den);
        ASSERT_TRUE(ok);

        cpp_int qv = nv / dv;
        cpp_int rv = nv % dv;

        mp_uint_t qref{}, rref{};
        from_cpp_int(qref, MP_N64, qv);
        from_cpp_int(rref, MP_N64, rv);
        expect_eq_mp(q, qref);
        expect_eq_mp(r, rref);
    }
}

TEST(mp, pow_mod_and_set_mod_i64) {
    mp_uint_t mod{};
    set_u64(mod, 97);
    cpp_int mv = to_cpp_int(mod, MP_N64);

    mp_uint_t base{}, exp{}, out{};
    set_u64(base, 5);
    set_u64(exp, 117); // 5^117 mod 97
    mp_pow_mod(out, base, exp, mod);

    cpp_int ref = mod_pow(cpp_int(5), cpp_int(117), mv);
    mp_uint_t ref_mp{};
    from_cpp_int(ref_mp, MP_N64, ref);
    expect_eq_mp(out, ref_mp);

    // mp_set_mod(int64_t)
    mp_uint_t x{};
    set_i64_mod(x, -5, mod);
    mp_uint_t expect{};
    set_u64(expect, 92);
    expect_eq_mp(x, expect);

    set_i64_mod(x, 194, mod); // 194 mod 97 = 0
    set_u64(expect, 0);
    expect_eq_mp(x, expect);
}

TEST(mp, mul_addmul_and_bitwise) {
    std::mt19937_64 rng(333);
    for (int t = 0; t < 300; t++) {
        mp_uint_t a{}, r{}, r2{};
        for (int i = 0; i < MP_N64; i++) a[i] = rng();
        uint64_t b = rng();

        cpp_int av = to_cpp_int(a, MP_N64);
        cpp_int prod = av * b;

        uint64_t carry = mp_mul(r, a, b);
        cpp_int mask = (cpp_int(1) << (64 * MP_N64)) - 1;
        cpp_int prod_lo = prod & mask;
        uint64_t ref_carry = ((prod >> (64 * MP_N64)) & cpp_int(0xFFFFFFFFFFFFFFFFULL)).convert_to<uint64_t>();

        mp_uint_t ref_mp{};
        from_cpp_int(ref_mp, MP_N64, prod_lo);
        expect_eq_mp(r, ref_mp);
        EXPECT_EQ(carry, ref_carry);

        // addmul: r2 = r + a*b
        mp_copy(r2, r);
        uint64_t carry2 = mp_addmul(r2, a, b);
        cpp_int sum = prod + prod_lo; // r == low(prod), addmul adds full (a*b)
        uint64_t ref_carry2 = ((sum >> (64 * MP_N64)) & cpp_int(0xFFFFFFFFFFFFFFFFULL)).convert_to<uint64_t>();
        cpp_int sum_lo = sum & mask;

        mp_uint_t sum_ref{};
        from_cpp_int(sum_ref, MP_N64, sum_lo);
        expect_eq_mp(r2, sum_ref);
        EXPECT_EQ(carry2, ref_carry2);

        // bitwise ops
        mp_uint_t x{}, y{}, z{};
        for (int i = 0; i < MP_N64; i++) { x[i] = rng(); y[i] = rng(); }

        mp_and(z, x, y);
        for (int i = 0; i < MP_N64; i++) EXPECT_EQ(z[i], (x[i] & y[i]));

        mp_or(z, x, y);
        for (int i = 0; i < MP_N64; i++) EXPECT_EQ(z[i], (x[i] | y[i]));

        mp_xor(z, x, y);
        for (int i = 0; i < MP_N64; i++) EXPECT_EQ(z[i], (x[i] ^ y[i]));

        mp_not(z, x);
        for (int i = 0; i < MP_N64; i++) EXPECT_EQ(z[i], ~x[i]);
    }
}

TEST(mp, add_variable_limbs) {
    // test n=1..4 with reference
    std::mt19937_64 rng(1001);
    for (size_t an = 1; an <= (size_t)MP_N64; an++) {
        for (size_t bn = 1; bn <= (size_t)MP_N64; bn++) {
            for (int t = 0; t < 100; t++) {
                std::vector<uint64_t> a(an), b(bn), r(std::max(an, bn));
                for (size_t i = 0; i < an; i++) a[i] = rng();
                for (size_t i = 0; i < bn; i++) b[i] = rng();

                uint64_t carry = mp_add(r.data(), a.data(), an, b.data(), bn);

                cpp_int av = to_cpp_int(a.data(), an);
                cpp_int bv = to_cpp_int(b.data(), bn);
                cpp_int sum = av + bv;
                cpp_int mask = (cpp_int(1) << (64 * std::max(an, bn))) - 1;
                cpp_int sum_lo = sum & mask;
                uint64_t ref_carry = ((sum >> (64 * std::max(an, bn))) & cpp_int(0xFFFFFFFFFFFFFFFFULL)).convert_to<uint64_t>();

                std::vector<uint64_t> ref(std::max(an, bn));
                from_cpp_int(ref.data(), std::max(an, bn), sum_lo);

                ASSERT_EQ(carry, ref_carry);
                for (size_t i = 0; i < ref.size(); i++) EXPECT_EQ(r[i], ref[i]);
            }
        }
    }
}

TEST(mp, inv_mod) {
    // prime modulus 97
    mp_uint_t mod{};
    set_u64(mod, 97);
    cpp_int mv = to_cpp_int(mod, MP_N64);

    for (uint64_t a_u : {1ULL, 2ULL, 3ULL, 5ULL, 96ULL}) {
        mp_uint_t a{}, inv{};
        set_u64(a, a_u);
        bool ok = mp_inv_mod(inv, a, mod);
        ASSERT_TRUE(ok);

        cpp_int av = cpp_int(a_u);
        cpp_int ref = mod_pow(av, mv - 2, mv);

        mp_uint_t ref_mp{};
        from_cpp_int(ref_mp, MP_N64, ref);
        expect_eq_mp(inv, ref_mp);

        // check (a*inv) % mod == 1
        cpp_int check = (av * ref) % mv;
        EXPECT_EQ(check, 1);
    }

    // non-invertible: a=0
    mp_uint_t zero{}, inv{};
    set_u64(zero, 0);
    EXPECT_FALSE(mp_inv_mod(inv, zero, mod));
}

} // namespace
