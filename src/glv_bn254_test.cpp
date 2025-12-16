// src/glv_bn254_test.cpp
#include <gtest/gtest.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <random>
#include <gmp.h>
#include "alt_bn128.hpp"
#include "glv_bn254_nogmp.hpp"

// ---------------- helpers ----------------
static inline void mpz_from_limbs_le(mpz_t out, const mp_limb_t* limbs, size_t n)
{
    mpz_import(
        out,
        n,
        -1,
        sizeof(mp_limb_t),
        0,
        0,
        limbs
    );
}

static inline void mpz_from_bytes_le(mpz_t out, const uint8_t *bytes, size_t len)
{
    mpz_import(
        out,
        len,
        -1,
        1,
        0,
        0,
        bytes
    );
}

static inline void mpz_to_bytes32_le(uint8_t out[32], const mpz_t in)
{
    std::memset(out, 0, 32);

    size_t count = 0;
    uint8_t tmp[64];
    std::memset(tmp, 0, sizeof(tmp));

    mpz_export(tmp, &count, -1, 1, 0, 0, in);
    if (count > 32) count = 32;
    std::memcpy(out, tmp, count);
}

static inline int signed_bitlen(const mpz_t x)
{
    if (mpz_sgn(x) == 0) return 0;
    mpz_t tmp;
    mpz_init(tmp);
    mpz_abs(tmp, x);
    int bl = mpz_sizeinbase(tmp, 2);
    mpz_clear(tmp);
    return bl;
}

static inline bool g1_eq_affine(
    AltBn128::Engine& E,
    const AltBn128::Engine::G1::PointAffine& a,
    const AltBn128::Engine::G1::PointAffine& b)
{
    auto aa = a;
    auto bb = b;
    return E.g1.eq(aa, bb);
}

static inline bool g1_eq(
    AltBn128::Engine& E,
    AltBn128::Engine::G1::Point& a,
    AltBn128::Engine::G1::Point& b)
{
    return E.g1.eq(a, b);
}

static inline AltBn128::Engine::G1::Point g1_mul_le(
    AltBn128::Engine& E,
    AltBn128::Engine::G1::PointAffine base,
    const uint8_t* s_le,
    uint64_t s_len)
{
    AltBn128::Engine::G1::Point r;
    E.g1.mulByScalar(r, base, (uint8_t*)s_le, s_len);
    return r;
}

// Calculate lambda from (a1,b1) using the formula: a1 + b1*lambda ≡ 0 (mod N)
// We only have abs(b1), so try both signs and choose the one
// that actually matches the endomorphism on the generator
static inline void compute_lambda_le32(uint8_t lambda_le[32])
{
    using namespace glv_bn254;

    mpz_t N, absb1, inv, a1, lam_pos, lam_neg;
    mpz_inits(N, absb1, inv, a1, lam_pos, lam_neg, nullptr);

    mpz_from_limbs_le(N, FR_N, N64);
    mpz_from_limbs_le(absb1, ABS_B1, 2);
    mpz_set_ui(a1, (unsigned long)A1);

    // inv = abs(b1)^(-1) mod N
    ASSERT_NE(mpz_invert(inv, absb1, N), 0) << "mpz_invert(abs(b1), N) failed";

    // b1 = -abs(b1)  => lambda =  a1 * inv mod N
    mpz_mul(lam_pos, a1, inv);
    mpz_mod(lam_pos, lam_pos, N);

    // b1 = +abs(b1)  => lambda = -a1 * inv mod N = N - (a1*inv mod N)
    mpz_sub(lam_neg, N, lam_pos);
    mpz_mod(lam_neg, lam_neg, N);

    // [lambda]G == phi(G)
    AltBn128::Engine E;

    AltBn128::Engine::G1::PointAffine G = E.g1.oneAffine();
    AltBn128::Engine::G1::PointAffine phiG = G;
    glv_bn254::apply_phi_inplace_g1(phiG);

    // [lam_pos]G
    uint8_t lam_pos_le[32];
    mpz_to_bytes32_le(lam_pos_le, lam_pos);
    auto lamPosG = g1_mul_le(E, G, lam_pos_le, 32);

    // [lam_neg]G
    uint8_t lam_neg_le[32];
    mpz_to_bytes32_le(lam_neg_le, lam_neg);
    auto lamNegG = g1_mul_le(E, G, lam_neg_le, 32);

    // affine
    AltBn128::Engine::G1::PointAffine a1_aff, a2_aff;
    E.g1.copy(a1_aff, lamPosG);
    E.g1.copy(a2_aff, lamNegG);

    bool ok_pos = g1_eq_affine(E, a1_aff, phiG);
    bool ok_neg = g1_eq_affine(E, a2_aff, phiG);

    ASSERT_TRUE(ok_pos || ok_neg) << "Neither lambda candidate matches phi(G).";

    if (ok_pos) {
        std::memcpy(lambda_le, lam_pos_le, 32);
    } else {
        std::memcpy(lambda_le, lam_neg_le, 32);
    }

    mpz_clears(N, absb1, inv, a1, lam_pos, lam_neg, nullptr);
}

// ---------------- tests ----------------

TEST(GlvBn254, PhiMatchesLambdaOnGenerator)
{
    AltBn128::Engine E;

    AltBn128::Engine::G1::PointAffine G = E.g1.oneAffine();

    AltBn128::Engine::G1::PointAffine phiG = G;
    glv_bn254::apply_phi_inplace_g1(phiG);

    uint8_t lambda_le[32];
    compute_lambda_le32(lambda_le);

    auto lamG = g1_mul_le(E, G, lambda_le, 32);

    AltBn128::Engine::G1::PointAffine lamG_aff;
    E.g1.copy(lamG_aff, lamG);

    EXPECT_TRUE(g1_eq_affine(E, lamG_aff, phiG));
}

static inline void check_decomp_reconstruct(AltBn128::Engine& E, const uint8_t k_le[32])
{
    AltBn128::Engine::G1::PointAffine G = E.g1.oneAffine();

    // kG
    auto kG = g1_mul_le(E, G, k_le, 32);

    // decompose(k) -> (k1,k2,signs)
    auto d = glv_bn254::decompose_fr_le_32(k_le);

    // phi(G)
    AltBn128::Engine::G1::PointAffine phiG = G;
    glv_bn254::apply_phi_inplace_g1(phiG);

    // p1 = (+/-)k1 * G
    auto p1 = g1_mul_le(E, G, d.k1, 16);
    if (d.neg1) E.g1.neg(p1, p1);

    // p2 = (+/-)k2 * phi(G)
    auto p2 = g1_mul_le(E, phiG, d.k2, 16);
    if (d.neg2) E.g1.neg(p2, p2);

    // sum = p1 + p2
    AltBn128::Engine::G1::Point sum;
    E.g1.add(sum, p1, p2);

    EXPECT_TRUE(g1_eq(E, sum, kG));
}

TEST(GlvBn254, DecomposeReconstruct_FixedScalars)
{
    AltBn128::Engine E;

    // 0
    {
        uint8_t k[32] = {0};
        check_decomp_reconstruct(E, k);
    }

    // 1
    {
        uint8_t k[32] = {0};
        k[0] = 1;
        check_decomp_reconstruct(E, k);
    }

    {
        uint8_t k[32] = {
            0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,
            0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x10,
            0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,
            0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x20
        };
        check_decomp_reconstruct(E, k);
    }

    {
        uint8_t k[32] = {
            0xff,0xee,0xdd,0xcc,0xbb,0xaa,0x99,0x88,
            0x77,0x66,0x55,0x44,0x33,0x22,0x11,0x00,
            0x10,0x32,0x54,0x76,0x98,0xba,0xdc,0xfe,
            0x01,0x23,0x45,0x67,0x89,0xab,0xcd,0xef
        };
        check_decomp_reconstruct(E, k);
    }
}

TEST(GlvBn254, DecomposeReconstruct_RandomDeterministic)
{
    AltBn128::Engine E;

    std::mt19937_64 rng(0xC0FFEE123456789ULL);

    for (int t = 0; t < 64; t++) {
        uint8_t k[32];
        for (int i = 0; i < 32; i++) {
            k[i] = (uint8_t)(rng() & 0xff);
        }
        check_decomp_reconstruct(E, k);
    }
}

TEST(GlvBn254, LambdaMatchesReferenceFromPython)
{
    using namespace glv_bn254;

    mpz_t N;
    mpz_init(N);
    mpz_from_limbs_le(N, FR_N, N64);

    uint8_t lambda_le[32];
    compute_lambda_le32(lambda_le);

    mpz_t lambda_cpp;
    mpz_init(lambda_cpp);
    mpz_from_bytes_le(lambda_cpp, lambda_le, 32);
    mpz_mod(lambda_cpp, lambda_cpp, N);

    // lamb = 4407920970296243842393367215006156084916469457145843978461
    mpz_t lambda_ref;
    mpz_init(lambda_ref);
    mpz_set_str(lambda_ref,
        "4407920970296243842393367215006156084916469457145843978461", 10);
    mpz_mod(lambda_ref, lambda_ref, N);

    EXPECT_EQ(mpz_cmp(lambda_cpp, lambda_ref), 0)
        << "Lambda from C++ does not match reference Python lambda";

    mpz_clear(lambda_cpp);
    mpz_clear(lambda_ref);
    mpz_clear(N);
}

static void check_lattice_example(const char* name, const char* k_dec_str)
{
    using namespace glv_bn254;

    mpz_t N, lambda_mpz, k;
    mpz_inits(N, lambda_mpz, k, nullptr);

    mpz_from_limbs_le(N, FR_N, N64);

    uint8_t lambda_le[32];
    compute_lambda_le32(lambda_le);
    mpz_from_bytes_le(lambda_mpz, lambda_le, 32);
    mpz_mod(lambda_mpz, lambda_mpz, N);

    ASSERT_EQ(mpz_set_str(k, k_dec_str, 10), 0) << "Bad decimal string for k";
    mpz_mod(k, k, N);

    uint8_t k_le[32];
    mpz_to_bytes32_le(k_le, k);

    auto d = glv_bn254::decompose_fr_le_32(k_le);

    mpz_t k1, k2;
    mpz_inits(k1, k2, nullptr);

    mpz_from_bytes_le(k1, d.k1, 16);
    mpz_from_bytes_le(k2, d.k2, 16);

    if (d.neg1) mpz_neg(k1, k1);
    if (d.neg2) mpz_neg(k2, k2);

    // lhs = k1 + k2*lambda (mod N)
    mpz_t lhs;
    mpz_init(lhs);
    mpz_mul(lhs, k2, lambda_mpz);
    mpz_add(lhs, lhs, k1);
    mpz_mod(lhs, lhs, N);

    EXPECT_EQ(mpz_cmp(lhs, k), 0) << name << ": (k1 + k2*lambda) != k (mod N)";

    int bl1 = signed_bitlen(k1);
    int bl2 = signed_bitlen(k2);

    EXPECT_LE(bl1, 130) << name << ": k1 too large in bits";
    EXPECT_LE(bl2, 130) << name << ": k2 too large in bits";

    mpz_clears(N, lambda_mpz, k, k1, k2, lhs, nullptr);
}

TEST(GlvBn254, LatticeExamples_FromPythonScript)
{
    // Example 1
    check_lattice_example(
        "Example 1",
        "21877113871839275222246405745257275088548364400416034343698204186575808497898"
    );

    // Example 2
    check_lattice_example(
        "Example 2",
        "78328712387328791256126782316736712671278938271389021"
    );

    // Example 3
    check_lattice_example(
        "Example 3",
        "563178213881238123"
    );
}

TEST(GlvBn254, PhiMatchesLambdaOn2G)
{
    AltBn128::Engine E;

    AltBn128::Engine::G1::PointAffine G = E.g1.oneAffine();
    AltBn128::Engine::G1::Point P2;
    AltBn128::Engine::G1::PointAffine P2_aff;

    // P2 = 2G
    E.g1.mulByScalar(P2, G, (uint8_t*)"\x02", 1);
    E.g1.copy(P2_aff, P2);

    // phi(P2)
    AltBn128::Engine::G1::PointAffine phiP2 = P2_aff;
    glv_bn254::apply_phi_inplace_g1(phiP2);

    // lambda
    uint8_t lambda_le[32];
    compute_lambda_le32(lambda_le);

    // [lambda]P2
    auto lamP2 = g1_mul_le(E, P2_aff, lambda_le, 32);
    AltBn128::Engine::G1::PointAffine lamP2_aff;
    E.g1.copy(lamP2_aff, lamP2);

    EXPECT_TRUE(g1_eq_affine(E, lamP2_aff, phiP2));
}
