#include "fq_element.hpp"
#include "u256.hpp"

#include <cstdint>
#include <cstring>
#include <cstddef>

#if !defined(__SIZEOF_INT128__)
#error "montgomery helpers require __int128"
#endif

static inline uint64_t limbs_mul_1(uint64_t *r, const uint64_t *a, size_t n, uint64_t b) {
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t t = (__uint128_t)a[i] * b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t limbs_addmul_1(uint64_t *r, const uint64_t *a, size_t n, uint64_t b) {
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t t = (__uint128_t)a[i] * b + r[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline void limbs_add_inplace(uint64_t *r, size_t rlen, const uint64_t *a, size_t alen) {
    __uint128_t carry = 0;
    size_t i = 0;

    for (; i < alen; i++) {
        __uint128_t t = (__uint128_t)r[i] + a[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    while (carry && i < rlen) {
        __uint128_t t = (__uint128_t)r[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
        i++;
    }
}

static inline int limbs_cmp(const uint64_t *a, const uint64_t *b, size_t n) {
    for (size_t i = n; i-- > 0;) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

static inline uint64_t limbs_sub_n(uint64_t *r, const uint64_t *a, const uint64_t *b, size_t n) {
    uint64_t borrow = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t ai = a[i];
        uint64_t bi = b[i];
        uint64_t ri = ai - bi - borrow;
        borrow = (borrow ? (ai <= bi) : (ai < bi));
        r[i] = ri;
    }
    return borrow;
}

static inline void mont_mul_4(uint64_t out[4], const uint64_t a[4], const uint64_t b[4],
                             const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[4] = limbs_mul_1(p0, b, 4, a[0]);

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);

    p1[4] = limbs_addmul_1(p1, b, 4, a[1]);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);

    p2[4] = limbs_addmul_1(p2, b, 4, a[2]);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);

    p3[4] = limbs_addmul_1(p3, b, 4, a[3]);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

static inline void mont_mul_4_u64(uint64_t out[4], const uint64_t a[4], uint64_t b,
                                 const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[4] = limbs_mul_1(p0, a, 4, b);

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

static inline void mont_reduce_4(uint64_t out[4], const uint64_t in_mont[4],
                                const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[0] = in_mont[0];
    p0[1] = in_mont[1];
    p0[2] = in_mont[2];
    p0[3] = in_mont[3];
    p0[4] = 0;

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

static uint64_t Fq_rawq[] = {
    0x3c208c16d87cfd47ULL,
    0x97816a916871ca8dULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL,
    0ULL
};

static constexpr uint64_t Fq_np   = 0x87d20782e4866389ULL;
static constexpr uint64_t lboMask = 0x3fffffffffffffffULL;

static const U256 Fq_qU256 = {{
    0x3c208c16d87cfd47ULL,
    0x97816a916871ca8dULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
}};

void Fq_rawAdd(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    uint64_t carry = mpn_add_n(pRawResult, pRawA, pRawB, 4);
    if (carry || mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawAddLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t carry = mpn_add_1(pRawResult, pRawA, 4, rawB);
    if (carry || mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawSub(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    uint64_t borrow = mpn_sub_n(pRawResult, pRawA, pRawB, 4);
    if (borrow) {
        (void)mpn_add_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawSubRegular(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    (void)mpn_sub_n(pRawResult, pRawA, pRawB, 4);
}

void Fq_rawSubSL(FqRawElement pRawResult, uint64_t rawA, FqRawElement pRawB)
{
    const uint64_t a[4] = { rawA, 0, 0, 0 };
    uint64_t borrow = mpn_sub_n(pRawResult, a, pRawB, 4);
    if (borrow) {
        (void)mpn_add_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawSubLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t borrow = mpn_sub_1(pRawResult, pRawA, 4, rawB);
    if (borrow) {
        (void)mpn_add_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawNeg(FqRawElement pRawResult, const FqRawElement pRawA)
{
    if (!mpn_zero_p(pRawA, 4)) {
        (void)mpn_sub_n(pRawResult, Fq_qU256.limb, pRawA, 4);
    } else {
        pRawResult[0] = pRawResult[1] = pRawResult[2] = pRawResult[3] = 0;
    }
}

void Fq_rawNegLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t t[4];
    (void)mpn_sub_1(t, Fq_qU256.limb, 4, rawB);    // t = q - rawB

    if (mpn_cmp(t, pRawA, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, t, pRawA, 4);
    } else {
        uint64_t tt[4];
        (void)mpn_add_n(tt, t, Fq_qU256.limb, 4);  // tt = t + q (без переполнения для BN254)
        (void)mpn_sub_n(pRawResult, tt, pRawA, 4);
        if (mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
            (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
        }
    }
}

void Fq_rawCopy(FqRawElement pRawResult, const FqRawElement pRawA)
{
    pRawResult[0] = pRawA[0];
    pRawResult[1] = pRawA[1];
    pRawResult[2] = pRawA[2];
    pRawResult[3] = pRawA[3];
}

int Fq_rawIsEq(const FqRawElement pRawA, const FqRawElement pRawB)
{
    return (pRawA[0] == pRawB[0] &&
            pRawA[1] == pRawB[1] &&
            pRawA[2] == pRawB[2] &&
            pRawA[3] == pRawB[3]) ? 1 : 0;
}

void Fq_rawMMul(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    mont_mul_4(pRawResult, pRawA, pRawB, Fq_rawq, Fq_np);
}

void Fq_rawMMul1(FqRawElement pRawResult, const FqRawElement pRawA, uint64_t pRawB)
{
    mont_mul_4_u64(pRawResult, pRawA, pRawB, Fq_rawq, Fq_np);
}

void Fq_rawFromMontgomery(FqRawElement pRawResult, const FqRawElement &pRawA)
{
    mont_reduce_4(pRawResult, pRawA, Fq_rawq, Fq_np);
}

int Fq_rawIsZero(const FqRawElement rawA)
{
    return ((rawA[0] | rawA[1] | rawA[2] | rawA[3]) == 0) ? 1 : 0;
}

int Fq_rawCmp(FqRawElement pRawA, FqRawElement pRawB)
{
    return mpn_cmp(pRawA, pRawB, 4);
}

void Fq_rawSwap(FqRawElement pRawResult, FqRawElement pRawA)
{
    uint64_t t0 = pRawResult[0], t1 = pRawResult[1], t2 = pRawResult[2], t3 = pRawResult[3];
    pRawResult[0] = pRawA[0]; pRawResult[1] = pRawA[1]; pRawResult[2] = pRawA[2]; pRawResult[3] = pRawA[3];
    pRawA[0] = t0; pRawA[1] = t1; pRawA[2] = t2; pRawA[3] = t3;
}

void Fq_rawCopyS2L(FqRawElement pRawResult, int64_t val)
{
    uint64_t tmp[4];
    tmp[0] = (uint64_t)val;
    tmp[1] = 0;
    tmp[2] = 0;
    tmp[3] = 0;

    if (val < 0) {
        tmp[1] = ~0ULL;
        tmp[2] = ~0ULL;
        tmp[3] = ~0ULL;
        (void)mpn_add_n(pRawResult, tmp, Fq_qU256.limb, 4);
        return;
    }
    pRawResult[0] = tmp[0]; pRawResult[1] = tmp[1]; pRawResult[2] = tmp[2]; pRawResult[3] = tmp[3];
}

void Fq_rawAnd(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    pRawResult[0] = pRawA[0] & pRawB[0];
    pRawResult[1] = pRawA[1] & pRawB[1];
    pRawResult[2] = pRawA[2] & pRawB[2];
    pRawResult[3] = (pRawA[3] & pRawB[3]) & lboMask;

    if (mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawOr(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    pRawResult[0] = pRawA[0] | pRawB[0];
    pRawResult[1] = pRawA[1] | pRawB[1];
    pRawResult[2] = pRawA[2] | pRawB[2];
    pRawResult[3] = (pRawA[3] | pRawB[3]) & lboMask;

    if (mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawXor(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    pRawResult[0] = pRawA[0] ^ pRawB[0];
    pRawResult[1] = pRawA[1] ^ pRawB[1];
    pRawResult[2] = pRawA[2] ^ pRawB[2];
    pRawResult[3] = (pRawA[3] ^ pRawB[3]) & lboMask;

    if (mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}

void Fq_rawShl(FqRawElement r, FqRawElement a, uint64_t b)
{
    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    uint64_t tmp[4];
    (void)mpn_lshift(tmp, a, 4, (unsigned)b);
    tmp[3] &= lboMask;
    if (mpn_cmp(tmp, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(tmp, tmp, Fq_qU256.limb, 4);
    }
    r[0] = tmp[0]; r[1] = tmp[1]; r[2] = tmp[2]; r[3] = tmp[3];
}

void Fq_rawShr(FqRawElement r, FqRawElement a, uint64_t b)
{
    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    uint64_t tmp[4];
    (void)mpn_rshift(tmp, a, 4, (unsigned)b);
    r[0] = tmp[0]; r[1] = tmp[1]; r[2] = tmp[2]; r[3] = tmp[3];
}

void Fq_rawNot(FqRawElement pRawResult, FqRawElement pRawA)
{
    pRawResult[0] = ~pRawA[0];
    pRawResult[1] = ~pRawA[1];
    pRawResult[2] = ~pRawA[2];
    pRawResult[3] = (~pRawA[3]) & lboMask;

    if (mpn_cmp(pRawResult, Fq_qU256.limb, 4) >= 0) {
        (void)mpn_sub_n(pRawResult, pRawResult, Fq_qU256.limb, 4);
    }
}
