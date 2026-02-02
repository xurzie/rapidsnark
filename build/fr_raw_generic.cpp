#include "fr_element.hpp"
#include "mp.hpp"

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
    limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        limbs_sub_n(out, out, mod5, 4);
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
    limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        limbs_sub_n(out, out, mod5, 4);
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
    limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        limbs_sub_n(out, out, mod5, 4);
    }
}

static uint64_t Fr_rawq[] = {
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL,
    0ULL
};

static constexpr uint64_t Fr_np   = 0xc2e1f593efffffffULL;
static constexpr uint64_t lboMask = 0x3fffffffffffffffULL;

static const U256 Fr_qU256 = {{
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
}};

void Fr_rawAdd(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    uint64_t carry = mp_add(pRawResult, pRawA, pRawB);
    if (carry || mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

 void Fr_rawAddLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
 {
    uint64_t carry = mp_add(pRawResult, pRawA, rawB);
    if (carry || mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
     }
 }

void Fr_rawSub(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    uint64_t borrow = mp_sub(pRawResult, pRawA, pRawB);
    if (borrow) {
        mp_add(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawSubRegular(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    mp_sub(pRawResult, pRawA, pRawB);
}

void Fr_rawSubSL(FrRawElement pRawResult, uint64_t rawA, FrRawElement pRawB)
{
    const uint64_t a[4] = { rawA, 0, 0, 0 };
    uint64_t borrow = mp_sub(pRawResult, a, pRawB);
    if (borrow) {
        mp_add(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawSubLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    uint64_t borrow = mp_sub(pRawResult, pRawA, rawB);
    if (borrow) {
        mp_add(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawNeg(FrRawElement pRawResult, const FrRawElement pRawA)
{
    if (!mp_is_zero(pRawA)) {
        mp_sub(pRawResult, Fr_qU256.limb, pRawA);
    } else {
        pRawResult[0] = pRawResult[1] = pRawResult[2] = pRawResult[3] = 0;
    }
}

void Fr_rawNegLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    uint64_t t[4];
    mp_sub(t, Fr_qU256.limb, rawB);

    if (mp_cmp(t, pRawA) >= 0) {
        mp_sub(pRawResult, t, pRawA);
    } else {
        uint64_t tt[4];
        mp_add(tt, t, Fr_qU256.limb);
        mp_sub(pRawResult, tt, pRawA);
        if (mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
            mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
        }
    }
}

void Fr_rawCopy(FrRawElement pRawResult, const FrRawElement pRawA)
{
    pRawResult[0] = pRawA[0];
    pRawResult[1] = pRawA[1];
    pRawResult[2] = pRawA[2];
    pRawResult[3] = pRawA[3];
}

int Fr_rawIsEq(const FrRawElement pRawA, const FrRawElement pRawB)
{
    return (pRawA[0] == pRawB[0] &&
            pRawA[1] == pRawB[1] &&
            pRawA[2] == pRawB[2] &&
            pRawA[3] == pRawB[3]) ? 1 : 0;
}

void Fr_rawMMul(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    mont_mul_4(pRawResult, pRawA, pRawB, Fr_rawq, Fr_np);
}

void Fr_rawMMul1(FrRawElement pRawResult, const FrRawElement pRawA, uint64_t pRawB)
{
    mont_mul_4_u64(pRawResult, pRawA, pRawB, Fr_rawq, Fr_np);
}

void Fr_rawFromMontgomery(FrRawElement pRawResult, const FrRawElement &pRawA)
{
    mont_reduce_4(pRawResult, pRawA, Fr_rawq, Fr_np);
}

int Fr_rawIsZero(const FrRawElement rawA)
{
    return ((rawA[0] | rawA[1] | rawA[2] | rawA[3]) == 0) ? 1 : 0;
}

int Fr_rawCmp(FrRawElement pRawA, FrRawElement pRawB)
{
    return mp_cmp(pRawA, pRawB);
}

void Fr_rawSwap(FrRawElement pRawResult, FrRawElement pRawA)
{
    uint64_t t0 = pRawResult[0], t1 = pRawResult[1], t2 = pRawResult[2], t3 = pRawResult[3];
    pRawResult[0] = pRawA[0]; pRawResult[1] = pRawA[1]; pRawResult[2] = pRawA[2]; pRawResult[3] = pRawA[3];
    pRawA[0] = t0; pRawA[1] = t1; pRawA[2] = t2; pRawA[3] = t3;
}

void Fr_rawCopyS2L(FrRawElement pRawResult, int64_t val)
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
        mp_add(pRawResult, tmp, Fr_qU256.limb);
        return;
    }
    pRawResult[0] = tmp[0]; pRawResult[1] = tmp[1]; pRawResult[2] = tmp[2]; pRawResult[3] = tmp[3];
}

void Fr_rawAnd(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    pRawResult[0] = pRawA[0] & pRawB[0];
    pRawResult[1] = pRawA[1] & pRawB[1];
    pRawResult[2] = pRawA[2] & pRawB[2];
    pRawResult[3] = (pRawA[3] & pRawB[3]) & lboMask;

    if (mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawOr(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    pRawResult[0] = pRawA[0] | pRawB[0];
    pRawResult[1] = pRawA[1] | pRawB[1];
    pRawResult[2] = pRawA[2] | pRawB[2];
    pRawResult[3] = (pRawA[3] | pRawB[3]) & lboMask;

    if (mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawXor(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    pRawResult[0] = pRawA[0] ^ pRawB[0];
    pRawResult[1] = pRawA[1] ^ pRawB[1];
    pRawResult[2] = pRawA[2] ^ pRawB[2];
    pRawResult[3] = (pRawA[3] ^ pRawB[3]) & lboMask;

    if (mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
    }
}

void Fr_rawShl(FrRawElement r, FrRawElement a, uint64_t b)
{
    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }
    uint64_t tmp[4];
    mp_shl_2exp(tmp, a, (uint32_t)b);
    tmp[3] &= lboMask;
    if (mp_cmp(tmp, Fr_qU256.limb) >= 0) {
        mp_sub(tmp, tmp, Fr_qU256.limb);
    }
    r[0] = tmp[0]; r[1] = tmp[1]; r[2] = tmp[2]; r[3] = tmp[3];
}

void Fr_rawShr(FrRawElement r, FrRawElement a, uint64_t b)
{
    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }
    uint64_t tmp[4];
    mp_shr_2exp(tmp, a, (uint32_t)b);
    r[0] = tmp[0]; r[1] = tmp[1]; r[2] = tmp[2]; r[3] = tmp[3];
}

void Fr_rawNot(FrRawElement pRawResult, FrRawElement pRawA)
{
    pRawResult[0] = ~pRawA[0];
    pRawResult[1] = ~pRawA[1];
    pRawResult[2] = ~pRawA[2];
    pRawResult[3] = (~pRawA[3]) & lboMask;

    if (mp_cmp(pRawResult, Fr_qU256.limb) >= 0) {
        mp_sub(pRawResult, pRawResult, Fr_qU256.limb);
    }
}
