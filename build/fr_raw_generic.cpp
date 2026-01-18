#include "fr_element.hpp"
#include "u256.hpp"

#include <cstdint>
#include <cstring>

static uint64_t Fr_rawq[] = {
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL,
    0ULL
};

static constexpr uint64_t Fr_np   = 0xc2e1f593efffffffULL;
static constexpr uint64_t lboMask = 0x3fffffffffffffffULL;

static const U256 Fr_q_u256 = {{
    0x43e1f593f0000001ULL,
    0x2833e84879b97091ULL,
    0xb85045b68181585dULL,
    0x30644e72e131a029ULL
}};

static inline void load_u256(U256* x, const FrRawElement a) {
    x->limb[0] = a[0];
    x->limb[1] = a[1];
    x->limb[2] = a[2];
    x->limb[3] = a[3];
}

static inline void store_u256(FrRawElement r, const U256* x) {
    r[0] = x->limb[0];
    r[1] = x->limb[1];
    r[2] = x->limb[2];
    r[3] = x->limb[3];
}

void Fr_rawAdd(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);

    uint64_t carry = u256_add(&r, &a, &b);
    if (carry || u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawAddLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    U256 a, r;
    load_u256(&a, pRawA);

    uint64_t carry = u256_add_ui(&r, &a, rawB);
    if (carry || u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawSub(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);

    uint64_t borrow = u256_sub(&r, &a, &b);
    if (borrow) {
        (void)u256_add(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawSubRegular(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);
    (void)u256_sub(&r, &a, &b);
    store_u256(pRawResult, &r);
}

void Fr_rawSubSL(FrRawElement pRawResult, uint64_t rawA, FrRawElement pRawB)
{
    U256 a, b, r;
    u256_set_ui(&a, rawA);
    load_u256(&b, pRawB);

    uint64_t borrow = u256_sub(&r, &a, &b);
    if (borrow) {
        (void)u256_add(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawSubLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    U256 a, r;
    load_u256(&a, pRawA);

    uint64_t borrow = u256_sub_ui(&r, &a, rawB);
    if (borrow) {
        (void)u256_add(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawNeg(FrRawElement pRawResult, const FrRawElement pRawA)
{
    U256 a, r;
    load_u256(&a, pRawA);

    if (!u256_is_zero(&a)) {
        (void)u256_sub(&r, &Fr_q_u256, &a);
        store_u256(pRawResult, &r);
    } else {
        pRawResult[0] = pRawResult[1] = pRawResult[2] = pRawResult[3] = 0;
    }
}

void Fr_rawNegLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    U256 a, t, r;
    load_u256(&a, pRawA);

    (void)u256_sub_ui(&t, &Fr_q_u256, rawB);

    if (u256_cmp(&t, &a) >= 0) {
        (void)u256_sub(&r, &t, &a);
    } else {
        U256 tt;
        (void)u256_add(&tt, &t, &Fr_q_u256);
        (void)u256_sub(&r, &tt, &a);
        if (u256_cmp(&r, &Fr_q_u256) >= 0) {
            (void)u256_sub(&r, &r, &Fr_q_u256);
        }
    }

    store_u256(pRawResult, &r);
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
    u256_mont_mul_4(pRawResult, pRawA, pRawB, Fr_rawq, Fr_np);
}

void Fr_rawMMul1(FrRawElement pRawResult, const FrRawElement pRawA, uint64_t pRawB)
{
    u256_mont_mul_4_u64(pRawResult, pRawA, pRawB, Fr_rawq, Fr_np);
}

void Fr_rawFromMontgomery(FrRawElement pRawResult, const FrRawElement &pRawA)
{
    u256_mont_reduce_4(pRawResult, pRawA, Fr_rawq, Fr_np);
}

int Fr_rawIsZero(const FrRawElement rawA)
{
    return ((rawA[0] | rawA[1] | rawA[2] | rawA[3]) == 0) ? 1 : 0;
}

int Fr_rawCmp(FrRawElement pRawA, FrRawElement pRawB)
{
    U256 a, b;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);
    return u256_cmp(&a, &b);
}

void Fr_rawSwap(FrRawElement pRawResult, FrRawElement pRawA)
{
    uint64_t t0 = pRawResult[0], t1 = pRawResult[1], t2 = pRawResult[2], t3 = pRawResult[3];
    pRawResult[0] = pRawA[0]; pRawResult[1] = pRawA[1]; pRawResult[2] = pRawA[2]; pRawResult[3] = pRawA[3];
    pRawA[0] = t0; pRawA[1] = t1; pRawA[2] = t2; pRawA[3] = t3;
}

void Fr_rawCopyS2L(FrRawElement pRawResult, int64_t val)
{
    U256 r;
    r.limb[0] = (uint64_t)val;
    r.limb[1] = 0;
    r.limb[2] = 0;
    r.limb[3] = 0;

    if (val < 0) {
        r.limb[1] = ~0ULL;
        r.limb[2] = ~0ULL;
        r.limb[3] = ~0ULL;
        (void)u256_add(&r, &r, &Fr_q_u256);
    }

    store_u256(pRawResult, &r);
}

void Fr_rawAnd(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] & pRawB[0];
    r.limb[1] = pRawA[1] & pRawB[1];
    r.limb[2] = pRawA[2] & pRawB[2];
    r.limb[3] = (pRawA[3] & pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawOr(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] | pRawB[0];
    r.limb[1] = pRawA[1] | pRawB[1];
    r.limb[2] = pRawA[2] | pRawB[2];
    r.limb[3] = (pRawA[3] | pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawXor(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] ^ pRawB[0];
    r.limb[1] = pRawA[1] ^ pRawB[1];
    r.limb[2] = pRawA[2] ^ pRawB[2];
    r.limb[3] = (pRawA[3] ^ pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}

void Fr_rawShl(FrRawElement r, FrRawElement a, uint64_t b)
{
    U256 A, R;
    load_u256(&A, a);

    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    u256_shl_2exp(&R, &A, (uint32_t)b);
    R.limb[3] &= lboMask;

    if (u256_cmp(&R, &Fr_q_u256) >= 0) {
        (void)u256_sub(&R, &R, &Fr_q_u256);
    }
    store_u256(r, &R);
}

void Fr_rawShr(FrRawElement r, FrRawElement a, uint64_t b)
{
    U256 A, R;
    load_u256(&A, a);

    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    u256_fdiv_q_2exp(&R, &A, (uint32_t)b);
    store_u256(r, &R);
}

void Fr_rawNot(FrRawElement pRawResult, FrRawElement pRawA)
{
    U256 r;
    r.limb[0] = ~pRawA[0];
    r.limb[1] = ~pRawA[1];
    r.limb[2] = ~pRawA[2];
    r.limb[3] = (~pRawA[3]) & lboMask;

    if (u256_cmp(&r, &Fr_q_u256) >= 0) {
        (void)u256_sub(&r, &r, &Fr_q_u256);
    }
    store_u256(pRawResult, &r);
}
