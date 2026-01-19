#include "fq_element.hpp"
#include "u256.hpp"

#include <cstdint>
#include <cstring>

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

static inline void load_u256(U256* out, const FqRawElement in) {
    out->limb[0] = in[0];
    out->limb[1] = in[1];
    out->limb[2] = in[2];
    out->limb[3] = in[3];
}

static inline void store_u256(FqRawElement out, const U256* a) {
    out[0] = a->limb[0];
    out[1] = a->limb[1];
    out[2] = a->limb[2];
    out[3] = a->limb[3];
}

void Fq_rawAdd(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);

    uint64_t carry = u256_add(&r, &a, &b);
    if (carry || u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawAddLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    U256 a, r;
    load_u256(&a, pRawA);

    uint64_t carry = u256_add_ui(&r, &a, rawB);
    if (carry || u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawSub(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);

    uint64_t borrow = u256_sub(&r, &a, &b);
    if (borrow) {
        (void)u256_add(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawSubRegular(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    U256 a, b, r;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);
    (void)u256_sub(&r, &a, &b);
    store_u256(pRawResult, &r);
}

void Fq_rawSubSL(FqRawElement pRawResult, uint64_t rawA, FqRawElement pRawB)
{
    U256 a, b, r;
    u256_set_ui(&a, rawA);
    load_u256(&b, pRawB);

    uint64_t borrow = u256_sub(&r, &a, &b);
    if (borrow) {
        (void)u256_add(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawSubLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    U256 a, r;
    load_u256(&a, pRawA);

    uint64_t borrow = u256_sub_ui(&r, &a, rawB);
    if (borrow) {
        (void)u256_add(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawNeg(FqRawElement pRawResult, const FqRawElement pRawA)
{
    U256 a, r;
    load_u256(&a, pRawA);

    if (!u256_is_zero(&a)) {
        (void)u256_sub(&r, &Fq_qU256, &a);
        store_u256(pRawResult, &r);
    } else {
        pRawResult[0] = pRawResult[1] = pRawResult[2] = pRawResult[3] = 0;
    }
}

void Fq_rawNegLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    U256 a, t, r;
    load_u256(&a, pRawA);

    (void)u256_sub_ui(&t, &Fq_qU256, rawB);

    if (u256_cmp(&t, &a) >= 0) {
        (void)u256_sub(&r, &t, &a);
    } else {
        U256 tt;
        (void)u256_add(&tt, &t, &Fq_qU256);
        (void)u256_sub(&r, &tt, &a);
        if (u256_cmp(&r, &Fq_qU256) >= 0) {
            (void)u256_sub(&r, &r, &Fq_qU256);
        }
    }

    store_u256(pRawResult, &r);
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
    u256_mont_mul_4(pRawResult, pRawA, pRawB, Fq_rawq, Fq_np);
}

void Fq_rawMMul1(FqRawElement pRawResult, const FqRawElement pRawA, uint64_t pRawB)
{
    u256_mont_mul_4_u64(pRawResult, pRawA, pRawB, Fq_rawq, Fq_np);
}

void Fq_rawFromMontgomery(FqRawElement pRawResult, const FqRawElement &pRawA)
{
    u256_mont_reduce_4(pRawResult, pRawA, Fq_rawq, Fq_np);
}

int Fq_rawIsZero(const FqRawElement rawA)
{
    return ((rawA[0] | rawA[1] | rawA[2] | rawA[3]) == 0) ? 1 : 0;
}

int Fq_rawCmp(FqRawElement pRawA, FqRawElement pRawB)
{
    U256 a, b;
    load_u256(&a, pRawA);
    load_u256(&b, pRawB);
    return u256_cmp(&a, &b);
}

void Fq_rawSwap(FqRawElement pRawResult, FqRawElement pRawA)
{
    uint64_t t0 = pRawResult[0], t1 = pRawResult[1], t2 = pRawResult[2], t3 = pRawResult[3];
    pRawResult[0] = pRawA[0]; pRawResult[1] = pRawA[1]; pRawResult[2] = pRawA[2]; pRawResult[3] = pRawA[3];
    pRawA[0] = t0; pRawA[1] = t1; pRawA[2] = t2; pRawA[3] = t3;
}

void Fq_rawCopyS2L(FqRawElement pRawResult, int64_t val)
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
        (void)u256_add(&r, &r, &Fq_qU256);
    }

    store_u256(pRawResult, &r);
}

void Fq_rawAnd(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] & pRawB[0];
    r.limb[1] = pRawA[1] & pRawB[1];
    r.limb[2] = pRawA[2] & pRawB[2];
    r.limb[3] = (pRawA[3] & pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawOr(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] | pRawB[0];
    r.limb[1] = pRawA[1] | pRawB[1];
    r.limb[2] = pRawA[2] | pRawB[2];
    r.limb[3] = (pRawA[3] | pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawXor(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    U256 r;
    r.limb[0] = pRawA[0] ^ pRawB[0];
    r.limb[1] = pRawA[1] ^ pRawB[1];
    r.limb[2] = pRawA[2] ^ pRawB[2];
    r.limb[3] = (pRawA[3] ^ pRawB[3]) & lboMask;

    if (u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}

void Fq_rawShl(FqRawElement r, FqRawElement a, uint64_t b)
{
    U256 A, R;
    load_u256(&A, a);

    if (b >= 256) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    u256_shl_2exp(&R, &A, (uint32_t)b);
    R.limb[3] &= lboMask;

    if (u256_cmp(&R, &Fq_qU256) >= 0) {
        (void)u256_sub(&R, &R, &Fq_qU256);
    }
    store_u256(r, &R);
}

void Fq_rawShr(FqRawElement r, FqRawElement a, uint64_t b)
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

void Fq_rawNot(FqRawElement pRawResult, FqRawElement pRawA)
{
    U256 r;
    r.limb[0] = ~pRawA[0];
    r.limb[1] = ~pRawA[1];
    r.limb[2] = ~pRawA[2];
    r.limb[3] = (~pRawA[3]) & lboMask;

    if (u256_cmp(&r, &Fq_qU256) >= 0) {
        (void)u256_sub(&r, &r, &Fq_qU256);
    }
    store_u256(pRawResult, &r);
}
