#include "fr_element.hpp"
#include "mp.hpp"
#include <cstring>

static uint64_t     Fr_rawq[] = {0x43e1f593f0000001,0x2833e84879b97091,0xb85045b68181585d,0x30644e72e131a029, 0};
static const uint64_t Fr_np   =  0xc2e1f593efffffffULL;
static uint64_t     lboMask   =  0x3fffffffffffffff;


void Fr_rawAdd(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    uint64_t carry = mp_add(pRawResult, pRawA, pRawB);

    if(carry || mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawAddLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    uint64_t carry = mp_add(pRawResult, pRawA, rawB);

    if(carry || mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawSub(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    uint64_t carry = mp_sub(pRawResult, pRawA, pRawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawSubRegular(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    mp_sub(pRawResult, pRawA, pRawB);
}

void Fr_rawSubSL(FrRawElement pRawResult, uint64_t rawA, FrRawElement pRawB)
{
    FrRawElement pRawA = {rawA, 0, 0, 0};

    uint64_t carry = mp_sub(pRawResult, pRawA, pRawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawSubLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    uint64_t carry = mp_sub(pRawResult, pRawA, rawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawNeg(FrRawElement pRawResult, const FrRawElement pRawA)
{
    if (!mp_is_zero(pRawA))
    {
        mp_sub(pRawResult, Fr_rawq, pRawA);
    }
    else
    {
        mp_set(pRawResult, 0);
    }
}

//  Substracts a long element and a short element form 0
void Fr_rawNegLS(FrRawElement pRawResult, FrRawElement pRawA, uint64_t rawB)
{
    uint64_t carry1 = mp_sub(pRawResult, Fr_rawq, rawB);
    uint64_t carry2 = mp_sub(pRawResult, pRawResult, pRawA);

    if (carry1 || carry2)
    {
        mp_add(pRawResult, pRawResult, Fr_rawq);
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
    return mp_cmp(pRawA, pRawB) == 0;
}

void Fr_rawMMul(FrRawElement pRawResult, const FrRawElement pRawA, const FrRawElement pRawB)
{
    constexpr size_t N = Fr_N64 + 1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    product0[4] = mp_mul(product0, pRawB, pRawA[0]);

    np0 = Fr_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);

    product1[4] = mp_addmul(product1, pRawB, Fr_N64, pRawA[1]);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fr_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);

    product2[4] = mp_addmul(product2, pRawB, Fr_N64, pRawA[2]);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fr_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);

    product3[4] = mp_addmul(product3, pRawB, Fr_N64, pRawA[3]);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fr_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

void Fr_rawMMul1(FrRawElement pRawResult, const FrRawElement pRawA, uint64_t pRawB)
{
    constexpr size_t N = Fr_N64 + 1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    product0[4] = mp_mul(product0, pRawA, pRawB);

    np0 = Fr_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fr_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fr_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fr_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

void Fr_rawFromMontgomery(FrRawElement pRawResult, const FrRawElement &pRawA)
{
    constexpr size_t N = Fr_N64 + 1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N];
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    mp_copy(product0, pRawA); product0[4] = 0;

    np0 = Fr_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fr_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fr_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fr_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

int Fr_rawIsZero(const FrRawElement rawA)
{
    return mp_is_zero(rawA) ? 1 : 0;
}

int Fr_rawCmp(FrRawElement pRawA, FrRawElement pRawB)
{
    return mp_cmp(pRawA, pRawB);
}

void Fr_rawSwap(FrRawElement pRawResult, FrRawElement pRawA)
{
    FrRawElement temp;

    temp[0] = pRawResult[0];
    temp[1] = pRawResult[1];
    temp[2] = pRawResult[2];
    temp[3] = pRawResult[3];

    pRawResult[0] = pRawA[0];
    pRawResult[1] = pRawA[1];
    pRawResult[2] = pRawA[2];
    pRawResult[3] = pRawA[3];

    pRawA[0] = temp[0];
    pRawA[1] = temp[1];
    pRawA[2] = temp[2];
    pRawA[3] = temp[3];
}

void Fr_rawCopyS2L(FrRawElement pRawResult, int64_t val)
{
    pRawResult[0] = val;
    pRawResult[1] = 0;
    pRawResult[2] = 0;
    pRawResult[3] = 0;

    if (val < 0)
    {
        pRawResult[1] = -1;
        pRawResult[2] = -1;
        pRawResult[3] = -1;

        mp_add(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawAnd(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    mp_and(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawOr(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    mp_or(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawXor(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    mp_xor(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}

void Fr_rawShl(FrRawElement r, FrRawElement a, uint64_t b)
{
    mp_shl(r, a, (uint32_t)b);

    r[3] &= lboMask;

    if (mp_cmp(r, Fr_rawq) >= 0)
        mp_sub(r, r, Fr_rawq);
}

void Fr_rawShr(FrRawElement r, FrRawElement a, uint64_t b)
{
    mp_shr(r, a, (uint32_t)b);
}

void Fr_rawNot(FrRawElement pRawResult, FrRawElement pRawA)
{
    mp_not(pRawResult, pRawA);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fr_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fr_rawq);
    }
}