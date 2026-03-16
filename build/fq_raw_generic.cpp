#include "fq_element.hpp"
#include "mp.hpp"
#include <cstring>

static uint64_t     Fq_rawq[] = {0x3c208c16d87cfd47,0x97816a916871ca8d,0xb85045b68181585d,0x30644e72e131a029, 0};
static const uint64_t Fq_np   =  0x87d20782e4866389ULL;
static uint64_t     lboMask   =  0x3fffffffffffffff;


void Fq_rawAdd(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    uint64_t carry = mp_add(pRawResult, pRawA, pRawB);

    if(carry || mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawAddLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t carry = mp_add(pRawResult, pRawA, rawB);

    if(carry || mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawSub(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    uint64_t carry = mp_sub(pRawResult, pRawA, pRawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawSubRegular(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    mp_sub(pRawResult, pRawA, pRawB);
}

void Fq_rawSubSL(FqRawElement pRawResult, uint64_t rawA, FqRawElement pRawB)
{
    FqRawElement pRawA = {rawA, 0, 0, 0};

    uint64_t carry = mp_sub(pRawResult, pRawA, pRawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawSubLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t carry = mp_sub(pRawResult, pRawA, rawB);

    if(carry)
    {
        mp_add(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawNeg(FqRawElement pRawResult, const FqRawElement pRawA)
{
    if (!mp_is_zero(pRawA))
    {
        mp_sub(pRawResult, Fq_rawq, pRawA);
    }
    else
    {
        mp_set(pRawResult, 0);
    }
}

//  Substracts a long element and a short element form 0
void Fq_rawNegLS(FqRawElement pRawResult, FqRawElement pRawA, uint64_t rawB)
{
    uint64_t carry1 = mp_sub(pRawResult, Fq_rawq, rawB);
    uint64_t carry2 = mp_sub(pRawResult, pRawResult, pRawA);

    if (carry1 || carry2)
    {
        mp_add(pRawResult, pRawResult, Fq_rawq);
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
    return mp_cmp(pRawA, pRawB) == 0;
}

void Fq_rawMMul(FqRawElement pRawResult, const FqRawElement pRawA, const FqRawElement pRawB)
{
    constexpr size_t N = Fq_N64 + 1;
    const uint64_t  *mq = Fq_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    product0[4] = mp_mul(product0, pRawB, pRawA[0]);

    np0 = Fq_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);

    product1[4] = mp_addmul(product1, pRawB, Fq_N64, pRawA[1]);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fq_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);

    product2[4] = mp_addmul(product2, pRawB, Fq_N64, pRawA[2]);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fq_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);

    product3[4] = mp_addmul(product3, pRawB, Fq_N64, pRawA[3]);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fq_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

void Fq_rawMMul1(FqRawElement pRawResult, const FqRawElement pRawA, uint64_t pRawB)
{
    constexpr size_t N = Fq_N64 + 1;
    const uint64_t  *mq = Fq_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    product0[4] = mp_mul(product0, pRawA, pRawB);

    np0 = Fq_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fq_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fq_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fq_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

void Fq_rawFromMontgomery(FqRawElement pRawResult, const FqRawElement &pRawA)
{
    constexpr size_t N = Fq_N64 + 1;
    const uint64_t  *mq = Fq_rawq;

    uint64_t  np0;

    uint64_t  product0[N];
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    mp_copy(product0, pRawA); product0[4] = 0;

    np0 = Fq_np * product0[0];
    product1[1] = mp_addmul(product0, mq, N, np0);
    mp_add(product1, product1, N, product0+1, N-1);

    np0 = Fq_np * product1[0];
    product2[1] = mp_addmul(product1, mq, N, np0);
    mp_add(product2, product2, N, product1+1, N-1);

    np0 = Fq_np * product2[0];
    product3[1] = mp_addmul(product2, mq, N, np0);
    mp_add(product3, product3, N, product2+1, N-1);

    np0 = Fq_np * product3[0];
    mp_addmul(product3, mq, N, np0);

    mp_copy(pRawResult, product3+1);

    if (mp_cmp(pRawResult, mq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, mq);
    }
}

int Fq_rawIsZero(const FqRawElement rawA)
{
    return mp_is_zero(rawA) ? 1 : 0;
}

int Fq_rawCmp(FqRawElement pRawA, FqRawElement pRawB)
{
    return mp_cmp(pRawA, pRawB);
}

void Fq_rawSwap(FqRawElement pRawResult, FqRawElement pRawA)
{
    FqRawElement temp;

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

void Fq_rawCopyS2L(FqRawElement pRawResult, int64_t val)
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

        mp_add(pRawResult, pRawResult, Fq_rawq);
    }
}


void Fq_rawAnd(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    mp_and(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawOr(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    mp_or(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawXor(FqRawElement pRawResult, FqRawElement pRawA, FqRawElement pRawB)
{
    mp_xor(pRawResult, pRawA, pRawB);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}

void Fq_rawShl(FqRawElement r, FqRawElement a, uint64_t b)
{
    mp_shl(r, a, (uint32_t)b);

    r[3] &= lboMask;

    if (mp_cmp(r, Fq_rawq) >= 0)
        mp_sub(r, r, Fq_rawq);
}

void Fq_rawShr(FqRawElement r, FqRawElement a, uint64_t b)
{
    mp_shr(r, a, (uint32_t)b);
}

void Fq_rawNot(FqRawElement pRawResult, FqRawElement pRawA)
{
    mp_not(pRawResult, pRawA);

    pRawResult[3] &= lboMask;

    if (mp_cmp(pRawResult, Fq_rawq) >= 0)
    {
        mp_sub(pRawResult, pRawResult, Fq_rawq);
    }
}