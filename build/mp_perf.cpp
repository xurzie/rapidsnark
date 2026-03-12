#include "mp.hpp"
#include "fq.hpp"

#include <chrono>
#include <cstdio>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>

#include "gtest/gtest.h"

static constexpr int LIMB_BITS = 2 * MP_N;

static inline uint64_t now_ns() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
}

static inline void set4(uint64_t r[MP_N64], uint64_t l0, uint64_t l1, uint64_t l2, uint64_t l3) {
    r[0] = l0;
    r[1] = l1;
    r[2] = l2;
    r[3] = l3;
}

static inline uint64_t splitmix64_next(uint64_t &state) {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t add_carry_ref(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
    uint64_t t0 = a + b;
    uint64_t carry1 = (t0 < a);
    uint64_t t1 = t0 + c;
    uint64_t carry2 = (t1 < t0);
    *out = t1;
    return carry1 | carry2;
}

static inline uint64_t sub_borrow_ref(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
    uint64_t t0 = a - b;
    uint64_t borrow1 = (a < b);
    uint64_t t1 = t0 - c;
    uint64_t borrow2 = (t0 < c);
    *out = t1;
    return borrow1 | borrow2;
}

static inline void mp_mul64_ref(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
    const uint64_t a0 = (uint32_t)a;
    const uint64_t a1 = a >> MP_N;
    const uint64_t b0 = (uint32_t)b;
    const uint64_t b1 = b >> MP_N;

    const uint64_t p00 = a0 * b0;
    const uint64_t p01 = a0 * b1;
    const uint64_t p10 = a1 * b0;
    const uint64_t p11 = a1 * b1;

    const uint64_t mid = (p00 >> MP_N) + (uint32_t)p01 + (uint32_t)p10;
    lo = (p00 & 0xFFFFFFFFULL) | (mid << MP_N);
    hi = p11 + (p01 >> MP_N) + (p10 >> MP_N) + (mid >> MP_N);
}

static inline unsigned mp_clz64_ref(uint64_t x) {
    return x ? (unsigned)__builtin_clzll(x) : (unsigned)LIMB_BITS;
}

static inline int mp_num_limbs_ref(const uint64_t *a) {
    for (int i = MP_N64 - 1; i >= 0; --i) {
        if (a[i] != 0) return i + 1;
    }
    return 0;
}

static inline int mp_bitlen_ref(const uint64_t *a) {
    const int n = mp_num_limbs_ref(a);
    if (n == 0) return 0;
    return (n - 1) * LIMB_BITS + (LIMB_BITS - (int)mp_clz64_ref(a[n - 1]));
}

static inline uint64_t div_2by1_ref(uint64_t *q, uint64_t hi, uint64_t lo, uint64_t v) {
    uint64_t qq  = 0;
    uint64_t rem = hi;

    for (int bit = LIMB_BITS - 1; bit >= 0; --bit) {
        uint64_t overflow = rem >> (LIMB_BITS - 1);
        rem = (rem << 1) | ((lo >> bit) & 1ULL);
        if (overflow || rem >= v) {
            rem -= v;
            qq |= (1ULL << bit);
        }
    }

    *q = qq;
    return rem;
}

static inline uint64_t div_1word_no128_ref(uint64_t *q, const uint64_t *u, uint64_t v) {
    uint64_t rem = 0;
    for (int i = MP_N64 - 1; i >= 0; --i) {
        uint64_t qi;
        rem = div_2by1_ref(&qi, rem, u[i], v);
        q[i] = qi;
    }
    return rem;
}

static inline uint64_t div_1word_int128_ref(uint64_t *q, const uint64_t *u, uint64_t v) {
    __uint128_t rem = 0;
    for (int i = MP_N64 - 1; i >= 0; --i) {
        __uint128_t cur = (rem << LIMB_BITS) | (__uint128_t)u[i];
        q[i] = (uint64_t)(cur / v);
        rem  = (uint64_t)(cur % v);
    }
    return (uint64_t)rem;
}

static inline uint64_t mul_sub_knuth_no128_ref(uint64_t *u, const uint64_t *v, int n, uint64_t qhat) {
    uint64_t carry  = 0;
    uint64_t borrow = 0;

    for (int i = 0; i < n; ++i) {
        uint64_t hi, lo;
        mp_mul64_ref(qhat, v[i], lo, hi);

        uint64_t c = add_carry_ref(&lo, lo, carry, 0);
        hi += c;

        borrow = sub_borrow_ref(&u[i], u[i], lo, borrow);
        carry = hi;
    }

    return sub_borrow_ref(&u[n], u[n], carry, borrow);
}

static inline uint64_t mul_sub_knuth_int128_ref(uint64_t *u, const uint64_t *v, int n, uint64_t qhat) {
    __uint128_t carry = 0;
    uint64_t borrow = 0;

    for (int i = 0; i < n; ++i) {
        __uint128_t prod = (__uint128_t)qhat * (__uint128_t)v[i] + carry;
        uint64_t pl = (uint64_t)prod;
        carry = prod >> LIMB_BITS;
        __int128_t t = (__int128_t)u[i] - (__int128_t)pl - (__int128_t)borrow;
        u[i] = (uint64_t)t;
        borrow = t < 0 ? 1u : 0u;
    }

    __int128_t ttop = (__int128_t)u[n] - (__int128_t)carry - (__int128_t)borrow;
    u[n] = (uint64_t)ttop;
    return ttop < 0 ? 1u : 0u;
}

static inline uint64_t add_back_knuth_no128_ref(uint64_t *u, const uint64_t *v, int n) {
    uint64_t carry = 0;
    for (int i = 0; i < n; ++i) {
        carry = add_carry_ref(&u[i], u[i], v[i], carry);
    }
    return add_carry_ref(&u[n], u[n], 0, carry);
}

static inline uint64_t add_back_knuth_int128_ref(uint64_t *u, const uint64_t *v, int n) {
    __uint128_t carry = 0;
    for (int i = 0; i < n; ++i) {
        __uint128_t t = (__uint128_t)u[i] + (__uint128_t)v[i] + carry;
        u[i] = (uint64_t)t;
        carry = t >> LIMB_BITS;
    }
    __uint128_t ttop = (__uint128_t)u[n] + carry;
    u[n] = (uint64_t)ttop;
    return (uint64_t)(ttop >> LIMB_BITS);
}

static inline void mp_add_mod_ref(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    uint64_t carry = mp_add(r, a, b);
    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
}

static inline void mp_div_no128_ref(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den) {
    if (mp_is_zero(den)) {
        mp_zero(q);
        mp_zero(r);
        return;
    }

    if (mp_cmp(num, den) < 0) {
        mp_zero(q);
        mp_copy(r, num);
        return;
    }

    const int m = mp_num_limbs_ref(num);
    const int n = mp_num_limbs_ref(den);

    if (n == 1) {
        mp_uint_t qq = {0};
        uint64_t rem = div_1word_no128_ref(qq, num, den[0]);
        mp_copy(q, qq);
        mp_set(r, rem);
        return;
    }

    mp_uint_t vnorm = {0};
    uint64_t  unorm[MP_N64 + 1] = {0};
    mp_uint_t vraw;
    uint64_t  uraw[MP_N64 + 1] = {0};

    mp_copy(vraw, den);
    for (int i = 0; i < MP_N64; ++i) uraw[i] = num[i];
    uraw[MP_N64] = 0;

    const unsigned s = mp_clz64_ref(vraw[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; ++i) vnorm[i] = vraw[i];
        for (int i = n; i < MP_N64; ++i) vnorm[i] = 0;
        for (int i = 0; i < MP_N64 + 1; ++i) unorm[i] = uraw[i];
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; ++i) {
            uint64_t x = vraw[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
        for (int i = n; i < MP_N64; ++i) vnorm[i] = 0;

        carry = 0;
        for (int i = 0; i < MP_N64 + 1; ++i) {
            uint64_t x = uraw[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
    }

    const int qn = m - n + 1;
    mp_uint_t qlimb = {0};
    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];

    for (int j = qn - 1; j >= 0; --j) {
        uint64_t qhat, rhat;
        rhat = div_2by1_ref(&qhat, unorm[j + n], unorm[j + n - 1], v1);

        for (;;) {
            uint64_t lo, hi;
            mp_mul64_ref(qhat, v2, lo, hi);

            bool too_big = false;
            if (hi > rhat) {
                too_big = true;
            } else if (hi == rhat && lo > unorm[j + n - 2]) {
                too_big = true;
            }

            if (!too_big) break;

            --qhat;
            uint64_t old_rhat = rhat;
            rhat += v1;
            if (rhat < old_rhat) break;
        }

        uint64_t borrow_out = mul_sub_knuth_no128_ref(&unorm[j], vnorm, n, qhat);
        if (borrow_out) {
            add_back_knuth_no128_ref(&unorm[j], vnorm, n);
            --qhat;
        }

        qlimb[j] = qhat;
    }

    mp_uint_t rlimb = {0};
    if (s == 0) {
        for (int i = 0; i < n; ++i) rlimb[i] = unorm[i];
    } else {
        uint64_t carry = 0;
        for (int i = n - 1; i >= 0; --i) {
            uint64_t x = unorm[i];
            rlimb[i] = (x >> s) | carry;
            carry = x << (LIMB_BITS - s);
        }
    }

    mp_copy(q, qlimb);
    mp_copy(r, rlimb);
}

static inline void mp_div_int128_ref(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den) {
    if (mp_is_zero(den)) {
        mp_zero(q);
        mp_zero(r);
        return;
    }

    if (mp_cmp(num, den) < 0) {
        mp_zero(q);
        mp_copy(r, num);
        return;
    }

    const int m = mp_num_limbs_ref(num);
    const int n = mp_num_limbs_ref(den);

    if (n == 1) {
        mp_uint_t qq = {0};
        uint64_t rem = div_1word_int128_ref(qq, num, den[0]);
        mp_copy(q, qq);
        mp_set(r, rem);
        return;
    }

    mp_uint_t vnorm = {0};
    uint64_t  unorm[MP_N64 + 1] = {0};
    mp_uint_t vraw;
    uint64_t  uraw[MP_N64 + 1] = {0};

    mp_copy(vraw, den);
    for (int i = 0; i < MP_N64; ++i) uraw[i] = num[i];
    uraw[MP_N64] = 0;

    const unsigned s = mp_clz64_ref(vraw[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; ++i) vnorm[i] = vraw[i];
        for (int i = n; i < MP_N64; ++i) vnorm[i] = 0;
        for (int i = 0; i < MP_N64 + 1; ++i) unorm[i] = uraw[i];
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; ++i) {
            uint64_t x = vraw[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
        for (int i = n; i < MP_N64; ++i) vnorm[i] = 0;

        carry = 0;
        for (int i = 0; i < MP_N64 + 1; ++i) {
            uint64_t x = uraw[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
    }

    const int qn = m - n + 1;
    mp_uint_t qlimb = {0};
    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];

    for (int j = qn - 1; j >= 0; --j) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + n] << LIMB_BITS) |
                          (__uint128_t)unorm[j + n - 1];
        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << LIMB_BITS) |
                                (__uint128_t)unorm[j + n - 2];
            if (left <= right) break;

            --qhat;
            uint64_t old_rhat = rhat;
            rhat += v1;
            if (rhat < old_rhat) break;
        }

        uint64_t borrow_out = mul_sub_knuth_int128_ref(&unorm[j], vnorm, n, qhat);
        if (borrow_out) {
            add_back_knuth_int128_ref(&unorm[j], vnorm, n);
            --qhat;
        }

        qlimb[j] = qhat;
    }

    mp_uint_t rlimb = {0};
    if (s == 0) {
        for (int i = 0; i < n; ++i) rlimb[i] = unorm[i];
    } else {
        uint64_t carry = 0;
        for (int i = n - 1; i >= 0; --i) {
            uint64_t x = unorm[i];
            rlimb[i] = (x >> s) | carry;
            carry = x << (LIMB_BITS - s);
        }
    }

    mp_copy(q, qlimb);
    mp_copy(r, rlimb);
}

static inline void mp_mul_full_no128_ref(uint64_t *out2n, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < 2 * MP_N64; ++i) out2n[i] = 0;

    for (int i = 0; i < MP_N64; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < MP_N64; ++j) {
            uint64_t lo, hi;
            mp_mul64_ref(a[i], b[j], lo, hi);

            uint64_t x;
            uint64_t c0 = add_carry_ref(&x, out2n[i + j], lo, 0);
            uint64_t c1 = add_carry_ref(&x, x, carry, 0);
            out2n[i + j] = x;
            carry = hi + c0 + c1;
        }

        int k = i + MP_N64;
        while (carry && k < 2 * MP_N64) {
            uint64_t x;
            uint64_t c = add_carry_ref(&x, out2n[k], carry, 0);
            out2n[k] = x;
            carry = c;
            ++k;
        }
    }
}

static inline void mp_mul_full_int128_ref(uint64_t *out2n, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < 2 * MP_N64; ++i) out2n[i] = 0;

    for (int i = 0; i < MP_N64; ++i) {
        __uint128_t carry = 0;
        for (int j = 0; j < MP_N64; ++j) {
            __uint128_t cur = (__uint128_t)a[i] * (__uint128_t)b[j]
                            + (__uint128_t)out2n[i + j]
                            + carry;
            out2n[i + j] = (uint64_t)cur;
            carry = cur >> LIMB_BITS;
        }

        int k = i + MP_N64;
        while (carry && k < 2 * MP_N64) {
            __uint128_t cur = (__uint128_t)out2n[k] + carry;
            out2n[k] = (uint64_t)cur;
            carry = cur >> LIMB_BITS;
            ++k;
        }
    }
}

static inline void mp_mod_2n_n_no128_ref(uint64_t *r, const uint64_t *num2n, const uint64_t *den) {
    const int n = mp_num_limbs_ref(den);

    if (n == 0) {
        mp_zero(r);
        return;
    }

    if (n == 1) {
        uint64_t rem = 0;
        for (int i = 2 * MP_N64 - 1; i >= 0; --i) {
            uint64_t qdummy;
            rem = div_2by1_ref(&qdummy, rem, num2n[i], den[0]);
        }
        mp_set(r, rem);
        return;
    }

    uint64_t vnorm[MP_N64] = {0};
    uint64_t unorm[2 * MP_N64 + 1] = {0};

    const unsigned s = mp_clz64_ref(den[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; ++i) vnorm[i] = den[i];
        for (int i = 0; i < 2 * MP_N64; ++i) unorm[i] = num2n[i];
        unorm[2 * MP_N64] = 0;
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; ++i) {
            uint64_t x = den[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }

        carry = 0;
        for (int i = 0; i < 2 * MP_N64; ++i) {
            uint64_t x = num2n[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
        unorm[2 * MP_N64] = carry;
    }

    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];
    const int qn = 2 * MP_N64 - n + 1;

    for (int j = qn - 1; j >= 0; --j) {
        uint64_t qhat, rhat;
        rhat = div_2by1_ref(&qhat, unorm[j + n], unorm[j + n - 1], v1);

        for (;;) {
            uint64_t lo, hi;
            mp_mul64_ref(qhat, v2, lo, hi);

            bool too_big = false;
            if (hi > rhat) {
                too_big = true;
            } else if (hi == rhat && lo > unorm[j + n - 2]) {
                too_big = true;
            }

            if (!too_big) break;

            --qhat;
            uint64_t old_rhat = rhat;
            rhat += v1;
            if (rhat < old_rhat) break;
        }

        uint64_t borrow_out = mul_sub_knuth_no128_ref(&unorm[j], vnorm, n, qhat);
        if (borrow_out) {
            add_back_knuth_no128_ref(&unorm[j], vnorm, n);
        }
    }

    mp_zero(r);
    if (s == 0) {
        for (int i = 0; i < n; ++i) r[i] = unorm[i];
    } else {
        uint64_t carry = 0;
        for (int i = n - 1; i >= 0; --i) {
            uint64_t x = unorm[i];
            r[i] = (x >> s) | carry;
            carry = x << (LIMB_BITS - s);
        }
    }

    while (mp_cmp(r, den) >= 0) {
        mp_sub(r, r, den);
    }
}

static inline void mp_mod_2n_n_int128_ref(uint64_t *r, const uint64_t *num2n, const uint64_t *den) {
    const int n = mp_num_limbs_ref(den);

    if (n == 0) {
        mp_zero(r);
        return;
    }

    if (n == 1) {
        uint64_t rem = 0;
        for (int i = 2 * MP_N64 - 1; i >= 0; --i) {
            __uint128_t cur = ((__uint128_t)rem << LIMB_BITS) | (__uint128_t)num2n[i];
            rem = (uint64_t)(cur % den[0]);
        }
        mp_set(r, rem);
        return;
    }

    uint64_t vnorm[MP_N64] = {0};
    uint64_t unorm[2 * MP_N64 + 1] = {0};

    const unsigned s = mp_clz64_ref(den[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; ++i) vnorm[i] = den[i];
        for (int i = 0; i < 2 * MP_N64; ++i) unorm[i] = num2n[i];
        unorm[2 * MP_N64] = 0;
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; ++i) {
            uint64_t x = den[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }

        carry = 0;
        for (int i = 0; i < 2 * MP_N64; ++i) {
            uint64_t x = num2n[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (LIMB_BITS - s);
        }
        unorm[2 * MP_N64] = carry;
    }

    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];
    const int qn = 2 * MP_N64 - n + 1;

    for (int j = qn - 1; j >= 0; --j) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + n] << LIMB_BITS) |
                          (__uint128_t)unorm[j + n - 1];

        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << LIMB_BITS) |
                                (__uint128_t)unorm[j + n - 2];
            if (left <= right) break;

            --qhat;
            uint64_t old_rhat = rhat;
            rhat += v1;
            if (rhat < old_rhat) break;
        }

        uint64_t borrow_out = mul_sub_knuth_int128_ref(&unorm[j], vnorm, n, qhat);
        if (borrow_out) {
            add_back_knuth_int128_ref(&unorm[j], vnorm, n);
        }
    }

    mp_zero(r);
    if (s == 0) {
        for (int i = 0; i < n; ++i) r[i] = unorm[i];
    } else {
        uint64_t carry = 0;
        for (int i = n - 1; i >= 0; --i) {
            uint64_t x = unorm[i];
            r[i] = (x >> s) | carry;
            carry = x << (LIMB_BITS - s);
        }
    }

    while (mp_cmp(r, den) >= 0) {
        mp_sub(r, r, den);
    }
}

static inline void mp_mulmod_no128_ref(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    if (mp_is_zero(mod)) {
        mp_zero(r);
        return;
    }

    uint64_t prod[2 * MP_N64];
    mp_mul_full_no128_ref(prod, a, b);
    mp_mod_2n_n_no128_ref(r, prod, mod);
}

static inline void mp_mulmod_int128_ref(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    if (mp_is_zero(mod)) {
        mp_zero(r);
        return;
    }

    uint64_t prod[2 * MP_N64];
    mp_mul_full_int128_ref(prod, a, b);
    mp_mod_2n_n_int128_ref(r, prod, mod);
}

static inline void mp_pow_mod_no128_ref(uint64_t *r, const uint64_t *base, const uint64_t *exp, const uint64_t *mod) {
    mp_uint_t one;
    mp_set(one, 1u);

    if (mp_cmp(mod, one) == 0) {
        mp_zero(r);
        return;
    }

    mp_uint_t bcur;
    if (mp_cmp(base, mod) >= 0) {
        mp_uint_t q, rem;
        mp_div_no128_ref(q, rem, base, mod);
        mp_copy(bcur, rem);
    } else {
        mp_copy(bcur, base);
    }

    int topBit = -1;
    for (int limb = MP_N64 - 1; limb >= 0 && topBit < 0; --limb) {
        uint64_t w = exp[limb];
        if (!w) continue;
        for (int bit = LIMB_BITS - 1; bit >= 0; --bit) {
            if ((w >> bit) & 1u) {
                topBit = limb * LIMB_BITS + bit;
                break;
            }
        }
    }

    if (topBit < 0) {
        if (mp_cmp(one, mod) >= 0) {
            mp_uint_t q, rem;
            mp_div_no128_ref(q, rem, one, mod);
            mp_copy(r, rem);
        } else {
            mp_copy(r, one);
        }
        return;
    }

    mp_uint_t acc;
    mp_copy(acc, bcur);

    for (int i = topBit - 1; i >= 0; --i) {
        mp_uint_t sq;
        mp_mulmod_no128_ref(sq, acc, acc, mod);
        mp_copy(acc, sq);

        const int limb = i / LIMB_BITS;
        const int bit  = i % LIMB_BITS;
        if ((exp[limb] >> bit) & 1u) {
            mp_uint_t tmp;
            mp_mulmod_no128_ref(tmp, acc, bcur, mod);
            mp_copy(acc, tmp);
        }
    }

    mp_copy(r, acc);
}

static inline void mp_pow_mod_int128_ref(uint64_t *r, const uint64_t *base, const uint64_t *exp, const uint64_t *mod) {
    mp_uint_t one;
    mp_set(one, 1u);

    if (mp_cmp(mod, one) == 0) {
        mp_zero(r);
        return;
    }

    mp_uint_t bcur;
    if (mp_cmp(base, mod) >= 0) {
        mp_uint_t q, rem;
        mp_div_int128_ref(q, rem, base, mod);
        mp_copy(bcur, rem);
    } else {
        mp_copy(bcur, base);
    }

    int topBit = -1;
    for (int limb = MP_N64 - 1; limb >= 0 && topBit < 0; --limb) {
        uint64_t w = exp[limb];
        if (!w) continue;
        for (int bit = LIMB_BITS - 1; bit >= 0; --bit) {
            if ((w >> bit) & 1u) {
                topBit = limb * LIMB_BITS + bit;
                break;
            }
        }
    }

    if (topBit < 0) {
        if (mp_cmp(one, mod) >= 0) {
            mp_uint_t q, rem;
            mp_div_int128_ref(q, rem, one, mod);
            mp_copy(r, rem);
        } else {
            mp_copy(r, one);
        }
        return;
    }

    mp_uint_t acc;
    mp_copy(acc, bcur);

    for (int i = topBit - 1; i >= 0; --i) {
        mp_uint_t sq;
        mp_mulmod_int128_ref(sq, acc, acc, mod);
        mp_copy(acc, sq);

        const int limb = i / LIMB_BITS;
        const int bit  = i % LIMB_BITS;
        if ((exp[limb] >> bit) & 1u) {
            mp_uint_t tmp;
            mp_mulmod_int128_ref(tmp, acc, bcur, mod);
            mp_copy(acc, tmp);
        }
    }

    mp_copy(r, acc);
}

static inline void prep_reduced_pair(mp_uint_t a, mp_uint_t b, size_t i) {
    set4(a,
         0x0123456789ABCDEFULL + i,
         0x0FEDCBA987654321ULL ^ (i * 23),
         0xAAAAAAAAAAAAAAAAULL + (i * 29),
         0x5555555555555555ULL ^ (i * 31));

    set4(b,
         0x9E3779B97F4A7C15ULL + i,
         0xBF58476D1CE4E5B9ULL ^ (i << 1),
         0x94D049BB133111EBULL + (i * 3),
         0xD6E8FEB86659FD93ULL ^ (i * 7));

    while (mp_cmp(a, Fq_q.longVal) >= 0) mp_sub(a, a, Fq_q.longVal);
    while (mp_cmp(b, Fq_q.longVal) >= 0) mp_sub(b, b, Fq_q.longVal);
}

TEST(mp_perf_int128_vs_no128, correctness_smoke) {
    uint64_t rng = 0x123456789ABCDEF0ULL;

    for (int i = 0; i < 200; ++i) {
        mp_uint_t a{}, b{}, q1{}, q2{}, r1{}, r2{}, m1{}, m2{}, p1{}, p2{};
        uint64_t prod1[2 * MP_N64], prod2[2 * MP_N64];

        for (int k = 0; k < MP_N64; ++k) {
            a[k] = splitmix64_next(rng);
            b[k] = splitmix64_next(rng);
        }

        while (mp_cmp(a, Fq_q.longVal) >= 0) mp_sub(a, a, Fq_q.longVal);
        while (mp_cmp(b, Fq_q.longVal) >= 0) mp_sub(b, b, Fq_q.longVal);
        if (mp_is_zero(b)) mp_set(b, 17u);

        mp_div_int128_ref(q1, r1, a, b);
        mp_div_no128_ref(q2, r2, a, b);
        ASSERT_EQ(mp_cmp(q1, q2), 0);
        ASSERT_EQ(mp_cmp(r1, r2), 0);

        mp_mul_full_int128_ref(prod1, a, b);
        mp_mul_full_no128_ref(prod2, a, b);
        ASSERT_EQ(std::memcmp(prod1, prod2, sizeof(prod1)), 0);

        mp_mod_2n_n_int128_ref(m1, prod1, Fq_q.longVal);
        mp_mod_2n_n_no128_ref(m2, prod2, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(m1, m2), 0);

        mp_mulmod_int128_ref(m1, a, b, Fq_q.longVal);
        mp_mulmod_no128_ref(m2, a, b, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(m1, m2), 0);

        mp_set(p1, 0u);
        p1[0] = splitmix64_next(rng) | 1ULL;
        mp_pow_mod_int128_ref(m1, a, p1, Fq_q.longVal);
        mp_pow_mod_no128_ref(m2, a, p1, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(m1, m2), 0);
    }
}

TEST(mp_perf_int128_vs_no128, div_1m) {
    constexpr size_t N = 1000000;
    volatile uint64_t sink = 0;
    mp_uint_t a{}, b{}, q{}, r{};

    for (int i = 0; i < 1000; ++i) {
        prep_reduced_pair(a, b, (size_t)i);
        if (mp_is_zero(b)) mp_set(b, 3u);
        mp_div_int128_ref(q, r, a, b);
        sink ^= q[0] ^ r[0];
        mp_div_no128_ref(q, r, a, b);
        sink ^= q[1] ^ r[1];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        if (mp_is_zero(b)) mp_set(b, 3u);
        mp_div_int128_ref(q, r, a, b);
        sink ^= q[0] ^ q[3] ^ r[0];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        if (mp_is_zero(b)) mp_set(b, 3u);
        mp_div_no128_ref(q, r, a, b);
        sink ^= q[1] ^ q[2] ^ r[1];
    }
    uint64_t t3 = now_ns();

    const double ms_i128 = (double)(t1 - t0) / 1e6;
    const double ms_no128 = (double)(t3 - t2) / 1e6;
    std::printf("[perf] mp_div int128 : %.3f ms (N=%zu) => %.2f ns/op\n", ms_i128, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_div no128  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_no128, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] speedup       : %.3fx\n", ms_no128 / ms_i128);
    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(mp_perf_int128_vs_no128, mul_full_1m) {
    constexpr size_t N = 1000000;
    volatile uint64_t sink = 0;
    mp_uint_t a{}, b{};
    uint64_t prod[2 * MP_N64];

    for (int i = 0; i < 1000; ++i) {
        prep_reduced_pair(a, b, (size_t)i);
        mp_mul_full_int128_ref(prod, a, b);
        sink ^= prod[0] ^ prod[7];
        mp_mul_full_no128_ref(prod, a, b);
        sink ^= prod[1] ^ prod[6];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mul_full_int128_ref(prod, a, b);
        sink ^= prod[0] ^ prod[7];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mul_full_no128_ref(prod, a, b);
        sink ^= prod[1] ^ prod[6];
    }
    uint64_t t3 = now_ns();

    const double ms_i128 = (double)(t1 - t0) / 1e6;
    const double ms_no128 = (double)(t3 - t2) / 1e6;
    std::printf("[perf] mp_mul_full int128 : %.3f ms (N=%zu) => %.2f ns/op\n", ms_i128, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_mul_full no128  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_no128, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] speedup            : %.3fx\n", ms_no128 / ms_i128);
    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(mp_perf_int128_vs_no128, mod_2n_n_200k) {
    constexpr size_t N = 200000;
    volatile uint64_t sink = 0;
    mp_uint_t a{}, b{}, r{};
    uint64_t prod[2 * MP_N64];

    for (int i = 0; i < 1000; ++i) {
        prep_reduced_pair(a, b, (size_t)i);
        mp_mul_full_int128_ref(prod, a, b);
        mp_mod_2n_n_int128_ref(r, prod, Fq_q.longVal);
        sink ^= r[0];
        mp_mod_2n_n_no128_ref(r, prod, Fq_q.longVal);
        sink ^= r[1];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mul_full_int128_ref(prod, a, b);
        mp_mod_2n_n_int128_ref(r, prod, Fq_q.longVal);
        sink ^= r[0] ^ r[3];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mul_full_no128_ref(prod, a, b);
        mp_mod_2n_n_no128_ref(r, prod, Fq_q.longVal);
        sink ^= r[1] ^ r[2];
    }
    uint64_t t3 = now_ns();

    const double ms_i128 = (double)(t1 - t0) / 1e6;
    const double ms_no128 = (double)(t3 - t2) / 1e6;
    std::printf("[perf] mp_mod_2n_n int128 : %.3f ms (N=%zu) => %.2f ns/op\n", ms_i128, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_mod_2n_n no128  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_no128, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] speedup            : %.3fx\n", ms_no128 / ms_i128);
    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(mp_perf_int128_vs_no128, mulmod_100k) {
    constexpr size_t N = 100000;
    volatile uint64_t sink = 0;
    mp_uint_t a{}, b{}, r{};

    for (int i = 0; i < 1000; ++i) {
        prep_reduced_pair(a, b, (size_t)i);
        mp_mulmod_int128_ref(r, a, b, Fq_q.longVal);
        sink ^= r[0];
        mp_mulmod_no128_ref(r, a, b, Fq_q.longVal);
        sink ^= r[1];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mulmod_int128_ref(r, a, b, Fq_q.longVal);
        sink ^= r[0] ^ r[3];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        prep_reduced_pair(a, b, i);
        mp_mulmod_no128_ref(r, a, b, Fq_q.longVal);
        sink ^= r[1] ^ r[2];
    }
    uint64_t t3 = now_ns();

    const double ms_i128 = (double)(t1 - t0) / 1e6;
    const double ms_no128 = (double)(t3 - t2) / 1e6;
    std::printf("[perf] mp_mulmod int128 : %.3f ms (N=%zu) => %.2f ns/op\n", ms_i128, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_mulmod no128  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_no128, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] speedup          : %.3fx\n", ms_no128 / ms_i128);
    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(mp_perf_int128_vs_no128, pow_mod_5k) {
    constexpr size_t N = 5000;
    volatile uint64_t sink = 0;
    mp_uint_t a{}, e{}, r{};
    mp_set(e, 0u);
    e[0] = 0xD6E8FEB86659FD93ULL;

    for (int i = 0; i < 50; ++i) {
        mp_uint_t dummy{};
        prep_reduced_pair(a, dummy, (size_t)i);
        mp_pow_mod_int128_ref(r, a, e, Fq_q.longVal);
        sink ^= r[0];
        mp_pow_mod_no128_ref(r, a, e, Fq_q.longVal);
        sink ^= r[1];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_uint_t dummy{};
        prep_reduced_pair(a, dummy, i);
        mp_pow_mod_int128_ref(r, a, e, Fq_q.longVal);
        sink ^= r[0] ^ r[3];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_uint_t dummy{};
        prep_reduced_pair(a, dummy, i);
        mp_pow_mod_no128_ref(r, a, e, Fq_q.longVal);
        sink ^= r[1] ^ r[2];
    }
    uint64_t t3 = now_ns();

    const double ms_i128 = (double)(t1 - t0) / 1e6;
    const double ms_no128 = (double)(t3 - t2) / 1e6;
    std::printf("[perf] mp_pow_mod int128 : %.3f ms (N=%zu) => %.2f ns/op\n", ms_i128, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_pow_mod no128  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_no128, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] speedup           : %.3fx\n", ms_no128 / ms_i128);
    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(fq_perf_inv, rawfq_inv_1m) {
    RawFq Fq;

    constexpr size_t N = 1'000'000;

    RawFq::Element a, inv;

    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) {
        Fq.inv(inv, a);
    }

    {
        RawFq::Element prod;
        RawFq::Element one = Fq.one();
        Fq.set(a, 1);
        Fq.inv(inv, a);
        Fq.mul(prod, a, inv);
        ASSERT_EQ(Fq.toString(prod, 16), Fq.toString(one, 16));

        Fq.set(a, 2);
        Fq.inv(inv, a);
        Fq.mul(prod, a, inv);
        ASSERT_EQ(Fq.toString(prod, 16), Fq.toString(one, 16));
    }

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
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

TEST(fq_perf_inv, rawfq_inv_1m_random) {
    RawFq Fq;
    constexpr size_t N = 1'000'000;

    std::vector<RawFq::Element> inputs(N);

    uint64_t rng = 0x123456789ABCDEF0ULL;

    for (size_t i = 0; i < N; i++) {
        uint64_t t[MP_N64];
        t[0] = splitmix64_next(rng);
        t[1] = splitmix64_next(rng);
        t[2] = splitmix64_next(rng);
        t[3] = splitmix64_next(rng);

        while (mp_cmp(t, Fq_q.longVal) >= 0) {
            mp_sub(t, t, Fq_q.longVal);
        }

        mp_copy(inputs[i].v, t);
        Fq_rawMMul(inputs[i].v, inputs[i].v, Fq_R2.longVal);
    }

    RawFq::Element a, inv;
    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) Fq.inv(inv, a);

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
        Fq.inv(inv, inputs[i]);
        sink ^= inv.v[0];
    }
    uint64_t t1 = now_ns();

    std::printf("[perf] RawFq::inv random: %.3f ms (N=%zu) => %.2f ns/op\n",
                (double)(t1 - t0) / 1e6, N, (double)(t1 - t0) / (double)N);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

TEST(fq_perf_inv, rawfq_inv_1m_random_nonzero) {
    RawFq Fq;
    constexpr size_t N = 1'000'000;

    std::vector<RawFq::Element> inputs(N);

    uint64_t rng = 0x0FEDCBA987654321ULL;

    for (size_t i = 0; i < N; i++) {
        uint64_t t[MP_N64];
        do {
            t[0] = splitmix64_next(rng);
            t[1] = splitmix64_next(rng);
            t[2] = splitmix64_next(rng);
            t[3] = splitmix64_next(rng);

            while (mp_cmp(t, Fq_q.longVal) >= 0) {
                mp_sub(t, t, Fq_q.longVal);
            }
        } while (mp_is_zero(t));

        mp_copy(inputs[i].v, t);
        Fq_rawMMul(inputs[i].v, inputs[i].v, Fq_R2.longVal);
    }

    RawFq::Element a, inv;
    Fq.set(a, 1);
    for (int i = 0; i < 1000; i++) Fq.inv(inv, a);

    volatile uint64_t sink = 0;

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; i++) {
        Fq.inv(inv, inputs[i]);
        sink ^= inv.v[0];
    }
    uint64_t t1 = now_ns();

    std::printf("[perf] RawFq::inv random nonzero: %.3f ms (N=%zu) => %.2f ns/op\n",
                (double)(t1 - t0) / 1e6, N, (double)(t1 - t0) / (double)N);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}

// ============================================================
// local inversion refs for mp_perf.cpp only
// ============================================================

static inline int mp_cmp_ext_ref(const uint64_t *a, const uint64_t *b, int n) {
    for (int i = n - 1; i >= 0; --i) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

static inline bool mp_is_zero_ext_ref(const uint64_t *a, int n) {
    uint64_t acc = 0;
    for (int i = 0; i < n; ++i) acc |= a[i];
    return acc == 0;
}

static inline bool mp_is_one_ext_ref(const uint64_t *a, int n) {
    if (a[0] != 1) return false;
    for (int i = 1; i < n; ++i) {
        if (a[i] != 0) return false;
    }
    return true;
}

static inline bool mp_is_even_ext_ref(const uint64_t *a) {
    return (a[0] & 1ULL) == 0;
}

static inline void mp_copy_ext_ref(uint64_t *dst, const uint64_t *src, int n) {
    std::memcpy(dst, src, sizeof(uint64_t) * n);
}

static inline void mp_zero_ext_ref(uint64_t *a, int n) {
    std::memset(a, 0, sizeof(uint64_t) * n);
}

static inline void mp_shr1_ext_ref(uint64_t *a, int n) {
    uint64_t carry = 0;
    for (int i = n - 1; i >= 0; --i) {
        uint64_t new_carry = a[i] << 63;
        a[i] = (a[i] >> 1) | carry;
        carry = new_carry;
    }
}

static inline void mp_add_ext_ref(uint64_t *r, const uint64_t *a, const uint64_t *b, int n) {
    uint64_t carry = 0;
    for (int i = 0; i < n; ++i) {
        carry = add_carry_ref(&r[i], a[i], b[i], carry);
    }
}

static inline void mp_sub_ext_ref(uint64_t *r, const uint64_t *a, const uint64_t *b, int n) {
    uint64_t borrow = 0;
    for (int i = 0; i < n; ++i) {
        borrow = sub_borrow_ref(&r[i], a[i], b[i], borrow);
    }
}

static inline void mp_reduce_once_ref(uint64_t *a, const uint64_t *mod) {
    while (mp_cmp(a, mod) >= 0) mp_sub(a, a, mod);
}

static inline void mp_reduce_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    mp_copy(r, a);
    while (mp_cmp(r, mod) >= 0) mp_sub(r, r, mod);
}

// ------------------------------------------------------------
// signed bigint as {mag[MP_N64+1], sign}
// sign: -1 / 0 / +1
// ------------------------------------------------------------
struct mp_signed_ref {
    uint64_t mag[MP_N64 + 1];
    int sign;
};

static inline void ms_norm_ref(mp_signed_ref &x) {
    if (mp_is_zero_ext_ref(x.mag, MP_N64 + 1)) x.sign = 0;
}

static inline void ms_zero_ref(mp_signed_ref &x) {
    mp_zero_ext_ref(x.mag, MP_N64 + 1);
    x.sign = 0;
}

static inline void ms_set_u_ref(mp_signed_ref &x, const uint64_t *a) {
    for (int i = 0; i < MP_N64; ++i) x.mag[i] = a[i];
    x.mag[MP_N64] = 0;
    x.sign = mp_is_zero(a) ? 0 : 1;
}

static inline void ms_set_si_ref(mp_signed_ref &x, int64_t v) {
    ms_zero_ref(x);
    if (v == 0) return;

#if defined(__SIZEOF_INT128__)
    uint64_t absv = (v < 0) ? (uint64_t)(-((__int128_t)v)) : (uint64_t)v;
#else
    uint64_t absv = (v < 0) ? (uint64_t)(-(uint64_t)v) : (uint64_t)v;
#endif

    x.mag[0] = absv;
    x.sign = (v < 0) ? -1 : 1;
}

static inline void ms_copy_ref(mp_signed_ref &dst, const mp_signed_ref &src) {
    std::memcpy(dst.mag, src.mag, sizeof(dst.mag));
    dst.sign = src.sign;
}

static inline int ms_mag_cmp_ref(const mp_signed_ref &a, const mp_signed_ref &b) {
    return mp_cmp_ext_ref(a.mag, b.mag, MP_N64 + 1);
}

static inline void ms_add_ref(mp_signed_ref &r, const mp_signed_ref &a, const mp_signed_ref &b) {
    if (a.sign == 0) { ms_copy_ref(r, b); return; }
    if (b.sign == 0) { ms_copy_ref(r, a); return; }

    if (a.sign == b.sign) {
        mp_add_ext_ref(r.mag, a.mag, b.mag, MP_N64 + 1);
        r.sign = a.sign;
        ms_norm_ref(r);
        return;
    }

    int c = ms_mag_cmp_ref(a, b);
    if (c == 0) {
        ms_zero_ref(r);
        return;
    }

    if (c > 0) {
        mp_sub_ext_ref(r.mag, a.mag, b.mag, MP_N64 + 1);
        r.sign = a.sign;
    } else {
        mp_sub_ext_ref(r.mag, b.mag, a.mag, MP_N64 + 1);
        r.sign = b.sign;
    }
    ms_norm_ref(r);
}

static inline void ms_sub_ref(mp_signed_ref &r, const mp_signed_ref &a, const mp_signed_ref &b) {
    mp_signed_ref nb = b;
    nb.sign = -nb.sign;
    ms_add_ref(r, a, nb);
}

static inline void ms_div2_ref(mp_signed_ref &x) {
    if (x.sign == 0) return;
    mp_shr1_ext_ref(x.mag, MP_N64 + 1);
    ms_norm_ref(x);
}

static inline void ms_add_mod_if_odd_then_div2_ref(mp_signed_ref &x, const uint64_t *mod) {
    if ((x.mag[0] & 1ULL) == 0) {
        ms_div2_ref(x);
        return;
    }

    mp_signed_ref m, t;
    ms_zero_ref(m);
    for (int i = 0; i < MP_N64; ++i) m.mag[i] = mod[i];
    m.sign = 1;

    ms_add_ref(t, x, m);
    ms_div2_ref(t);
    ms_copy_ref(x, t);
}

static inline void ms_to_mod_ref(uint64_t *r, const mp_signed_ref &x, const uint64_t *mod) {
    mp_uint_t tmp{};
    for (int i = 0; i < MP_N64; ++i) tmp[i] = x.mag[i];
    mp_reduce_once_ref(tmp, mod);

    if (x.sign < 0 && !mp_is_zero(tmp)) {
        mp_sub(r, mod, tmp);
    } else {
        mp_copy(r, tmp);
    }
}

// ------------------------------------------------------------
// signword packed as x[MP_N64] = sign bit
// 0 => positive / zero
// 1 => negative
// ------------------------------------------------------------
static inline void msw_zero_ref(uint64_t *x) {
    mp_zero_ext_ref(x, MP_N64 + 1);
}

static inline bool msw_is_zero_ref(const uint64_t *x) {
    return mp_is_zero_ext_ref(x, MP_N64);
}

static inline int msw_sign_ref(const uint64_t *x) {
    if (msw_is_zero_ref(x)) return 0;
    return x[MP_N64] ? -1 : 1;
}

static inline void msw_norm_ref(uint64_t *x) {
    if (msw_is_zero_ref(x)) x[MP_N64] = 0;
}

static inline void msw_set_si_ref(uint64_t *x, int64_t v) {
    msw_zero_ref(x);
    if (v == 0) return;

#if defined(__SIZEOF_INT128__)
    uint64_t absv = (v < 0) ? (uint64_t)(-((__int128_t)v)) : (uint64_t)v;
#else
    uint64_t absv = (v < 0) ? (uint64_t)(-(uint64_t)v) : (uint64_t)v;
#endif

    x[0] = absv;
    x[MP_N64] = (v < 0) ? 1 : 0;
}

static inline void msw_copy_ref(uint64_t *dst, const uint64_t *src) {
    mp_copy_ext_ref(dst, src, MP_N64 + 1);
}

static inline int msw_mag_cmp_ref(const uint64_t *a, const uint64_t *b) {
    return mp_cmp_ext_ref(a, b, MP_N64);
}

static inline void msw_add_ref(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    const int sa = msw_sign_ref(a);
    const int sb = msw_sign_ref(b);

    if (sa == 0) { msw_copy_ref(r, b); return; }
    if (sb == 0) { msw_copy_ref(r, a); return; }

    if (sa == sb) {
        mp_add_ext_ref(r, a, b, MP_N64);
        r[MP_N64] = (sa < 0) ? 1 : 0;
        msw_norm_ref(r);
        return;
    }

    int c = msw_mag_cmp_ref(a, b);
    if (c == 0) {
        msw_zero_ref(r);
        return;
    }

    if (c > 0) {
        mp_sub_ext_ref(r, a, b, MP_N64);
        r[MP_N64] = (sa < 0) ? 1 : 0;
    } else {
        mp_sub_ext_ref(r, b, a, MP_N64);
        r[MP_N64] = (sb < 0) ? 1 : 0;
    }
    msw_norm_ref(r);
}

static inline void msw_sub_ref(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t nb[MP_N64 + 1];
    msw_copy_ref(nb, b);
    if (!msw_is_zero_ref(nb)) nb[MP_N64] ^= 1ULL;
    msw_add_ref(r, a, nb);
}

static inline void msw_div2_ref(uint64_t *x) {
    mp_shr1_ext_ref(x, MP_N64);
    msw_norm_ref(x);
}

static inline void msw_add_mod_if_odd_then_div2_ref(uint64_t *x, const uint64_t *mod) {
    if ((x[0] & 1ULL) == 0) {
        msw_div2_ref(x);
        return;
    }

    uint64_t m[MP_N64 + 1];
    uint64_t t[MP_N64 + 1];
    msw_zero_ref(m);
    for (int i = 0; i < MP_N64; ++i) m[i] = mod[i];

    msw_add_ref(t, x, m);
    msw_div2_ref(t);
    msw_copy_ref(x, t);
}

static inline void msw_to_mod_ref(uint64_t *r, const uint64_t *x, const uint64_t *mod) {
    mp_uint_t tmp{};
    for (int i = 0; i < MP_N64; ++i) tmp[i] = x[i];
    mp_reduce_once_ref(tmp, mod);

    if (msw_sign_ref(x) < 0 && !mp_is_zero(tmp)) {
        mp_sub(r, mod, tmp);
    } else {
        mp_copy(r, tmp);
    }
}

// ------------------------------------------------------------
// inversion refs
// ------------------------------------------------------------

static inline bool mp_inv_mod_gmp_signed_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (mp_is_zero(mod) || mp_is_zero(a)) {
        mp_zero(r);
        return false;
    }

    mp_uint_t aa{};
    mp_reduce_ref(aa, a, mod);
    if (mp_is_zero(aa)) {
        mp_zero(r);
        return false;
    }

    uint64_t u[MP_N64 + 1]{};
    uint64_t v[MP_N64 + 1]{};
    for (int i = 0; i < MP_N64; ++i) {
        u[i] = aa[i];
        v[i] = mod[i];
    }

    mp_signed_ref A, B, C, D;
    ms_set_si_ref(A, 1);
    ms_set_si_ref(B, 0);
    ms_set_si_ref(C, 0);
    ms_set_si_ref(D, 1);

    while (!mp_is_one_ext_ref(u, MP_N64 + 1) && !mp_is_one_ext_ref(v, MP_N64 + 1)) {
        while (mp_is_even_ext_ref(u)) {
            mp_shr1_ext_ref(u, MP_N64 + 1);
            ms_add_mod_if_odd_then_div2_ref(A, mod);
            ms_add_mod_if_odd_then_div2_ref(B, mod);
        }

        while (mp_is_even_ext_ref(v)) {
            mp_shr1_ext_ref(v, MP_N64 + 1);
            ms_add_mod_if_odd_then_div2_ref(C, mod);
            ms_add_mod_if_odd_then_div2_ref(D, mod);
        }

        if (mp_cmp_ext_ref(u, v, MP_N64 + 1) >= 0) {
            mp_sub_ext_ref(u, u, v, MP_N64 + 1);
            mp_signed_ref t1, t2;
            ms_sub_ref(t1, A, C);
            ms_sub_ref(t2, B, D);
            ms_copy_ref(A, t1);
            ms_copy_ref(B, t2);
        } else {
            mp_sub_ext_ref(v, v, u, MP_N64 + 1);
            mp_signed_ref t1, t2;
            ms_sub_ref(t1, C, A);
            ms_sub_ref(t2, D, B);
            ms_copy_ref(C, t1);
            ms_copy_ref(D, t2);
        }
    }

    if (mp_is_one_ext_ref(u, MP_N64 + 1)) {
        ms_to_mod_ref(r, A, mod);
        return true;
    }
    if (mp_is_one_ext_ref(v, MP_N64 + 1)) {
        ms_to_mod_ref(r, C, mod);
        return true;
    }

    mp_zero(r);
    return false;
}

static inline bool mp_inv_mod_gmp_signword_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (mp_is_zero(mod) || mp_is_zero(a)) {
        mp_zero(r);
        return false;
    }

    mp_uint_t aa{};
    mp_reduce_ref(aa, a, mod);
    if (mp_is_zero(aa)) {
        mp_zero(r);
        return false;
    }

    uint64_t u[MP_N64 + 1]{};
    uint64_t v[MP_N64 + 1]{};
    for (int i = 0; i < MP_N64; ++i) {
        u[i] = aa[i];
        v[i] = mod[i];
    }

    uint64_t A[MP_N64 + 1], B[MP_N64 + 1], C[MP_N64 + 1], D[MP_N64 + 1];
    msw_set_si_ref(A, 1);
    msw_set_si_ref(B, 0);
    msw_set_si_ref(C, 0);
    msw_set_si_ref(D, 1);

    while (!mp_is_one_ext_ref(u, MP_N64 + 1) && !mp_is_one_ext_ref(v, MP_N64 + 1)) {
        while (mp_is_even_ext_ref(u)) {
            mp_shr1_ext_ref(u, MP_N64 + 1);
            msw_add_mod_if_odd_then_div2_ref(A, mod);
            msw_add_mod_if_odd_then_div2_ref(B, mod);
        }

        while (mp_is_even_ext_ref(v)) {
            mp_shr1_ext_ref(v, MP_N64 + 1);
            msw_add_mod_if_odd_then_div2_ref(C, mod);
            msw_add_mod_if_odd_then_div2_ref(D, mod);
        }

        if (mp_cmp_ext_ref(u, v, MP_N64 + 1) >= 0) {
            mp_sub_ext_ref(u, u, v, MP_N64 + 1);
            uint64_t t1[MP_N64 + 1], t2[MP_N64 + 1];
            msw_sub_ref(t1, A, C);
            msw_sub_ref(t2, B, D);
            msw_copy_ref(A, t1);
            msw_copy_ref(B, t2);
        } else {
            mp_sub_ext_ref(v, v, u, MP_N64 + 1);
            uint64_t t1[MP_N64 + 1], t2[MP_N64 + 1];
            msw_sub_ref(t1, C, A);
            msw_sub_ref(t2, D, B);
            msw_copy_ref(C, t1);
            msw_copy_ref(D, t2);
        }
    }

    if (mp_is_one_ext_ref(u, MP_N64 + 1)) {
        msw_to_mod_ref(r, A, mod);
        return true;
    }
    if (mp_is_one_ext_ref(v, MP_N64 + 1)) {
        msw_to_mod_ref(r, C, mod);
        return true;
    }

    mp_zero(r);
    return false;
}

static inline bool mp_inv_mod_qm2_int128_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (mp_is_zero(mod) || mp_is_zero(a)) {
        mp_zero(r);
        return false;
    }

    mp_uint_t aa{}, e{}, chk{}, one{};
    mp_reduce_ref(aa, a, mod);
    if (mp_is_zero(aa)) {
        mp_zero(r);
        return false;
    }

    mp_copy(e, mod);
    const uint64_t br = mp_sub(e, e, 2u);
    if (br != 0) {
        mp_zero(r);
        return false;
    }

    mp_pow_mod_int128_ref(r, aa, e, mod);
    mp_mulmod_int128_ref(chk, aa, r, mod);
    mp_set(one, 1u);
    return mp_cmp(chk, one) == 0;
}

static inline bool mp_inv_mod_qm2_no128_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (mp_is_zero(mod) || mp_is_zero(a)) {
        mp_zero(r);
        return false;
    }

    mp_uint_t aa{}, e{}, chk{}, one{};
    mp_reduce_ref(aa, a, mod);
    if (mp_is_zero(aa)) {
        mp_zero(r);
        return false;
    }

    mp_copy(e, mod);
    const uint64_t br = mp_sub(e, e, 2u);
    if (br != 0) {
        mp_zero(r);
        return false;
    }

    mp_pow_mod_no128_ref(r, aa, e, mod);
    mp_mulmod_no128_ref(chk, aa, r, mod);
    mp_set(one, 1u);
    return mp_cmp(chk, one) == 0;
}

// ============================================================
// narrow fixed-width xgcd inverse for perf comparison
// GMP-like semantics: compute coefficient, then normalize sign
// ============================================================

static inline uint64_t mp_add_n_ref(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t carry = 0;
    for (int i = 0; i < MP_N64; ++i) {
        carry = add_carry_ref(&r[i], a[i], b[i], carry);
    }
    return carry;
}

static inline uint64_t mp_sub_n_ref(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t borrow = 0;
    for (int i = 0; i < MP_N64; ++i) {
        borrow = sub_borrow_ref(&r[i], a[i], b[i], borrow);
    }
    return borrow;
}

static inline void mp_copy_n_ref(uint64_t *r, const uint64_t *a) {
    std::memcpy(r, a, sizeof(uint64_t) * MP_N64);
}

static inline void mp_zero_n_ref(uint64_t *r) {
    std::memset(r, 0, sizeof(uint64_t) * MP_N64);
}

static inline bool mp_is_one_n_ref(const uint64_t *a) {
    if (a[0] != 1) return false;
    for (int i = 1; i < MP_N64; ++i) {
        if (a[i] != 0) return false;
    }
    return true;
}

static inline void mp_shr1_n_ref(uint64_t *a) {
    uint64_t carry = 0;
    for (int i = MP_N64 - 1; i >= 0; --i) {
        uint64_t next = a[i] << 63;
        a[i] = (a[i] >> 1) | carry;
        carry = next;
    }
}

static inline void mp_div2_mod_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    mp_uint_t t{};
    if ((a[0] & 1ULL) == 0) {
        mp_copy(t, a);
    } else {
        mp_add_n_ref(t, a, mod);   // assume a < mod, result fits fixed width for our perf use
    }
    mp_shr1_n_ref(t);
    mp_copy(r, t);
}

// signed fixed-width value: sign in separate int, magnitude in MP_N64 limbs
// sign: -1 / 0 / +1
static inline void sfix_zero(int &sgn, uint64_t *mag) {
    sgn = 0;
    mp_zero_n_ref(mag);
}

static inline void sfix_set_u(int &sgn, uint64_t *mag, const uint64_t *a) {
    mp_copy_n_ref(mag, a);
    sgn = mp_is_zero(a) ? 0 : 1;
}

static inline void sfix_set_si(int &sgn, uint64_t *mag, int64_t v) {
    mp_zero_n_ref(mag);
    if (v == 0) {
        sgn = 0;
        return;
    }
#if defined(__SIZEOF_INT128__)
    uint64_t absv = (v < 0) ? (uint64_t)(-((__int128_t)v)) : (uint64_t)v;
#else
    uint64_t absv = (v < 0) ? (uint64_t)(-(uint64_t)v) : (uint64_t)v;
#endif
    mag[0] = absv;
    sgn = (v < 0) ? -1 : 1;
}

static inline void sfix_norm(int &sgn, uint64_t *mag) {
    if (mp_is_zero(mag)) sgn = 0;
}

static inline void sfix_copy(int &dsgn, uint64_t *dmag, int ssgn, const uint64_t *smag) {
    dsgn = ssgn;
    mp_copy_n_ref(dmag, smag);
}

static inline void sfix_add(
    int &rsgn, uint64_t *rmag,
    int asgn, const uint64_t *amag,
    int bsgn, const uint64_t *bmag
) {
    if (asgn == 0) {
        sfix_copy(rsgn, rmag, bsgn, bmag);
        return;
    }
    if (bsgn == 0) {
        sfix_copy(rsgn, rmag, asgn, amag);
        return;
    }

    if (asgn == bsgn) {
        mp_add_n_ref(rmag, amag, bmag);
        rsgn = asgn;
        sfix_norm(rsgn, rmag);
        return;
    }

    int c = mp_cmp(amag, bmag);
    if (c == 0) {
        sfix_zero(rsgn, rmag);
        return;
    }

    if (c > 0) {
        mp_sub_n_ref(rmag, amag, bmag);
        rsgn = asgn;
    } else {
        mp_sub_n_ref(rmag, bmag, amag);
        rsgn = bsgn;
    }
    sfix_norm(rsgn, rmag);
}

static inline void sfix_sub(
    int &rsgn, uint64_t *rmag,
    int asgn, const uint64_t *amag,
    int bsgn, const uint64_t *bmag
) {
    sfix_add(rsgn, rmag, asgn, amag, -bsgn, bmag);
}

static inline void sfix_div2_mod(
    int &sgn, uint64_t *mag,
    const uint64_t *mod
) {
    if (sgn == 0) return;

    if ((mag[0] & 1ULL) == 0) {
        mp_shr1_n_ref(mag);
        sfix_norm(sgn, mag);
        return;
    }

    if (sgn > 0) {
        mp_div2_mod_ref(mag, mag, mod);
        sgn = mp_is_zero(mag) ? 0 : 1;
        return;
    }

    // (-x)/2 mod-adapted update for coefficient tracking:
    // if negative odd, do (x + mod)/2 on magnitude and keep sign logic via add/sub semantics
    int tsgn;
    uint64_t tmag[MP_N64], mmag[MP_N64];
    sfix_set_u(tsgn, mmag, mod);
    sfix_add(tsgn, tmag, sgn, mag, +1, mmag);
    // t is now signed; for binary inverse use canonical "add mod then /2"
    // t should be non-negative in the typical path; still handle general case
    if (tsgn < 0) {
        // fallback: r = mod - ((|t|)/2 mod mod)
        mp_uint_t half{};
        mp_copy(half, tmag);
        mp_shr1_n_ref(half);

        if (mp_is_zero(half)) {
            sfix_zero(sgn, mag);
        } else {
            mp_sub_n_ref(mag, mod, half);
            sgn = 1;
        }
    } else {
        mp_shr1_n_ref(tmag);
        sfix_set_u(sgn, mag, tmag);
    }
}

static inline void sfix_to_mod(uint64_t *r, int sgn, const uint64_t *mag, const uint64_t *mod) {
    if (sgn == 0 || mp_is_zero(mag)) {
        mp_zero(r);
        return;
    }

    mp_uint_t t{};
    mp_copy(t, mag);
    while (mp_cmp(t, mod) >= 0) {
        mp_sub(t, t, mod);
    }

    if (sgn < 0) {
        if (mp_is_zero(t)) mp_zero(r);
        else mp_sub_n_ref(r, mod, t);
    } else {
        mp_copy(r, t);
    }
}

static inline bool mp_inv_mod_xgcd_narrow_ref(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (mp_is_zero(a) || mp_is_zero(mod)) {
        mp_zero(r);
        return false;
    }

    mp_uint_t u{}, v{};
    mp_copy(u, a);
    while (mp_cmp(u, mod) >= 0) mp_sub(u, u, mod);
    if (mp_is_zero(u)) {
        mp_zero(r);
        return false;
    }
    mp_copy(v, mod);

    int As = 1, Bs = 0, Cs = 0, Ds = 1;
    mp_uint_t A{}, B{}, C{}, D{}, T{};
    mp_set(A, 1u);
    mp_zero(B);
    mp_zero(C);
    mp_set(D, 1u);

    while (!mp_is_one_n_ref(u) && !mp_is_one_n_ref(v)) {
        while ((u[0] & 1ULL) == 0) {
            mp_shr1_n_ref(u);
            sfix_div2_mod(As, A, mod);
            sfix_div2_mod(Bs, B, mod);
        }

        while ((v[0] & 1ULL) == 0) {
            mp_shr1_n_ref(v);
            sfix_div2_mod(Cs, C, mod);
            sfix_div2_mod(Ds, D, mod);
        }

        if (mp_cmp(u, v) >= 0) {
            mp_sub_n_ref(u, u, v);

            int ts;
            mp_uint_t tm{};
            sfix_sub(ts, tm, As, A, Cs, C);
            As = ts; mp_copy(A, tm);

            sfix_sub(ts, tm, Bs, B, Ds, D);
            Bs = ts; mp_copy(B, tm);
        } else {
            mp_sub_n_ref(v, v, u);

            int ts;
            mp_uint_t tm{};
            sfix_sub(ts, tm, Cs, C, As, A);
            Cs = ts; mp_copy(C, tm);

            sfix_sub(ts, tm, Ds, D, Bs, B);
            Ds = ts; mp_copy(D, tm);
        }
    }

    if (mp_is_one_n_ref(u)) {
        sfix_to_mod(r, As, A, mod);
        return true;
    }
    if (mp_is_one_n_ref(v)) {
        sfix_to_mod(r, Cs, C, mod);
        return true;
    }

    mp_zero(r);
    return false;
}

TEST(mp_perf_inv_all, correctness_smoke) {
    constexpr size_t N = 1000;
    uint64_t rng = 0xD1B54A32D192ED03ULL;

    mp_uint_t one{};
    mp_set(one, 1u);

    for (size_t i = 0; i < N; ++i) {
        mp_uint_t a{}, inv_cur{}, inv_gmp_signed{}, inv_gmp_signword{}, inv_qm2_i128{}, inv_qm2_no128{}, inv_xgcd_narrow{}, chk{};

        do {
            for (int k = 0; k < MP_N64; ++k) {
                a[k] = splitmix64_next(rng);
            }
            while (mp_cmp(a, Fq_q.longVal) >= 0) {
                mp_sub(a, a, Fq_q.longVal);
            }
        } while (mp_is_zero(a));

        ASSERT_TRUE(mp_inv_mod(inv_cur, a, Fq_q.longVal));
        ASSERT_TRUE(mp_inv_mod_gmp_signed_ref(inv_gmp_signed, a, Fq_q.longVal));
        ASSERT_TRUE(mp_inv_mod_gmp_signword_ref(inv_gmp_signword, a, Fq_q.longVal));
        ASSERT_TRUE(mp_inv_mod_qm2_int128_ref(inv_qm2_i128, a, Fq_q.longVal));
        ASSERT_TRUE(mp_inv_mod_qm2_no128_ref(inv_qm2_no128, a, Fq_q.longVal));
        ASSERT_TRUE(mp_inv_mod_xgcd_narrow_ref(inv_xgcd_narrow, a, Fq_q.longVal));

        mp_mulmod_int128_ref(chk, a, inv_cur, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        mp_mulmod_int128_ref(chk, a, inv_gmp_signed, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        mp_mulmod_int128_ref(chk, a, inv_gmp_signword, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        mp_mulmod_int128_ref(chk, a, inv_qm2_i128, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        mp_mulmod_no128_ref(chk, a, inv_qm2_no128, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        mp_mulmod_int128_ref(chk, a, inv_xgcd_narrow, Fq_q.longVal);
        ASSERT_EQ(mp_cmp(chk, one), 0);

        ASSERT_EQ(mp_cmp(inv_cur, inv_gmp_signed), 0);
        ASSERT_EQ(mp_cmp(inv_cur, inv_gmp_signword), 0);
        ASSERT_EQ(mp_cmp(inv_cur, inv_qm2_i128), 0);
        ASSERT_EQ(mp_cmp(inv_cur, inv_qm2_no128), 0);
        ASSERT_EQ(mp_cmp(inv_cur, inv_xgcd_narrow), 0);
    }
}

TEST(mp_perf_inv_all, inv_compare_binary_1m) {
    constexpr size_t N = 10000;

    struct mp_word_vec {
        uint64_t v[MP_N64];
    };

    std::vector<mp_word_vec> inputs(N);
    uint64_t rng = 0x0FEDCBA987654321ULL;

    for (size_t i = 0; i < N; ++i) {
        do {
            for (int k = 0; k < MP_N64; ++k) {
                inputs[i].v[k] = splitmix64_next(rng);
            }
            while (mp_cmp(inputs[i].v, Fq_q.longVal) >= 0) {
                mp_sub(inputs[i].v, inputs[i].v, Fq_q.longVal);
            }
        } while (mp_is_zero(inputs[i].v));
    }

    volatile uint64_t sink = 0;
    mp_uint_t out{};

    for (int i = 0; i < 1000; ++i) {
        mp_inv_mod(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];

        mp_inv_mod_xgcd_narrow_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];

        mp_inv_mod_gmp_signed_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];

        mp_inv_mod_gmp_signword_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];
    }

    uint64_t t0 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_inv_mod(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];
    }
    uint64_t t1 = now_ns();

    uint64_t t2 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_inv_mod_xgcd_narrow_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];
    }
    uint64_t t3 = now_ns();

    uint64_t t4 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_inv_mod_gmp_signed_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];
    }
    uint64_t t5 = now_ns();

    uint64_t t6 = now_ns();
    for (size_t i = 0; i < N; ++i) {
        mp_inv_mod_gmp_signword_ref(out, inputs[i].v, Fq_q.longVal);
        sink ^= out[0];
    }
    uint64_t t7 = now_ns();

    const double ms_cur          = (double)(t1 - t0) / 1e6;
    const double ms_xgcd_narrow  = (double)(t3 - t2) / 1e6;
    const double ms_gmp_signed   = (double)(t5 - t4) / 1e6;
    const double ms_gmp_signword = (double)(t7 - t6) / 1e6;

    std::printf("[perf] mp_inv current           : %.3f ms (N=%zu) => %.2f ns/op\n", ms_cur, N, (double)(t1 - t0) / (double)N);
    std::printf("[perf] mp_inv xgcd_narrow_ref   : %.3f ms (N=%zu) => %.2f ns/op\n", ms_xgcd_narrow, N, (double)(t3 - t2) / (double)N);
    std::printf("[perf] mp_inv gmp_signed_ref    : %.3f ms (N=%zu) => %.2f ns/op\n", ms_gmp_signed, N, (double)(t5 - t4) / (double)N);
    std::printf("[perf] mp_inv gmp_signword_ref  : %.3f ms (N=%zu) => %.2f ns/op\n", ms_gmp_signword, N, (double)(t7 - t6) / (double)N);

    std::printf("[perf] speedup xgcd_narrow  vs current : %.3fx\n", ms_cur / ms_xgcd_narrow);
    std::printf("[perf] speedup gmp_signed   vs current : %.3fx\n", ms_cur / ms_gmp_signed);
    std::printf("[perf] speedup gmp_signword vs current : %.3fx\n", ms_cur / ms_gmp_signword);

    ASSERT_NE((uint64_t)sink, 0xFFFFFFFFFFFFFFFFULL ^ (uint64_t)sink);
}