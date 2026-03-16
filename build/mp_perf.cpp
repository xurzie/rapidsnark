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

