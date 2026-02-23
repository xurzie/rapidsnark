#include "mp.hpp"

#include <cstdint>
#include <string>
#include <climits>
#include <cctype>
#include <cstring>
#include <algorithm>

static inline uint64_t add_carry(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__clang__) || defined(__GNUC__)
    uint64_t t;
    unsigned c1 = __builtin_add_overflow(a, b, &t);
    unsigned c2 = __builtin_add_overflow(t, c, &t);
    *out = t;
    return (c1 | c2);
#else
    uint64_t t0 = a + b;
    uint64_t carry1 = (t0 < a);
    uint64_t t1 = t0 + c;
    uint64_t carry2 = (t1 < t0);
    *out = t1;
    return carry1 | carry2;
#endif
}

static inline uint64_t sub_borrow(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__clang__) || defined(__GNUC__)
    uint64_t t;
    unsigned b1 = __builtin_sub_overflow(a, b, &t);
    unsigned b2 = __builtin_sub_overflow(t, c, &t);
    *out = t;
    return (b1 | b2);
#else
    uint64_t t0 = a - b;
    uint64_t borrow1 = (a < b);
    uint64_t t1 = t0 - c;
    uint64_t borrow2 = (t0 < c);
    *out = t1;
    return borrow1 | borrow2;
#endif
}

static inline void mp_zero(uint64_t *r) {
    std::memset(r, 0, MP_N64 * sizeof(uint64_t));
}

void mp_set(uint64_t *r, uint64_t a) {
    r[0] = a;

    for (int i = 1; i < MP_N64; i++) r[i] = 0;
}

void mp_copy(uint64_t *r, const uint64_t *a) {
    if (r == a) return;

    int i = 0, end = MP_N64, step = 1;

    if (r > a && r < a + MP_N64) {
        i = MP_N64 - 1;
        end = -1;
        step = -1;
    }

    for (; i != end; i += step) r[i] = a[i];
}

int mp_cmp(const uint64_t *a, const uint64_t *b) {
    for (int i = MP_N64 - 1; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

bool mp_is_zero(const uint64_t *a) {
    return (a[0] | a[1] | a[2] | a[3]) == 0;
}

void mp_shl(uint64_t *r, const uint64_t *a, uint64_t k) {
    if (k >= 256u) { mp_zero(r); return; }
    if (k == 0)    { mp_copy(r, a); return; }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63u;

    if (bitShift == 0) {
        for (int i = 3; i >= 0; i--) {
            int si = i - (int)wordShift;
            r[i] = (si >= 0) ? a[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64u - bitShift;
    for (int i = MP_N64 - 1; i >= 0; i--) {
        int si0 = i - (int)wordShift;
        int si1 = si0 - 1;

        uint64_t lo = (si0 >= 0) ? a[si0] : 0;
        uint64_t hi = (si1 >= 0) ? a[si1] : 0;

        r[i] = (lo << bitShift) | (hi >> inv);
    }
}

void mp_shr(uint64_t *r, const uint64_t *a, uint64_t k) {
    if (k >= 256u) { mp_zero(r); return; }
    if (k == 0)    { mp_copy(r, a); return; }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63u;

    if (bitShift == 0) {
        for (uint32_t i = 0; i < MP_N64; i++) {
            uint32_t si = i + wordShift;
            r[i] = (si < 4u) ? a[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64u - bitShift;
    for (size_t i = 0; i < MP_N64; i++) {
        uint32_t si0 = i + wordShift;
        uint32_t si1 = si0 + 1;

        uint64_t lo = (si0 < 4u) ? a[si0] : 0;
        uint64_t hi = (si1 < 4u) ? a[si1] : 0;

        r[i] = (lo >> bitShift) | (hi << inv);
    }
}

uint64_t mp_add(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t c = 0;

    c = add_carry(&r[0], a[0], b[0], c);
    c = add_carry(&r[1], a[1], b[1], c);
    c = add_carry(&r[2], a[2], b[2], c);
    c = add_carry(&r[3], a[3], b[3], c);

    return c;
}

uint64_t mp_add(uint64_t *r, const uint64_t *a, uint64_t b) {
    uint64_t c = 0;

    c = add_carry(&r[0], a[0], b, 0);
    c = add_carry(&r[1], a[1], 0, c);
    c = add_carry(&r[2], a[2], 0, c);
    c = add_carry(&r[3], a[3], 0, c);

    return c;
}

uint64_t mp_sub(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t br = 0;

    br = sub_borrow(&r[0], a[0], b[0], br);
    br = sub_borrow(&r[1], a[1], b[1], br);
    br = sub_borrow(&r[2], a[2], b[2], br);
    br = sub_borrow(&r[3], a[3], b[3], br);

    return br;
}

uint64_t mp_sub(uint64_t *r, const uint64_t *a, uint64_t b) {
    uint64_t br = 0;

    br = sub_borrow(&r[0], a[0], b, 0);
    br = sub_borrow(&r[1], a[1], 0, br);
    br = sub_borrow(&r[2], a[2], 0, br);
    br = sub_borrow(&r[3], a[3], 0, br);

    return br;
}

bool mp_tstbit(const uint64_t *a, size_t bit) {
    if (bit >= 256u) return false;
    return (a[bit >> 6] >> (bit & 63u)) & 1ULL;
}

void mp_export_be(uint8_t *r, const uint64_t *a) {
    const uint8_t *p = (const uint8_t *)a;
    for (size_t i = 0; i < MP_N; i++) {
        r[i] = p[MP_N - 1 - i];
    }
}

void mp_import_be(uint64_t *r, const uint8_t *a) {
    mp_zero(r);

    for (int i = 0; i < MP_N; i++) {
        const uint8_t byte = a[MP_N - 1 - i];
        r[i / 8] |= (uint64_t)byte << (8 * (i % 8));
    }
}

bool mp_set(uint64_t *r, const char *str, uint32_t base) {
    if (base == 0u) base = 10u;
    if (base < 2u || base > 16u) return false;

    mp_zero(r);

    while (*str && std::isspace((unsigned char)*str)) str++;
    if (*str == '+') str++;
    if (*str == '-') return false;

    bool any = false;

#if !defined(__SIZEOF_INT128__)
    return false;
#else
    for (; *str; ++str) {
        if (std::isspace((unsigned char)*str)) break;

        const unsigned char uc = (unsigned char)*str;
        uint32_t digit;
        if (uc >= (unsigned char)'0' && uc <= (unsigned char)'9') {
            digit = (uint32_t)(uc - (unsigned char)'0');
        } else if (uc >= (unsigned char)'a' && uc <= (unsigned char)'f') {
            digit = (uint32_t)(uc - (unsigned char)'a') + 10u;
        } else if (uc >= (unsigned char)'A' && uc <= (unsigned char)'F') {
            digit = (uint32_t)(uc - (unsigned char)'A') + 10u;
        } else {
            return false;
        }
        if (digit >= base) return false;

        if (mp_mul(r, r, base) != 0) return false;
        if (mp_add(r, r, digit) != 0) return false;
        any = true;
    }
#endif

    if (!any) return false;
    while (*str && std::isspace((unsigned char)*str)) str++;
    return *str == '\0';
}

#if defined(__SIZEOF_INT128__)
static inline uint32_t mp_div(uint64_t *q, const uint64_t *a, uint32_t base) {
    uint64_t rem = 0;
    for (int i = 3; i >= 0; i--) {
        __uint128_t cur = (((__uint128_t)rem) << 2*MP_N) | (__uint128_t)a[i];
        q[i] = (uint64_t)(cur / base);
        rem  = (uint64_t)(cur % base);
    }
    return (uint32_t)rem;
}
#endif

std::string mp_get_str(const uint64_t *a, uint32_t base) {
    if (base == 0u) base = 10u;
    if (base < 2u || base > 16u) return std::string();

    if (mp_is_zero(a)) return std::string("0");

#if !defined(__SIZEOF_INT128__)
    return std::string();
#else
    uint64_t v[MP_N64];
    mp_copy(v, a);

    std::string out;

    while (!mp_is_zero(v)) {
        uint64_t q[MP_N64];
        uint32_t rem = mp_div(q, v, base);
        char digit = (rem < 10u) ? (char)('0' + rem) : (char)('a' + (rem - 10u));
        out.push_back(digit);
        mp_copy(v, q);
    }

    std::reverse(out.begin(), out.end());
    return out;
#endif
}

bool mp_fits_int32(const uint64_t *a) {
    if (a[1] || a[2] || a[3]) return false;
    return a[0] <= (uint64_t)INT_MAX;
}

int32_t mp_get_int32(const mp_uint_t a) {
    return (int32_t)a[0];
}

static inline void mp_add_mod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    uint64_t carry = mp_add(r, a, b);

    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
      mp_uint_t aa, bb, q, rem;

    if (mp_cmp(a, mod) >= 0) {
        mp_div(q, rem, a, mod);
        mp_copy(aa, rem);
    } else {
        mp_copy(aa, a);
    }

    if (mp_cmp(b, mod) >= 0) {
        mp_div(q, rem, b, mod);
        mp_copy(bb, rem);
    } else {
        mp_copy(bb, b);
    }

    const uint64_t car = mp_add(r, aa, bb);

    if (car || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
        if (mp_cmp(r, mod) >= 0) mp_sub(r, r, mod);
    }
}

static inline void mp_mul_mod(uint64_t *r, const uint64_t *a, uint32_t b, const uint64_t *mod) {
    uint64_t res[MP_N64]; mp_zero(res);

    uint64_t cur[MP_N64];
    if (mp_cmp(a, mod) >= 0) {
        mp_uint_t q, rem;
        mp_div(q, rem, a, mod);
        mp_copy(cur, rem);
    } else {
        mp_copy(cur, a);
    }

    uint32_t k = b;
    while (k) {
        if (k & 1u) {
            uint64_t tmp[MP_N64];
            mp_add_mod(tmp, res, cur, mod);
            mp_copy(res, tmp);
        }
        k >>= 1u;
        if (k) {
            uint64_t tmp[MP_N64];
            mp_add_mod(tmp, cur, cur, mod);
            mp_copy(cur, tmp);
        }
    }
    mp_copy(r, res);
}

bool mp_set_mod(uint64_t *r, const char *str, uint32_t base, const uint64_t *mod) {
    mp_zero(r);

    if (base == 0u) base = 10u;
    if (base < 2u || base > 16u) return false;

    while (*str && std::isspace((unsigned char)*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    uint64_t acc[MP_N64];
    mp_zero(acc);

    bool any = false;

#if !defined(__SIZEOF_INT128__)
    return false;
#else
    auto reduce_small = [&](uint64_t *x) {
        // x < mod*base + 15  => at most base subtractions (base <= 16)
        while (mp_cmp(x, mod) >= 0) {
            mp_sub(x, x, mod);
        }
    };

    for (; *str; str++) {
        if (std::isspace((unsigned char)*str)) break;

        const unsigned char uc = (unsigned char)*str;
        uint32_t d;
        if (uc >= (unsigned char)'0' && uc <= (unsigned char)'9') {
            d = (uint32_t)(uc - (unsigned char)'0');
        } else if (uc >= (unsigned char)'a' && uc <= (unsigned char)'f') {
            d = (uint32_t)(uc - (unsigned char)'a') + 10u;
        } else if (uc >= (unsigned char)'A' && uc <= (unsigned char)'F') {
            d = (uint32_t)(uc - (unsigned char)'A') + 10u;
        } else {
            return false;
        }

        // acc = (acc * base + d) mod mod
        uint64_t tmp[MP_N64];
        {
            __uint128_t carry = 0;
            for (size_t i = 0; i < MP_N64; i++) {
                __uint128_t prod = (__uint128_t)acc[i] * (__uint128_t)base + carry;
                tmp[i] = (uint64_t)prod;
                carry = prod >> 2*MP_N;
            }
            // safety: if overflow beyond 256-bit, fall back to the slow-but-correct path for this digit
            if (carry) {
                uint64_t t1[MP_N64];
                mp_mul_mod(t1, acc, base, mod);
                uint64_t dv[MP_N64];
                mp_set(dv, d);
                uint64_t t2[MP_N64];
                mp_add_mod(t2, t1, dv, mod);
                mp_copy(acc, t2);
                any = true;
                continue;
            }
        }

        // tmp += d
        {
            uint64_t c = d;
            for (size_t i = 0; i < MP_N64 && c; i++) {
                uint64_t before = tmp[i];
                tmp[i] += c;
                c = (tmp[i] < before) ? 1u : 0u;
            }
            if (c) {
                // overflow beyond 256-bit => fallback for this digit
                uint64_t t1[MP_N64];
                mp_mul_mod(t1, acc, base, mod);
                uint64_t dv[MP_N64];
                mp_set(dv, d);
                uint64_t t2[MP_N64];
                mp_add_mod(t2, t1, dv, mod);
                mp_copy(acc, t2);
                any = true;
                continue;
            }
        }

        mp_copy(acc, tmp);
        reduce_small(acc);
        any = true;
    }

    if (!any) return false;

    while (*str && std::isspace((unsigned char)*str)) str++;
    if (*str != '\0') return false;

    if (neg && !mp_is_zero(acc)) {
        uint64_t tmp[MP_N64];
        mp_sub(tmp, mod, acc);
        mp_copy(acc, tmp);
    }

    mp_copy(r, acc);
    return true;
#endif
}

static inline unsigned mp_clz64(uint64_t x) { return x ? __builtin_clzll(x) : 64; }

static inline int mp_num_limbs(const uint64_t *a) {
    if (a[3]) return 4;
    if (a[2]) return 3;
    if (a[1]) return 2;
    if (a[0]) return 1;

    return 0;
}

#if defined(__SIZEOF_INT128__)
static inline uint64_t div_1word(uint64_t *q, const uint64_t *u, uint64_t v) {
    __uint128_t rem = 0;

    for (int i = MP_N64 - 1; i >= 0; i--) {
        __uint128_t cur = (rem << 2*MP_N) | (__uint128_t)u[i];
        q[i] = (uint64_t)(cur / v);
        rem  = (uint64_t)(cur % v);
    }

    return (uint64_t)rem;
}

static inline uint64_t mul_sub_knuth(uint64_t *u, const uint64_t *v, int n, uint64_t qhat) {
    __uint128_t carry = 0;
    uint64_t borrow = 0;

    for (int i = 0; i < n; i++) {
        __uint128_t prod = (__uint128_t)qhat * (__uint128_t)v[i] + carry;
        uint64_t pl = (uint64_t)prod;
        carry = (prod >> 64);
        __int128_t t = (__int128_t)u[i] - (__int128_t)pl - (__int128_t)borrow;
        u[i] = (uint64_t)t;
        borrow = (t < 0) ? 1u : 0u;
    }

    __int128_t ttop = (__int128_t)u[n] - (__int128_t)carry - (__int128_t)borrow;
    u[n] = (uint64_t)ttop;

    return (ttop < 0) ? 1u : 0u;
}

static inline uint64_t add_back_knuth(uint64_t *u, const uint64_t *v, int n) {
    __uint128_t carry = 0;

    for (int i = 0; i < n; i++) {
        __uint128_t t = (__uint128_t)u[i] + (__uint128_t)v[i] + carry;
        u[i] = (uint64_t)t;
        carry = (t >> 2*MP_N);
    }

    __uint128_t ttop = (__uint128_t)u[n] + carry;
    u[n] = (uint64_t)ttop;

    return (uint64_t)(ttop >> 2*MP_N);
}
#endif

void mp_div(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den) {

    if (mp_cmp(num, den) < 0) {
        mp_zero(q);
        mp_copy(r, num);
        return;
    }

#if !defined(__SIZEOF_INT128__)
    mp_zero(q);
    mp_zero(r);
    return;
#else
    int m = mp_num_limbs(num);
    int n = mp_num_limbs(den);

    if (n == 1) {
        uint64_t v = den[0];
        mp_uint_t u  = { num[0], num[1], num[2], num[3] };
        mp_uint_t qq = {0};

        uint64_t rem = div_1word(qq, u, v);

        mp_copy(q, qq);
        mp_set(r, rem);
        return;
    }

    mp_uint_t vnorm             = {0};
    uint64_t  unorm[MP_N64 + 1] = {0};
    mp_uint_t vraw              = { den[0], den[1], den[2], den[3] };
    uint64_t  uraw[MP_N64 + 1]  = { num[0], num[1], num[2], num[3], 0 };

    unsigned s = mp_clz64(vraw[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; i++) vnorm[i] = vraw[i];
        for (int i = n; i < MP_N64; i++) vnorm[i] = 0;
        for (int i = 0; i < MP_N64 + 1; i++) unorm[i] = uraw[i];
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; i++) {
            uint64_t x = vraw[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (2*MP_N - s);
        }

        for (int i = n; i < MP_N64; i++) vnorm[i] = 0;

        carry = 0;

        for (int i = 0; i < MP_N64 + 1; i++) {
            uint64_t x = uraw[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (2*MP_N - s);
        }
    }

    int qn = m - n + 1;
    mp_uint_t qlimb = {0};

    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];

    for (int j = qn - 1; j >= 0; j--) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + n] << 2*MP_N) | (__uint128_t)unorm[j + n - 1];
        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 2*MP_N) | (__uint128_t)unorm[j + n - 2];
            if (left <= right) break;
            qhat--;
            rhat += v1;
            if (rhat < v1) break;
        }

        uint64_t *u_seg = &unorm[j];
        uint64_t borrow_out = mul_sub_knuth(u_seg, vnorm, n, qhat);

        if (borrow_out) {
            add_back_knuth(u_seg, vnorm, n);
            qhat--;
        }

        qlimb[j] = qhat;
    }

    mp_uint_t rlimb = {0};

    if (s == 0) {
        for (int i = 0; i < n; i++) rlimb[i] = unorm[i];
    } else {
        uint64_t carry = 0;

        for (int i = n - 1; i >= 0; i--) {
            uint64_t x = unorm[i];
            rlimb[i] = (x >> s) | carry;
            carry = x << (2*MP_N - s);
        }
    }

    mp_copy(q, qlimb);
    mp_copy(r, rlimb);
#endif
}

#if defined(__SIZEOF_INT128__)

static inline void mp_mul_full(uint64_t *out2n,
                               const uint64_t *a,
                               const uint64_t *b)
{
    for (int i = 0; i < 2 * MP_N64; i++) out2n[i] = 0;

    for (int i = 0; i < MP_N64; i++) {
        __uint128_t carry = 0;
        for (int j = 0; j < MP_N64; j++) {
            __uint128_t cur = (__uint128_t)a[i] * (__uint128_t)b[j]
                            + (__uint128_t)out2n[i + j]
                            + carry;
            out2n[i + j] = (uint64_t)cur;
            carry = cur >> 2*MP_N;
        }

        int k = i + MP_N64;
        while (carry) {
            __uint128_t cur = (__uint128_t)out2n[k] + carry;
            out2n[k] = (uint64_t)cur;
            carry = cur >> 2*MP_N;
            k++;
        }
    }
}

static inline void mp_mod_2n_n(uint64_t *r,
                              const uint64_t *num2n,
                              const uint64_t *den)
{
    if (den[MP_N64 - 1] == 0) { mp_zero(r); return; }

    uint64_t vnorm[MP_N64];
    uint64_t unorm[2 * MP_N64 + 1];

    unsigned s = mp_clz64(den[MP_N64 - 1]);

    if (s == 0) {
        for (int i = 0; i < MP_N64; i++) vnorm[i] = den[i];
        for (int i = 0; i < 2 * MP_N64; i++) unorm[i] = num2n[i];
        unorm[2 * MP_N64] = 0;
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < MP_N64; i++) {
            uint64_t x = den[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (2*MP_N - s);
        }

        carry = 0;
        for (int i = 0; i < 2 * MP_N64; i++) {
            uint64_t x = num2n[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (2*MP_N - s);
        }
        unorm[2 * MP_N64] = carry;
    }

    const uint64_t v1 = vnorm[MP_N64 - 1];
    const uint64_t v2 = vnorm[MP_N64 - 2];

    for (int j = MP_N64; j >= 0; j--) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + MP_N64] << 2*MP_N)
                        | (__uint128_t)unorm[j + MP_N64 - 1];

        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 2*MP_N) | (__uint128_t)unorm[j + MP_N64 - 2];
            if (left <= right) break;
            qhat--;
            rhat += v1;
            if (rhat < v1) break;
        }

        uint64_t borrow_out = mul_sub_knuth(&unorm[j], vnorm, MP_N64, qhat);

        if (borrow_out) {
            add_back_knuth(&unorm[j], vnorm, MP_N64);
        }
    }

    if (s == 0) {
        for (int i = 0; i < MP_N64; i++) r[i] = unorm[i];
    } else {
        uint64_t carry = 0;

        for (int i = MP_N64 - 1; i >= 0; i--) {
            uint64_t x = unorm[i];
            r[i] = (x >> s) | carry;
            carry = x << (2*MP_N - s);
        }
    }

    if (mp_cmp(r, den) >= 0) {
        mp_sub(r, r, den);
    }
}

static inline void mp_mulmod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod)
{
    if (mod[MP_N64 - 1] == 0) {
        // If modulus fits in 64-bit, do 64-bit reduction
        if (mod[1] == 0 && mod[2] == 0 && mod[3] == 0) {
            const uint64_t m = mod[0];
            if (m == 0) { mp_zero(r); return; }

            auto red_u64 = [&](const uint64_t *x) -> uint64_t {
                __uint128_t rem = 0;
                for (int i = MP_N64 - 1; i >= 0; --i) {
                    rem = ((rem << 2*MP_N) | x[i]) % m;
                }
                return (uint64_t)rem;
            };

            const uint64_t ra = red_u64(a);
            const uint64_t rb = red_u64(b);
            const uint64_t out = (uint64_t)(((__uint128_t)ra * (__uint128_t)rb) % m);
            mp_set(r, out);
            return;
        }

        mp_uint_t qa, ra, qb, rb;
        mp_div(qa, ra, a, mod);
        mp_div(qb, rb, b, mod);
        uint64_t prod[2 * MP_N64];
        mp_mul_full(prod, ra, rb);
        mp_mod_2n_n(r, prod, mod); // ra/rb < mod
        return;
    }

    uint64_t prod[2 * MP_N64];
    mp_mul_full(prod, a, b);
    mp_mod_2n_n(r, prod, mod);
}

void mp_pow_mod(uint64_t *r, const uint64_t *base, const uint64_t *exp, const uint64_t *mod)
{
    mp_uint_t one;
    mp_set(one, 1u);

    if (mp_cmp(mod, one) == 0) { mp_zero(r); return; }

    mp_uint_t bcur;

    if (mp_cmp(base, mod) >= 0) {
        mp_uint_t q, rem;
        mp_div(q, rem, base, mod);
        mp_copy(bcur, rem);
    } else {
        mp_copy(bcur, base);
    }

    int topBit = -1;
    for (int limb = MP_N64 - 1; limb >= 0 && topBit < 0; --limb) {
        uint64_t w = exp[limb];
        if (!w) continue;

        for (int bit = 2*MP_N - 1; bit >= 0; --bit) {
            if ((w >> bit) & 1u) {
                topBit = limb * 2*MP_N + bit;
                break;
            }
        }
    }

    if (topBit < 0) {
        if (mp_cmp(one, mod) >= 0) {
            mp_uint_t q, rem;
            mp_div(q, rem, one, mod);
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
        mp_mulmod(sq, acc, acc, mod);
        mp_copy(acc, sq);

        const int limb = i >> 6;
        const int bit  = i & 2*MP_N - 1;
        if ((exp[limb] >> bit) & 1u) {
            mp_uint_t tmp;
            mp_mulmod(tmp, acc, bcur, mod);
            mp_copy(acc, tmp);
        }
    }

    mp_copy(r, acc);
}

#endif // __SIZEOF_INT128__

void mp_set_mod(uint64_t *r, int64_t a, const uint64_t *mod) {
    if (a >= 0) {
        mp_set(r, (uint64_t)a);

        if (mp_cmp(r, mod) >= 0) {
            uint64_t q[MP_N64], rem[MP_N64];

             mp_div(q, rem, r, mod);
             mp_copy(r, rem);
        }

        return;
    }

#if defined(__SIZEOF_INT128__)
    uint64_t absv = (uint64_t)(-((__int128_t)a));
#else
    uint64_t absv = (uint64_t)(-(uint64_t)a);
#endif

    uint64_t av[MP_N64]; mp_set(av, absv);

    if (mp_is_zero(av)) { mp_zero(r); return; }

    if (mp_cmp(av, mod) >= 0) {
        uint64_t q[MP_N64], rem[MP_N64];
        mp_div(q, rem, r, mod);
        mp_copy(r, rem);
    }

    if (mp_is_zero(av)) { mp_zero(r); return; }

    mp_sub(r, mod, av);
}

static inline void mp_mul64(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
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

uint64_t mp_mul(uint64_t *r, const uint64_t *a, uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 2*MP_N;
    }
    return (uint64_t)carry;
#else
    uint64_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        uint64_t lo, hi;
        mp_mul64(a[i], b, lo, hi);

        uint64_t out = lo + carry;
        hi += (out < lo);

        r[i] = out;
        carry = hi;
    }
    return carry;
#endif
}

uint64_t mp_addmul(uint64_t *r, const uint64_t *a, uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)r[i] + (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 2*MP_N;
    }
    return (uint64_t)carry;
#else
    uint64_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        uint64_t lo, hi;
        mp_mul64(a[i], b, lo, hi);

        // r[i] + lo + carry
        uint64_t t = r[i] + lo;
        uint64_t c1 = (t < r[i]);
        uint64_t out = t + carry;
        uint64_t c2 = (out < t);

        r[i] = out;

        // carry = hi + c1 + c2
        uint64_t c = hi;
        c += (uint64_t)c1;
        c += (uint64_t)c2;
        carry = c;
    }
    return carry;
#endif
}

void mp_and(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] & b[i];
}
void mp_or(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] | b[i];
}
void mp_xor(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] ^ b[i];
}
void mp_not(uint64_t *r, const uint64_t *a) {
    for (int i = 0; i < MP_N64; i++) r[i] = ~a[i];
}

uint64_t mp_add(uint64_t *r, const uint64_t *a, size_t an, const uint64_t *b, size_t bn) {
    size_t n = (an > bn) ? an : bn;
    uint64_t carry = 0;

    for (size_t i = 0; i < n; i++) {
        uint64_t ai = (i < an) ? a[i] : 0;
        uint64_t bi = (i < bn) ? b[i] : 0;
        carry = add_carry(&r[i], ai, bi, carry);
    }

    return carry;
}

static inline int mp_is_even(const uint64_t *a) {
    return ((a[0] & 1ULL) == 0ULL);
}

static inline void mp_sub_mod(uint64_t *x, const uint64_t *y, const uint64_t *mod) {
    // x = x - y (mod mod), assumes 0 <= x,y < mod and mod odd
    uint64_t t0, t1, t2, t3;
    uint64_t br = 0;

    br = sub_borrow(&t0, x[0], y[0], br);
    br = sub_borrow(&t1, x[1], y[1], br);
    br = sub_borrow(&t2, x[2], y[2], br);
    br = sub_borrow(&t3, x[3], y[3], br);

    // if borrow, add mod back
    if (br) {
        uint64_t c = 0;
        c = add_carry(&t0, t0, mod[0], c);
        c = add_carry(&t1, t1, mod[1], c);
        c = add_carry(&t2, t2, mod[2], c);
        c = add_carry(&t3, t3, mod[3], c);
        (void)c;
    }

    x[0] = t0; x[1] = t1; x[2] = t2; x[3] = t3;
}

// x = x/2 mod mod, mod must be odd, 0<=x<mod
static inline void mp_div2_mod(uint64_t *x, const uint64_t *mod) {
    // if x is odd: x += mod  (still fits into 5 limbs, but we discard top carry after >>1)
    if (x[0] & 1u) {
        uint64_t c = 0;
        c = add_carry(&x[0], x[0], mod[0], c);
        c = add_carry(&x[1], x[1], mod[1], c);
        c = add_carry(&x[2], x[2], mod[2], c);
        c = add_carry(&x[3], x[3], mod[3], c);
        (void)c;
    }
    // x >>= 1
    uint64_t b3 = x[3];
    uint64_t b2 = x[2];
    uint64_t b1 = x[1];
    uint64_t b0 = x[0];
    x[0] = (b0 >> 1) | (b1 << (2*MP_N - 1));
    x[1] = (b1 >> 1) | (b2 << (2*MP_N - 1));
    x[2] = (b2 >> 1) | (b3 << (2*MP_N - 1));
    x[3] = (b3 >> 1);
}

static inline bool mp_is_one(const uint64_t *a) {
    return a[0] == 1u && a[1] == 0u && a[2] == 0u && a[3] == 0u;
}

bool mp_inv_mod(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    uint64_t u[MP_N64], v[MP_N64];
    mp_copy(u, a);
    mp_copy(v, mod);

    if (mp_is_zero(u)) { mp_zero(r); return false; }

    uint64_t x1[MP_N64], x2[MP_N64];
    mp_set(x1, 1u);
    mp_zero(x2);

    while (!mp_is_one(u) && !mp_is_one(v)) {

        while (mp_is_even(u)) {
            mp_shr(u, u, 1u);
            mp_div2_mod(x1, mod);
        }

        while (mp_is_even(v)) {
            mp_shr(v, v, 1u);
            mp_div2_mod(x2, mod);
        }

        if (mp_cmp(u, v) >= 0) {
            mp_sub(u, u, v);
            mp_sub_mod(x1, x2, mod);   // x1 = x1 - x2 (mod)
        } else {
            mp_sub(v, v, u);
            mp_sub_mod(x2, x1, mod);   // x2 = x2 - x1 (mod)
        }
    }

    if (mp_is_one(u)) mp_copy(r, x1);
    else              mp_copy(r, x2);
    return true;
}