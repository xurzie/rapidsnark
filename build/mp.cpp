#include "mp.hpp"

#include <cstdint>
#include <string>
#include <climits>
#include <cctype>
#include <cstring>

static inline uint64_t add_carry(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__SIZEOF_INT128__)
    __uint128_t t = ( (__uint128_t)a ) + b + c;
    *out = (uint64_t)t;
    return (uint64_t)(t >> 64);
#else
    uint64_t t = a + b;
    uint64_t carry1 = (t < a);
    uint64_t u = t + c;
    uint64_t carry2 = (u < t);
    *out = u;
    return carry1 | carry2;
#endif
}

static inline uint64_t sub_borrow(uint64_t *out, uint64_t a, uint64_t b, uint64_t c) {
#if defined(__SIZEOF_INT128__)
    __int128_t t = ( (__int128_t)a ) - (__int128_t)b - (__int128_t)c;
    *out = (uint64_t)t;
    return (t < 0) ? 1u : 0u;
#else
    uint64_t t = a - b;
    uint64_t borrow1 = (a < b);
    uint64_t u = t - c;
    uint64_t borrow2 = (t < c);
    *out = u;
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

    if (wordShift >= 4u) { mp_zero(r); return; }

    if (bitShift == 0) {
        for (int i = 3; i >= 0; i--) {
            int si = i - (int)wordShift;
            r[i] = (si >= 0) ? a[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64u - bitShift;
    for (int i = 3; i >= 0; i--) {
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

    if (wordShift >= 4u) { mp_zero(r); return; }

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

int mp_tstbit(const uint64_t *a, size_t bit) {
    if (bit >= 256u) return 0;
    const size_t w = bit >> 6;
    const size_t s = bit & 63u;
    return (int)((a[w] >> s) & 1ULL);
}

void mp_export_be(uint8_t *r, const uint64_t *a) {
    for (int i = 0; i < MP_N; i++) {
        uint8_t byte = (uint8_t)((a[i / 8] >> (8 * (i % 8))) & 0xFFu);
        r[MP_N - 1 - i] = byte;
    }
}

void mp_import_be(uint64_t *r, const uint8_t *a) {
    mp_zero(r);

    for (int i = 0; i < MP_N; i++) {
        const uint8_t byte = a[MP_N - 1 - i];
        r[i / 8] |= (uint64_t)byte << (8 * (i % 8));
    }
}

#if defined(__SIZEOF_INT128__)
static inline uint64_t mp_mul_small_dec(uint64_t *x, uint32_t m) {
    __uint128_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)x[i] * (uint64_t)m + carry;
        x[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t mp_add_small_dec(uint64_t *x, uint32_t add) {
    __uint128_t t = (__uint128_t)x[0] + (uint64_t)add;
    x[0] = (uint64_t)t;
    uint64_t c = (uint64_t)(t >> 64);

    for (int i = 1; i < MP_N64 && c; i++) {
        __uint128_t ti = (__uint128_t)x[i] + c;
        x[i] = (uint64_t)ti;
        c = (uint64_t)(ti >> 64);
    }
    return c;
}
#endif

bool mp_set(uint64_t *r, const char *str, uint32_t base) {
    if (!r || !str) return false;
    if (base == 0u) base = 10u;
    if (base != 10u) return false;

    mp_zero(r);

    while (*str && isspace((unsigned char)*str)) str++;
    if (*str == '+') str++;
    if (*str == '-') return false;

    bool any = false;
    for (; *str; str++) {
        if (*str < '0' || *str > '9') break;
        any = true;

#if !defined(__SIZEOF_INT128__)
        return false;
#else
        uint32_t digit = (uint32_t)(*str - '0');
        if (mp_mul_small_dec(r, 10u) != 0) return false;
        if (mp_add_small_dec(r, digit) != 0) return false;
#endif
    }

    if (!any) return false;

    while (*str && isspace((unsigned char)*str)) str++;
    if (*str != '\0') return false;

    return true;
}

#if defined(__SIZEOF_INT128__)
static inline uint32_t mp_div(uint64_t *q, const uint64_t *a, uint32_t base) {
    uint64_t rem = 0;
    for (int i = 3; i >= 0; i--) {
        __uint128_t cur = (((__uint128_t)rem) << 64) | (__uint128_t)a[i];
        q[i] = (uint64_t)(cur / base);
        rem  = (uint64_t)(cur % base);
    }
    return (uint32_t)rem;
}
#endif

std::string mp_set_str(const uint64_t *a, uint32_t base) {
    if (!a) return std::string();
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

    for (size_t i = 0, j = out.size() - 1; i < j; i++, j--) {
        char t = out[i]; out[i] = out[j]; out[j] = t;
    }
    return out;
#endif
}

bool mp_fits_int32(const uint64_t *a) {
    if (!a) return false;
    if (a[1] || a[2] || a[3]) return false;
    return a[0] <= (uint64_t)INT_MAX;
}

int32_t mp_get_int32(const mp_uint_t a) {
    return (int32_t)a[0];
}

void mp_add_mod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    uint64_t carry = mp_add(r, a, b);

    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
}

void mp_mul_mod(uint64_t *r, const uint64_t *a, uint32_t b, const uint64_t *mod) {
    uint64_t res[MP_N64]; mp_set(res, 0);
    uint64_t cur[MP_N64]; mp_copy(cur, a);

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
    if (!r || !str || !mod) return false;

    mp_set(r, 0u);

    if (base == 0u) base = 10u;
    if (base < 2u || base > 16u) return false;

    while (*str && std::isspace((unsigned char)*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    uint64_t acc[MP_N64];
    mp_set(acc, 0u);

    bool any = false;

    for (; *str; str++) {
        if (std::isspace((unsigned char)*str)) break;

        unsigned char uc = (unsigned char)*str;
        if (!std::isxdigit(uc)) return false;

        uint32_t d;
        if (std::isdigit(uc)) {
            d = (uint32_t)(uc - (unsigned char)'0');
        } else {
            d = (uint32_t)(std::tolower(uc) - (unsigned char)'a') + 10u;
        }

        if (d >= base) return false;

        uint64_t t1[MP_N64];
        mp_mul_mod(t1, acc, base, mod);

        uint64_t dv[MP_N64];
        mp_set(dv, (uint64_t)d);

        uint64_t t2[MP_N64];
        mp_add_mod(t2, t1, dv, mod);

        mp_copy(acc, t2);
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
static inline uint64_t div_1word(uint64_t *q, const uint64_t *u, int m, uint64_t v) {
    __uint128_t rem = 0;

    for (int i = m - 1; i >= 0; i--) {
        __uint128_t cur = (rem << 64) | (__uint128_t)u[i];
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
        carry = (t >> 64);
    }

    __uint128_t ttop = (__uint128_t)u[n] + carry;
    u[n] = (uint64_t)ttop;

    return (uint64_t)(ttop >> 64);
}
#endif

bool mp_div(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den) {
    if (!q || !r || !num || !den) return false;
    if (mp_is_zero(den)) return false;

    if (mp_cmp(num, den) < 0) {
        mp_set(q, 0);
        mp_copy(r, num);
        return true;
    }

#if !defined(__SIZEOF_INT128__)
    return false;
#else
    int m = mp_num_limbs(num);
    int n = mp_num_limbs(den);

    if (n == 1) {
        uint64_t v = den[0];
        mp_uint_t u  = { num[0], num[1], num[2], num[3] };
        mp_uint_t qq = {0};

        uint64_t rem = div_1word(qq, u, 4, v);

        mp_copy(q, qq);
        mp_set(r, rem);
        return true;
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
            carry = x >> (64 - s);
        }

        for (int i = n; i < MP_N64; i++) vnorm[i] = 0;

        carry = 0;

        for (int i = 0; i < MP_N64 + 1; i++) {
            uint64_t x = uraw[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (64 - s);
        }
    }

    int qn = m - n + 1;
    mp_uint_t qlimb = {0};

    const uint64_t v1 = vnorm[n - 1];
    const uint64_t v2 = vnorm[n - 2];

    for (int j = qn - 1; j >= 0; j--) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + n] << 64) | (__uint128_t)unorm[j + n - 1];
        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 64) | (__uint128_t)unorm[j + n - 2];
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
            carry = x << (64 - s);
        }
    }

    mp_copy(q, qlimb);
    mp_copy(r, rlimb);

    return true;
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
            carry = cur >> 64;
        }

        int k = i + MP_N64;
        while (carry) {
            __uint128_t cur = (__uint128_t)out2n[k] + carry;
            out2n[k] = (uint64_t)cur;
            carry = cur >> 64;
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
            carry = x >> (64 - s);
        }

        carry = 0;
        for (int i = 0; i < 2 * MP_N64; i++) {
            uint64_t x = num2n[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (64 - s);
        }
        unorm[2 * MP_N64] = carry;
    }

    const uint64_t v1 = vnorm[MP_N64 - 1];
    const uint64_t v2 = vnorm[MP_N64 - 2];

    for (int j = MP_N64; j >= 0; j--) {
        __uint128_t uj2 = ((__uint128_t)unorm[j + MP_N64] << 64)
                        | (__uint128_t)unorm[j + MP_N64 - 1];

        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for (;;) {
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 64) | (__uint128_t)unorm[j + MP_N64 - 2];
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
            carry = x << (64 - s);
        }
    }

    if (mp_cmp(r, den) >= 0) {
        mp_sub(r, r, den);
    }
}

static inline void mp_mulmod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod)
{
    uint64_t prod[2 * MP_N64];
    mp_mul_full(prod, a, b);
    mp_mod_2n_n(r, prod, mod);
}

void mp_pow_mod(uint64_t *r, const uint64_t *base, const uint64_t *exp, const uint64_t *mod)
{
    mp_uint_t one;
    mp_set(one, 1u);

    if (mp_cmp(mod, one) == 0) { mp_set(r, 0u); return; }

    mp_uint_t bcur;

    if (mp_cmp(base, mod) >= 0) {
        mp_uint_t q, rem;
        if (!mp_div(q, rem, base, mod)) { mp_set(r, 0u); return; }
        mp_copy(bcur, rem);
    } else {
        mp_copy(bcur, base);
    }

    int topBit = -1;
    for (int limb = MP_N64 - 1; limb >= 0 && topBit < 0; --limb) {
        uint64_t w = exp[limb];
        if (!w) continue;

        for (int bit = 63; bit >= 0; --bit) {
            if ((w >> bit) & 1u) {
                topBit = limb * 64 + bit;
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
        const int bit  = i & 63;
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

            if (mp_div(q, rem, r, mod)) mp_copy(r, rem);
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
        if (mp_div(q, rem, av, mod)) mp_copy(av, rem);
    }

    if (mp_is_zero(av)) { mp_zero(r); return; }

    mp_sub(r, mod, av);
}

uint64_t mp_mul(uint64_t *r, const uint64_t *a, uint64_t b) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;

    for (int i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)a[i] * (__uint128_t)b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }

    return (uint64_t)carry;
#else
    uint64_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        uint64_t a0 = (uint32_t)a[i];
        uint64_t a1 = a[i] >> MP_N;
        uint64_t m0 = (uint32_t)b;
        uint64_t m1 = b >> MP_N;

        uint64_t p00 = a0 * m0;
        uint64_t p01 = a0 * m1;
        uint64_t p10 = a1 * m0;
        uint64_t p11 = a1 * m1;

        uint64_t lo = p00 + (p01 << MP_N) + (p10 << MP_N);
        uint64_t hi = p11 + (p01 >> MP_N) + (p10 >> MP_N);

        hi += (lo < p00);

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
        __uint128_t prod = (__uint128_t)a[i] * (__uint128_t)b;
        __uint128_t t = (__uint128_t)r[i] + prod + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }

    return (uint64_t)carry;
#else
    uint64_t tmp[MP_N64];
    uint64_t c = mp_mul(tmp, a, b);
    uint64_t carry = mp_add(r, r, tmp);
    return (c || carry) ? 1u : 0u;
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
    return (int)((a[0] & 1ULL) == 0ULL);
}

static inline void mp_sub_mod(uint64_t *x, const uint64_t *y, const uint64_t *mod) {
    if (mp_cmp(x, y) >= 0) {
        mp_sub(x, x, y);
    } else {
        uint64_t t[MP_N64];
        mp_add(t, x, mod);
        mp_sub(x, t, y);
    }
}

static inline void mp_div2_mod(uint64_t *x, const uint64_t *mod) {
    if (mp_is_even(x)) {
        mp_shr(x, x, 1u);
    } else {
        uint64_t t[MP_N64];
        mp_add(t, x, mod);
        mp_copy(x, t);
        mp_shr(x, x, 1u);
    }
}

bool mp_inv_mod(uint64_t *r, const uint64_t *a, const uint64_t *mod) {
    if (!r || !a || !mod) return false;

    mp_set(r, 0);

    if (mp_is_zero(mod)) return false;

    uint64_t u[MP_N64], v[MP_N64];
    mp_copy(u, a);
    mp_copy(v, mod);

    if (mp_is_zero(u)) {
        mp_set(r, 0);
        return false;
    }

    if (mp_is_even(v)) return false;

    uint64_t x1[MP_N64], x2[MP_N64];
    mp_set(x1, 1u);
    mp_set(x2, 0u);

    while (!(u[0] == 1u && u[1] == 0 && u[2] == 0 && u[3] == 0) &&
           !(v[0] == 1u && v[1] == 0 && v[2] == 0 && v[3] == 0)) {

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
            mp_sub_mod(x1, x2, mod);
        } else {
            mp_sub(v, v, u);
            mp_sub_mod(x2, x1, mod);
        }
    }

    if (u[0] == 1u && u[1] == 0 && u[2] == 0 && u[3] == 0) {
        mp_copy(r, x1);
        return true;
    } else {
        mp_copy(r, x2);
        return true;
    }
}
