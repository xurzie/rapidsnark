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


void mp_set(uint64_t *r, uint64_t x) {
    r[0] = (uint64_t)x;
    for (int i = 1; i < MP_N64; i++) r[i] = 0;
}

void mp_copy(uint64_t *r, const uint64_t *a) {
    if (r == a) return;
    std::memcpy(r, a, MP_N64 * sizeof(uint64_t));
}

int mp_cmp(const uint64_t *a, const uint64_t *b) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

int mp_is_zero(const uint64_t *a) {
    return (a[0] | a[1] | a[2] | a[3]) == 0;
}

void mp_shr(uint64_t *r, const uint64_t *a, uint32_t k) {
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
    for (uint32_t i = 0; i < MP_N64; i++) {
        uint32_t si0 = i + wordShift;
        uint32_t si1 = si0 + 1;

        uint64_t lo = (si0 < 4u) ? a[si0] : 0;
        uint64_t hi = (si1 < 4u) ? a[si1] : 0;

        r[i] = (lo >> bitShift) | (hi << inv);
    }
}

void mp_shl(uint64_t *r, const uint64_t *a, uint32_t k) {
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

uint64_t mp_add(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    uint64_t c = 0;
    c = add_carry(&r[0], a[0], b[0], c);
    c = add_carry(&r[1], a[1], b[1], c);
    c = add_carry(&r[2], a[2], b[2], c);
    c = add_carry(&r[3], a[3], b[3], c);
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

uint64_t mp_add_1(uint64_t *r, const uint64_t *a, uint64_t b) {
    uint64_t c = 0;
    c = add_carry(&r[0], a[0], b, 0);
    c = add_carry(&r[1], a[1], 0, c);
    c = add_carry(&r[2], a[2], 0, c);
    c = add_carry(&r[3], a[3], 0, c);
    return c;
}

uint64_t mp_sub_1(uint64_t *r, const uint64_t *a, uint64_t b) {
    uint64_t br = 0;
    br = sub_borrow(&r[0], a[0], b,  0);
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

void mp_export(uint8_t *out, const uint64_t *a) {
    for (int i = 0; i < MP_N64; i++) {
        uint64_t w = (uint64_t)a[i];
        for (int j = 0; j < 8; j++) {
            out[i * 8 + j] = (uint8_t)(w & 0xFFu);
            w >>= 8;
        }
    }
}

void mp_export_be(uint8_t *out, const uint64_t *a) {
    for (int i = 0; i < MP_N; i++) {
        uint8_t byte = (uint8_t)((a[i / 8] >> (8 * (i % 8))) & 0xFFu);
        out[MP_N - 1 - i] = byte;
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

int mp_set(uint64_t *r, const char *str, int base) {
    if (!r || !str) return -1;
    if (base == 0) base = 10;
    if (base != 10) return -1;

    mp_zero(r);

    while (*str && std::isspace((unsigned char)*str)) str++;
    if (*str == '+') str++;
    if (*str == '-') return -1;

    bool any = false;
    for (; *str; str++) {
        if (*str < '0' || *str > '9') break;
        any = true;

#if !defined(__SIZEOF_INT128__)
        return -1;
#else
        uint32_t digit = (uint32_t)(*str - '0');
        if (mp_mul_small_dec(r, 10) != 0) return -1;
        if (mp_add_small_dec(r, digit) != 0) return -1;
#endif
    }

    if (!any) return -1;

    while (*str && std::isspace((unsigned char)*str)) str++;
    if (*str != '\0') return -1;

    return 0;
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

std::string mp_get_str(const uint64_t *a, int base) {
    if (!a) return std::string();
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return std::string();

    if (mp_is_zero(a)) return std::string("0");

#if !defined(__SIZEOF_INT128__)
    return std::string();
#else
    uint64_t v[4]; mp_copy(v, a);
    std::string out;

    while (!mp_is_zero(v)) {
        uint64_t q[4];
        uint32_t rem = mp_div(q, v, (uint32_t)base);
        char digit = (rem < 10) ? (char)('0' + rem) : (char)('a' + (rem - 10));
        out.push_back(digit);
        mp_copy(v, q);
    }

    for (size_t i = 0, j = out.size() - 1; i < j; i++, j--) {
        char t = out[i]; out[i] = out[j]; out[j] = t;
    }
    return out;
#endif
}

int mp_fits_int32(const uint64_t *a) {
    if (!a) return 0;
    if (a[1] || a[2] || a[3]) return 0;
    return a[0] <= (uint64_t)INT_MAX;
}

void mp_add_mod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod) {
    uint64_t carry = mp_add(r, a, b);
    if (carry || mp_cmp(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
}

void mp_mul_mod(uint64_t *r, const uint64_t *a, uint32_t m, const uint64_t *mod) {
    uint64_t res[4]; mp_set(res, 0);
    uint64_t cur[4]; mp_copy(cur, a);

    uint32_t k = m;
    while (k) {
        if (k & 1u) {
            uint64_t tmp[4];
            mp_add_mod(tmp, res, cur, mod);
            mp_copy(res, tmp);
        }
        k >>= 1u;
        if (k) {
            uint64_t tmp[4];
            mp_add_mod(tmp, cur, cur, mod);
            mp_copy(cur, tmp);
        }
    }

    mp_copy(r, res);
}

int mp_set_mod(uint64_t *r, const char *str, int base, const uint64_t *mod) {
    if (!r || !str || !mod) return -1;
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return -1;

    while (*str && std::isspace((unsigned char)*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    uint64_t acc[4]; mp_set(acc, 0);
    bool any = false;

    for (; *str; str++) {
        if (std::isspace((unsigned char)*str)) break;

        unsigned char uc = (unsigned char)*str;

        if (!std::isxdigit(uc)) return -1;

        int d;
        if (std::isdigit(uc)) {
            d = (int)(uc - (unsigned char)'0');
        } else {
            d = (int)(std::tolower(uc) - (unsigned char)'a') + 10;
        }

        if (d >= base) return -1;

        uint64_t t1[4];
        mp_mul_mod(t1, acc, (uint32_t)base, mod);

        uint64_t dv[4]; mp_set(dv, (uint64_t)d);

        uint64_t t2[4];
        mp_add_mod(t2, t1, dv, mod);

        mp_copy(acc, t2);
        any = true;
    }

    if (!any) return -1;

    while (*str && std::isspace((unsigned char)*str)) str++;
    if (*str != '\0') return -1;

    if (neg && !mp_is_zero(acc)) {
        uint64_t tmp[4];
        mp_sub(tmp, mod, acc);
        mp_copy(acc, tmp);
    }

    mp_copy(r, acc);
    return 0;
}

static inline int mp_clz64(uint64_t x) { return x ? __builtin_clzll(x) : 64; }

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

int mp_divmod(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den) {
    if (!q || !r || !num || !den) return -1;
    if (mp_is_zero(den)) return -1;

    if (mp_cmp(num, den) < 0) {
        mp_set(q, 0);
        mp_copy(r, num);
        return 0;
    }

#if !defined(__SIZEOF_INT128__)
    return -1;
#else
    int m = mp_num_limbs(num);
    int n = mp_num_limbs(den);

    if (n == 1) {
        uint64_t v = (uint64_t)den[0];
        uint64_t u[4]  = {(uint64_t)num[0], (uint64_t)num[1], (uint64_t)num[2], (uint64_t)num[3]};
        uint64_t qq[4] = {0, 0, 0, 0};
        uint64_t rem = div_1word(qq, u, 4, v);

        q[0] = (uint64_t)qq[0];
        q[1] = (uint64_t)qq[1];
        q[2] = (uint64_t)qq[2];
        q[3] = (uint64_t)qq[3];

        mp_set(r, rem);
        return 0;
    }

    uint64_t vnorm[4] = {0,0,0,0};
    uint64_t unorm[5] = {0,0,0,0,0};
    uint64_t vraw[4]  = {(uint64_t)den[0], (uint64_t)den[1], (uint64_t)den[2], (uint64_t)den[3]};
    uint64_t uraw[5]  = {(uint64_t)num[0], (uint64_t)num[1], (uint64_t)num[2], (uint64_t)num[3], 0};

    unsigned s = (unsigned)mp_clz64(vraw[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; i++) vnorm[i] = vraw[i];
        for (int i = n; i < MP_N64; i++) vnorm[i] = 0;
        for (int i = 0; i < 5; i++) unorm[i] = uraw[i];
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; i++) {
            uint64_t x = vraw[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (64 - s);
        }
        for (int i = n; i < MP_N64; i++) vnorm[i] = 0;

        carry = 0;
        for (int i = 0; i < 5; i++) {
            uint64_t x = uraw[i];
            unorm[i] = (x << s) | carry;
            carry = x >> (64 - s);
        }
    }

    int qn = m - n + 1;
    uint64_t qlimb[4] = {0,0,0,0};

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

    uint64_t rlimb[4] = {0,0,0,0};
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

    q[0] = (uint64_t)qlimb[0];
    q[1] = (uint64_t)qlimb[1];
    q[2] = (uint64_t)qlimb[2];
    q[3] = (uint64_t)qlimb[3];

    r[0] = (uint64_t)rlimb[0];
    r[1] = (uint64_t)rlimb[1];
    r[2] = (uint64_t)rlimb[2];
    r[3] = (uint64_t)rlimb[3];

    return 0;
#endif
}

void mp_set_mod(uint64_t *r, int64_t x, const uint64_t *mod) {
    if (x >= 0) {
        mp_set(r, (uint64_t)x);
        if (mp_cmp(r, mod) >= 0) {
            uint64_t q[4], rem[4];
            if (mp_divmod(q, rem, r, mod) == 0) mp_copy(r, rem);
        }
        return;
    }

#if defined(__SIZEOF_INT128__)
    uint64_t absv = (uint64_t)(-((__int128_t)x));
#else
    uint64_t absv = (uint64_t)(-(uint64_t)x);
#endif

    uint64_t a[4]; mp_set(a, absv);

    if (mp_is_zero(a)) { mp_zero(r); return; }

    if (mp_cmp(a, mod) >= 0) {
        uint64_t q[4], rem[4];
        if (mp_divmod(q, rem, a, mod) == 0) mp_copy(a, rem);
    }

    if (mp_is_zero(a)) { mp_zero(r); return; }

    mp_sub(r, mod, a);
}


uint64_t mp_mul_1(uint64_t *r, const uint64_t *a, uint64_t m) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        __uint128_t t = (__uint128_t)a[i] * (__uint128_t)m + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
#else
    uint64_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        uint64_t a0 = (uint32_t)a[i];
        uint64_t a1 = a[i] >> MP_N;
        uint64_t m0 = (uint32_t)m;
        uint64_t m1 = m >> MP_N;

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

uint64_t mp_addmul_1(uint64_t *r, const uint64_t *a, uint64_t m) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (int i = 0; i < MP_N64; i++) {
        __uint128_t prod = (__uint128_t)a[i] * (__uint128_t)m;
        __uint128_t t = (__uint128_t)r[i] + prod + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
#else
    uint64_t tmp[MP_N64];
    uint64_t c = mp_mul_1(tmp, a, m);
    uint64_t carry = mp_add(r, r, tmp);
    return (c || carry) ? 1u : 0u;
#endif
}

uint64_t mp_lshift(uint64_t *r, const uint64_t *a, unsigned int cnt) {
    if (cnt == 0) { if (r != a) mp_copy(r, a); return 0; }
    uint64_t carry = 0;
    unsigned int rcnt = 64u - cnt;
    for (int i = 0; i < MP_N64; i++) {
        uint64_t x = a[i];
        r[i] = (x << cnt) | carry;
        carry = (rcnt == 64u) ? 0 : (x >> rcnt);
    }
    return carry;
}

uint64_t mp_rshift(uint64_t *r, const uint64_t *a, unsigned int cnt) {
    if (cnt == 0) { if (r != a) mp_copy(r, a); return 0; }
    uint64_t carry = 0;
    unsigned int lcnt = 64u - cnt;
    for (int i = MP_N64 - 1; i >= 0; i--) {
        uint64_t x = a[i];
        r[i] = (x >> cnt) | carry;
        carry = (lcnt == 64u) ? 0 : (x << lcnt);
    }
    return carry;
}

void mp_and(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] & b[i];
}
void mp_ior(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] | b[i];
}
void mp_xor(uint64_t *r, const uint64_t *a, const uint64_t *b) {
    for (int i = 0; i < MP_N64; i++) r[i] = a[i] ^ b[i];
}
void mp_com(uint64_t *r, const uint64_t *a) {
    for (int i = 0; i < MP_N64; i++) r[i] = ~a[i];
}



uint64_t mp_add(uint64_t *rp, const uint64_t *ap, size_t an, const uint64_t *bp, size_t bn) {
    size_t n = (an > bn) ? an : bn;
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t a = (i < an) ? ap[i] : 0;
        uint64_t b = (i < bn) ? bp[i] : 0;
        carry = add_carry(&rp[i], a, b, carry);
    }
    return carry;
}


uint64_t mp_add(uint64_t *r, const uint64_t *a, uint64_t b) { return mp_add_1(r, a, b); }

uint64_t mp_sub(uint64_t *r, const uint64_t *a, uint64_t b) { return mp_sub_1(r, a, b); }
