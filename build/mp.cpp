#include "mp.hpp"

#include <cstdint>
#include <string>
#include <climits>

static inline uint64_t add_carry_u64(uint64_t a, uint64_t b, uint64_t c, uint64_t *out) {
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

static inline uint64_t sub_borrow_u64(uint64_t a, uint64_t b, uint64_t c, uint64_t *out) {
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

static inline void mp_zero_4(mp_limb_t *r) {
    r[0] = r[1] = r[2] = r[3] = 0;
}

static inline void mp_copy_4(mp_limb_t *r, const mp_limb_t *a) {
    if (r == a) return;
    r[0] = a[0]; r[1] = a[1]; r[2] = a[2]; r[3] = a[3];
}

static inline int mp_cmp_4(const mp_limb_t *a, const mp_limb_t *b) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

static inline int mp_is_zero_4(const mp_limb_t *a) {
    return (a[0] | a[1] | a[2] | a[3]) == 0;
}

static inline void mp_rshift_4(mp_limb_t *r, const mp_limb_t *a, uint32_t k) {
    if (k >= 256u) { mp_zero_4(r); return; }
    if (k == 0) { mp_copy_4(r, a); return; }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63u;

    if (wordShift >= 4u) { mp_zero_4(r); return; }

    if (bitShift == 0) {
        for (uint32_t i = 0; i < 4; i++) {
            uint32_t si = i + wordShift;
            r[i] = (si < 4u) ? a[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64u - bitShift;
    for (uint32_t i = 0; i < 4; i++) {
        uint32_t si0 = i + wordShift;
        uint32_t si1 = si0 + 1;

        uint64_t lo = (si0 < 4u) ? a[si0] : 0;
        uint64_t hi = (si1 < 4u) ? a[si1] : 0;

        r[i] = (lo >> bitShift) | (hi << inv);
    }
}

static inline void mp_lshift_4(mp_limb_t *r, const mp_limb_t *a, uint32_t k) {
    if (k >= 256u) { mp_zero_4(r); return; }
    if (k == 0) { mp_copy_4(r, a); return; }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63u;

    if (wordShift >= 4u) { mp_zero_4(r); return; }

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

/* ===== mp.hpp API ===== */

void mp_set_ui(mp_limb_t *r, uint64_t x) {
    r[0] = (mp_limb_t)x;
    r[1] = r[2] = r[3] = 0;
}

void mp_copy(mp_limb_t *r, const mp_limb_t *a) {
    mp_copy_4(r, a);
}

int mp_cmp(const mp_limb_t *a, const mp_limb_t *b) {
    return mp_cmp_4(a, b);
}

int mp_is_zero(const mp_limb_t *a) {
    return mp_is_zero_4(a);
}

uint64_t mp_add(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b) {
    uint64_t c = 0;
    c = add_carry_u64(a[0], b[0], c, &r[0]);
    c = add_carry_u64(a[1], b[1], c, &r[1]);
    c = add_carry_u64(a[2], b[2], c, &r[2]);
    c = add_carry_u64(a[3], b[3], c, &r[3]);
    return c;
}

uint64_t mp_sub(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b) {
    uint64_t br = 0;
    br = sub_borrow_u64(a[0], b[0], br, &r[0]);
    br = sub_borrow_u64(a[1], b[1], br, &r[1]);
    br = sub_borrow_u64(a[2], b[2], br, &r[2]);
    br = sub_borrow_u64(a[3], b[3], br, &r[3]);
    return br;
}

uint64_t mp_add(mp_limb_t *r, const mp_limb_t *a, uint64_t b) {
    uint64_t c = 0;
    c = add_carry_u64(a[0], b, 0, &r[0]);
    c = add_carry_u64(a[1], 0, c, &r[1]);
    c = add_carry_u64(a[2], 0, c, &r[2]);
    c = add_carry_u64(a[3], 0, c, &r[3]);
    return c;
}

uint64_t mp_sub(mp_limb_t *r, const mp_limb_t *a, uint64_t b) {
    uint64_t br = 0;
    br = sub_borrow_u64(a[0], b, 0, &r[0]);
    br = sub_borrow_u64(a[1], 0, br, &r[1]);
    br = sub_borrow_u64(a[2], 0, br, &r[2]);
    br = sub_borrow_u64(a[3], 0, br, &r[3]);
    return br;
}

int mp_tstbit(const mp_limb_t *a, size_t bit) {
    if (bit >= 256u) return 0;
    const size_t w = bit >> 6;
    const size_t s = bit & 63u;
    return (int)((a[w] >> s) & 1ULL);
}

void mp_fdiv_q_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k) {
    mp_rshift_4(r, a, k);
}

void mp_shl_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k) {
    mp_lshift_4(r, a, k);
}

void mp_shr_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k) {
    mp_rshift_4(r, a, k);
}

void mp_export(uint8_t *out, const mp_limb_t *a) {
    for (int i = 0; i < 4; i++) {
        uint64_t w = (uint64_t)a[i];
        for (int j = 0; j < 8; j++) {
            out[i * 8 + j] = (uint8_t)(w & 0xFFu);
            w >>= 8;
        }
    }
}

void mp_export_be(uint8_t *out, const mp_limb_t *a) {
    for (int i = 0; i < 32; i++) {
        uint8_t byte = (uint8_t)((a[i / 8] >> (8 * (i % 8))) & 0xFFu);
        out[31 - i] = byte;
    }
}

void mp_import_be(mp_limb_t *r, const uint8_t in[32]) {
    mp_zero_4(r);
    for (int i = 0; i < 32; i++) {
        uint8_t byte = in[31 - i];
        r[i / 8] |= (mp_limb_t)byte << (8 * (i % 8));
    }
}

/* ===== String parsing/formatting ===== */

static inline int is_space(char c) {
    return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}

#if defined(__SIZEOF_INT128__)
static inline uint64_t mp_mul_small_dec_4(mp_limb_t *x, uint32_t m) {
    __uint128_t carry = 0;
    for (int i = 0; i < 4; i++) {
        __uint128_t t = (__uint128_t)x[i] * (uint64_t)m + carry;
        x[i] = (mp_limb_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t mp_add_small_dec_4(mp_limb_t *x, uint32_t add) {
    __uint128_t t = (__uint128_t)x[0] + (uint64_t)add;
    x[0] = (mp_limb_t)t;
    uint64_t c = (uint64_t)(t >> 64);

    for (int i = 1; i < 4 && c; i++) {
        __uint128_t ti = (__uint128_t)x[i] + c;
        x[i] = (mp_limb_t)ti;
        c = (uint64_t)(ti >> 64);
    }
    return c;
}
#endif

int mp_set_str(mp_limb_t *r, const char *str, int base) {
    if (!r || !str) return -1;
    if (base == 0) base = 10;
    if (base != 10) return -1;

    mp_set_ui(r, 0);

    while (*str && is_space(*str)) str++;
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
        if (mp_mul_small_dec_4(r, 10) != 0) return -1;
        if (mp_add_small_dec_4(r, digit) != 0) return -1;
#endif
    }

    if (!any) return -1;

    while (*str && is_space(*str)) str++;
    if (*str != '\0') return -1;

    return 0;
}

#if defined(__SIZEOF_INT128__)
static inline uint32_t mp_div_ui_4(mp_limb_t *q, const mp_limb_t *a, uint32_t base) {
    uint64_t rem = 0;
    for (int i = 3; i >= 0; i--) {
        __uint128_t cur = (((__uint128_t)rem) << 64) | (__uint128_t)a[i];
        q[i] = (mp_limb_t)(cur / base);
        rem  = (uint64_t)(cur % base);
    }
    return (uint32_t)rem;
}
#endif

std::string mp_get_str(const mp_limb_t *a, int base) {
    if (!a) return std::string();
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return std::string();

    if (mp_is_zero(a)) return std::string("0");

#if !defined(__SIZEOF_INT128__)
    return std::string();
#else
    mp_limb_t v[4]; mp_copy_4(v, a);
    std::string out;

    while (!mp_is_zero_4(v)) {
        mp_limb_t q[4];
        uint32_t rem = mp_div_ui_4(q, v, (uint32_t)base);
        char digit = (rem < 10) ? (char)('0' + rem) : (char)('a' + (rem - 10));
        out.push_back(digit);
        mp_copy_4(v, q);
    }

    for (size_t i = 0, j = out.size() - 1; i < j; i++, j--) {
        char t = out[i]; out[i] = out[j]; out[j] = t;
    }
    return out;
#endif
}

/* ===== Mod helpers ===== */

int mp_fits_sint(const mp_limb_t *a) {
    if (!a) return 0;
    if (a[1] || a[2] || a[3]) return 0;
    return a[0] <= (uint64_t)INT_MAX;
}

void mp_add_mod(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b, const mp_limb_t *mod) {
    uint64_t carry = mp_add(r, a, b);
    if (carry || mp_cmp_4(r, mod) >= 0) {
        mp_sub(r, r, mod);
    }
}

void mp_mul_small_mod(mp_limb_t *r, const mp_limb_t *a, uint32_t m, const mp_limb_t *mod) {
    mp_limb_t res[4]; mp_set_ui(res, 0);
    mp_limb_t cur[4]; mp_copy_4(cur, a);

    uint32_t k = m;
    while (k) {
        if (k & 1u) {
            mp_limb_t tmp[4];
            mp_add_mod(tmp, res, cur, mod);
            mp_copy_4(res, tmp);
        }
        k >>= 1u;
        if (k) {
            mp_limb_t tmp[4];
            mp_add_mod(tmp, cur, cur, mod);
            mp_copy_4(cur, tmp);
        }
    }

    mp_copy_4(r, res);
}

static inline int digit_val(char c) {
    if (c >= '0' && c <= '9') return (int)(c - '0');
    if (c >= 'a' && c <= 'f') return 10 + (int)(c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (int)(c - 'A');
    return -1;
}

int mp_set_str_mod(mp_limb_t *r, const char *str, int base, const mp_limb_t *mod) {
    if (!r || !str || !mod) return -1;
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return -1;

    while (*str && is_space(*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    mp_limb_t acc[4]; mp_set_ui(acc, 0);
    bool any = false;

    for (; *str; str++) {
        if (is_space(*str)) break;

        int d = digit_val(*str);
        if (d < 0 || d >= base) return -1;

        mp_limb_t t1[4];
        mp_mul_small_mod(t1, acc, (uint32_t)base, mod);

        mp_limb_t dv[4]; mp_set_ui(dv, (uint64_t)d);

        mp_limb_t t2[4];
        mp_add_mod(t2, t1, dv, mod);

        mp_copy_4(acc, t2);
        any = true;
    }

    if (!any) return -1;

    while (*str && is_space(*str)) str++;
    if (*str != '\0') return -1;

    if (neg && !mp_is_zero_4(acc)) {
        mp_limb_t tmp[4];
        mp_sub(tmp, mod, acc);
        mp_copy_4(acc, tmp);
    }

    mp_copy_4(r, acc);
    return 0;
}

/* ===== Division ===== */

static inline int mp_clz64(uint64_t x) { return x ? __builtin_clzll(x) : 64; }

static inline int mp_num_limbs_4(const mp_limb_t *a) {
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

int mp_divmod(mp_limb_t *q, mp_limb_t *r, const mp_limb_t *num, const mp_limb_t *den) {
    if (!q || !r || !num || !den) return -1;
    if (mp_is_zero_4(den)) return -1;

    if (mp_cmp_4(num, den) < 0) {
        mp_set_ui(q, 0);
        mp_copy_4(r, num);
        return 0;
    }

#if !defined(__SIZEOF_INT128__)
    return -1;
#else
    int m = mp_num_limbs_4(num);
    int n = mp_num_limbs_4(den);

    if (n == 1) {
        uint64_t v = (uint64_t)den[0];
        uint64_t u[4]  = {(uint64_t)num[0], (uint64_t)num[1], (uint64_t)num[2], (uint64_t)num[3]};
        uint64_t qq[4] = {0, 0, 0, 0};
        uint64_t rem = div_1word(qq, u, 4, v);

        q[0] = (mp_limb_t)qq[0];
        q[1] = (mp_limb_t)qq[1];
        q[2] = (mp_limb_t)qq[2];
        q[3] = (mp_limb_t)qq[3];

        mp_set_ui(r, rem);
        return 0;
    }

    uint64_t vnorm[4] = {0,0,0,0};
    uint64_t unorm[5] = {0,0,0,0,0};
    uint64_t vraw[4]  = {(uint64_t)den[0], (uint64_t)den[1], (uint64_t)den[2], (uint64_t)den[3]};
    uint64_t uraw[5]  = {(uint64_t)num[0], (uint64_t)num[1], (uint64_t)num[2], (uint64_t)num[3], 0};

    unsigned s = (unsigned)mp_clz64(vraw[n - 1]);

    if (s == 0) {
        for (int i = 0; i < n; i++) vnorm[i] = vraw[i];
        for (int i = n; i < 4; i++) vnorm[i] = 0;
        for (int i = 0; i < 5; i++) unorm[i] = uraw[i];
    } else {
        uint64_t carry = 0;
        for (int i = 0; i < n; i++) {
            uint64_t x = vraw[i];
            vnorm[i] = (x << s) | carry;
            carry = x >> (64 - s);
        }
        for (int i = n; i < 4; i++) vnorm[i] = 0;

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

    q[0] = (mp_limb_t)qlimb[0];
    q[1] = (mp_limb_t)qlimb[1];
    q[2] = (mp_limb_t)qlimb[2];
    q[3] = (mp_limb_t)qlimb[3];

    r[0] = (mp_limb_t)rlimb[0];
    r[1] = (mp_limb_t)rlimb[1];
    r[2] = (mp_limb_t)rlimb[2];
    r[3] = (mp_limb_t)rlimb[3];

    return 0;
#endif
}

void mp_set_sint_mod(mp_limb_t *r, int64_t x, const mp_limb_t *mod) {
    if (x >= 0) {
        mp_set_ui(r, (uint64_t)x);
        if (mp_cmp_4(r, mod) >= 0) {
            mp_limb_t q[4], rem[4];
            if (mp_divmod(q, rem, r, mod) == 0) mp_copy_4(r, rem);
        }
        return;
    }

#if defined(__SIZEOF_INT128__)
    uint64_t absv = (uint64_t)(-((__int128_t)x));
#else
    uint64_t absv = (uint64_t)(-(uint64_t)x);
#endif

    mp_limb_t a[4]; mp_set_ui(a, absv);

    if (mp_is_zero_4(a)) { mp_set_ui(r, 0); return; }

    if (mp_cmp_4(a, mod) >= 0) {
        mp_limb_t q[4], rem[4];
        if (mp_divmod(q, rem, a, mod) == 0) mp_copy_4(a, rem);
    }

    if (mp_is_zero_4(a)) { mp_set_ui(r, 0); return; }

    mp_sub(r, mod, a);
}
