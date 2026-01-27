#include "u256.hpp"
#include <cstring>
//#include <cstdlib>
#include <string>
#include <climits>

void mp_set_ui(U256 *r, uint64_t x) {
    r->limb[0] = x;
    r->limb[1] = 0;
    r->limb[2] = 0;
    r->limb[3] = 0;
}

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
    __int128_t t = ( (__int128_t)a ) - b - c;
    *out = (uint64_t)t;
    return (t < 0) ? 1 : 0;
#else
    uint64_t t = a - b;
    uint64_t borrow1 = (a < b);
    uint64_t u = t - c;
    uint64_t borrow2 = (t < c);
    *out = u;
    return borrow1 | borrow2;
#endif
}

void mpn_zero(mp_limb_t *rp, mp_size_t n) {
    std::memset(rp, 0, n * sizeof(mp_limb_t));
}

void mpn_copyi(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n) {
    if (rp == ap) return;
    std::memcpy(rp, ap, n * sizeof(mp_limb_t));
}

int mpn_zero_p(const mp_limb_t *ap, mp_size_t n) {
    mp_limb_t acc = 0;
    for (mp_size_t i = 0; i < n; i++) acc |= ap[i];
    return acc == 0;
}

int mpn_cmp(const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n) {
    for (mp_size_t i = n; i-- > 0;) {
        if (ap[i] < bp[i]) return -1;
        if (ap[i] > bp[i]) return  1;
    }
    return 0;
}

int mpn_cmp_1(const mp_limb_t *ap, mp_size_t n, mp_limb_t b) {
    for (mp_size_t i = n; i-- > 1;) {
        if (ap[i] != 0) return 1;
    }
    if (n == 0) return (b == 0) ? 0 : -1;
    if (ap[0] < b) return -1;
    if (ap[0] > b) return 1;
    return 0;
}

mp_limb_t mpn_add_n(mp_limb_t *rp, const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n) {
    mp_limb_t c = 0;
    for (mp_size_t i = 0; i < n; i++) {
        c = add_carry_u64(ap[i], bp[i], c, &rp[i]);
    }
    return c;
}

mp_limb_t mpn_sub_n(mp_limb_t *rp, const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n) {
    mp_limb_t b = 0;
    for (mp_size_t i = 0; i < n; i++) {
        b = sub_borrow_u64(ap[i], bp[i], b, &rp[i]);
    }
    return b;
}

mp_limb_t mpn_add_1(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, mp_limb_t b) {
    mp_limb_t c = 0;
    if (n == 0) return b ? 1 : 0;

    c = add_carry_u64(ap[0], b, 0, &rp[0]);
    for (mp_size_t i = 1; i < n; i++) {
        c = add_carry_u64(ap[i], 0, c, &rp[i]);
    }
    return c;
}

mp_limb_t mpn_sub_1(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, mp_limb_t b) {
    mp_limb_t br = 0;
    if (n == 0) return b ? 1 : 0;

    br = sub_borrow_u64(ap[0], b, 0, &rp[0]);
    for (mp_size_t i = 1; i < n; i++) {
        br = sub_borrow_u64(ap[i], 0, br, &rp[i]);
    }
    return br;
}

int mpn_tstbit(const mp_limb_t *ap, mp_size_t n, unsigned bit) {
    const unsigned maxBits = (unsigned)(n * 64);
    if (bit >= maxBits) return 0;
    const unsigned w = bit >> 6;
    const unsigned s = bit & 63;
    return (int)((ap[w] >> s) & 1ULL);
}

void mpn_rshift(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, unsigned k) {
    const unsigned maxBits = (unsigned)(n * 64);
    if (k >= maxBits) { mpn_zero(rp, n); return; }
    if (k == 0) { mpn_copyi(rp, ap, n); return; }

    const unsigned wordShift = k >> 6;
    const unsigned bitShift  = k & 63;

    if (bitShift == 0) {
        for (mp_size_t i = 0; i < n; i++) {
            const mp_size_t si = i + wordShift;
            rp[i] = (si < n) ? ap[si] : 0;
        }
        return;
    }

    const unsigned inv = 64 - bitShift;
    for (mp_size_t i = 0; i < n; i++) {
        const mp_size_t si0 = i + wordShift;
        const mp_size_t si1 = si0 + 1;

        const mp_limb_t lo = (si0 < n) ? ap[si0] : 0;
        const mp_limb_t hi = (si1 < n) ? ap[si1] : 0;

        rp[i] = (lo >> bitShift) | (hi << inv);
    }
}

void mpn_lshift(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, unsigned k) {
    const unsigned maxBits = (unsigned)(n * 64);
    if (k >= maxBits) { mpn_zero(rp, n); return; }
    if (k == 0) { mpn_copyi(rp, ap, n); return; }

    const unsigned wordShift = k >> 6;
    const unsigned bitShift  = k & 63;

    if (bitShift == 0) {
        for (mp_size_t i = n; i-- > 0;) {
            const long si = (long)i - (long)wordShift;
            rp[i] = (si >= 0) ? ap[(mp_size_t)si] : 0;
        }
        return;
    }

    const unsigned inv = 64 - bitShift;
    for (mp_size_t i = n; i-- > 0;) {
        const long si0 = (long)i - (long)wordShift;
        const long si1 = si0 - 1;

        const mp_limb_t lo = (si0 >= 0) ? ap[(mp_size_t)si0] : 0;
        const mp_limb_t hi = (si1 >= 0) ? ap[(mp_size_t)si1] : 0;

        rp[i] = (lo << bitShift) | (hi >> inv);
    }
}

void mp_copy(U256 *r, const U256 *a) { mpn_copyi(r->limb, a->limb, 4); }

int mp_is_zero(const U256 *a) { return mpn_zero_p(a->limb, 4); }

int mp_cmp(const U256 *a, const U256 *b) { return mpn_cmp(a->limb, b->limb, 4); }

int mp_cmp_ui(const U256 *a, uint64_t x) { return mpn_cmp_1(a->limb, 4, x); }

uint64_t mp_add(U256 *r, const U256 *a, const U256 *b) { return mpn_add_n(r->limb, a->limb, b->limb, 4); }

uint64_t mp_sub(U256 *r, const U256 *a, const U256 *b) { return mpn_sub_n(r->limb, a->limb, b->limb, 4); }

uint64_t mp_add_ui(U256 *r, const U256 *a, uint64_t b) { return mpn_add_1(r->limb, a->limb, 4, b); }

uint64_t mp_sub_ui(U256 *r, const U256 *a, uint64_t b) { return mpn_sub_1(r->limb, a->limb, 4, b); }

int mp_tstbit(const U256 *a, uint32_t bit) { return mpn_tstbit(a->limb, 4, bit); }

void mp_fdiv_q_2exp(U256 *r, const U256 *a, uint32_t k) { mpn_rshift(r->limb, a->limb, 4, k); }

void mp_export(uint8_t out[32], const U256 *a) {
    for (int i = 0; i < 4; i++) {
        uint64_t w = a->limb[i];
        for (int j = 0; j < 8; j++) {
            out[i * 8 + j] = (uint8_t)(w & 0xFF);
            w >>= 8;
        }
    }
}

static inline int is_space(char c) {
    return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}

static inline uint64_t u256_mul_small_dec(U256 *x, uint32_t m) {
#if defined(__SIZEOF_INT128__)
    __uint128_t carry = 0;
    for (int i = 0; i < 4; i++) {
        __uint128_t t = (__uint128_t)x->limb[i] * m + carry;
        x->limb[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
#else
    (void)x; (void)m;
    return 1;
#endif
}

static inline uint64_t u256_add_small_dec(U256 *x, uint32_t add) {
#if defined(__SIZEOF_INT128__)
    __uint128_t t = (__uint128_t)x->limb[0] + (uint64_t)add;
    x->limb[0] = (uint64_t)t;
    uint64_t c = (uint64_t)(t >> 64);

    for (int i = 1; i < 4 && c; i++) {
        __uint128_t ti = (__uint128_t)x->limb[i] + c;
        x->limb[i] = (uint64_t)ti;
        c = (uint64_t)(ti >> 64);
    }
    return c;
#else
    (void)x; (void)add;
    return 1;
#endif
}

int mp_set_str(U256 *r, const char *str, int base) {
    if (!r || !str) return -1;
    if (base == 0) base = 10;
    if (base != 10) return -1;

    r->limb[0] = r->limb[1] = r->limb[2] = r->limb[3] = 0;

    while (*str && is_space(*str)) str++;
    if (*str == '+') str++;
    if (*str == '-') return -1;

    bool any = false;
    for (; *str; str++) {
        if (*str < '0' || *str > '9') break;
        any = true;

        uint32_t digit = (uint32_t)(*str - '0');

        if (u256_mul_small_dec(r, 10) != 0) return -1;
        if (u256_add_small_dec(r, digit) != 0) return -1;
    }

    if (!any) return -1;

    while (*str && is_space(*str)) str++;
    if (*str != '\0') return -1;

    return 0;
}

int mp_fits_sint(const U256 *a) {
    if (!a) return 0;
    if (a->limb[1] || a->limb[2] || a->limb[3]) return 0;
    return a->limb[0] <= (uint64_t)INT_MAX;
}

void mp_add_mod(U256 *r, const U256 *a, const U256 *b, const U256 *mod) {
    uint64_t carry = mp_add(r, a, b);
    if (carry || mp_cmp(r, mod) >= 0) {
        (void)mp_sub(r, r, mod);
    }
}

void mp_add_ui_mod(U256 *r, const U256 *a, uint64_t b, const U256 *mod) {
    uint64_t carry = mp_add_ui(r, a, b);
    if (carry || mp_cmp(r, mod) >= 0) {
        (void)mp_sub(r, r, mod);
    }
}

void mp_mul_small_mod(U256 *r, const U256 *a, uint32_t m, const U256 *mod) {
    U256 res; mp_set_ui(&res, 0);
    U256 cur; mp_copy(&cur, a);

    uint32_t k = m;
    while (k) {
        if (k & 1u) {
            U256 tmp;
            mp_add_mod(&tmp, &res, &cur, mod);
            res = tmp;
        }
        k >>= 1u;
        if (k) {
            U256 tmp;
            mp_add_mod(&tmp, &cur, &cur, mod);
            cur = tmp;
        }
    }
    *r = res;
}

static inline int digit_val(char c) {
    if (c >= '0' && c <= '9') return (int)(c - '0');
    if (c >= 'a' && c <= 'f') return 10 + (int)(c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (int)(c - 'A');
    return -1;
}

int mp_set_str_mod(U256 *r, const char *str, int base, const U256 *mod) {
    if (!r || !str || !mod) return -1;
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return -1;

    while (*str && is_space(*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    U256 acc; mp_set_ui(&acc, 0);
    bool any = false;

    for (; *str; str++) {
        if (is_space(*str)) break;

        int d = digit_val(*str);
        if (d < 0 || d >= base) return -1;

        U256 t1;
        mp_mul_small_mod(&t1, &acc, (uint32_t)base, mod);

        U256 t2;
        mp_add_ui_mod(&t2, &t1, (uint64_t)d, mod);

        acc = t2;
        any = true;
    }

    if (!any) return -1;

    while (*str && is_space(*str)) str++;
    if (*str != '\0') return -1;

    if (neg && !mp_is_zero(&acc)) {
        U256 tmp;
        (void)mp_sub(&tmp, mod, &acc);
        acc = tmp;
    }

    *r = acc;
    return 0;
}

#if defined(__SIZEOF_INT128__)
static inline uint32_t u256_div_ui(U256 *q, const U256 *a, uint32_t base) {
    uint64_t rem = 0;
    for (int i = 3; i >= 0; i--) {
        __uint128_t cur = (((__uint128_t)rem) << 64) | (__uint128_t)a->limb[i];
        q->limb[i] = (uint64_t)(cur / base);
        rem        = (uint64_t)(cur % base);
    }
    return (uint32_t)rem;
}
#endif

std::string mp_get_str(const U256 *a, int base) {
    if (!a) return std::string();
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return std::string();

    if (mp_is_zero(a)) {
        return std::string("0");
    }

#if !defined(__SIZEOF_INT128__)
    return std::string();
#else
    U256 v; mp_copy(&v, a);
    std::string out;

    while (!mp_is_zero(&v)) {
        U256 q;
        uint32_t rem = u256_div_ui(&q, &v, (uint32_t)base);

        char digit = (rem < 10) ? (char)('0' + rem) : (char)('a' + (rem - 10));
        out.push_back(digit);
        v = q;
    }

    for (size_t i = 0, j = out.size() - 1; i < j; i++, j--) {
        char t = out[i]; out[i] = out[j]; out[j] = t;
    }
    return out;
#endif
}

void mp_export_be(uint8_t out[32], const U256 *a) {
    for (int i = 0; i < 32; i++) {
        uint8_t byte = (uint8_t)((a->limb[i / 8] >> (8 * (i % 8))) & 0xFF);
        out[31 - i] = byte;
    }
}

void mp_import_be(U256 *r, const uint8_t in[32]) {
    r->limb[0] = r->limb[1] = r->limb[2] = r->limb[3] = 0;
    for (int i = 0; i < 32; i++) {
        uint8_t byte = in[31 - i];
        r->limb[i / 8] |= (uint64_t)byte << (8 * (i % 8));
    }
}

static inline int u256_clz64(uint64_t x){ return x?__builtin_clzll(x):64; }
static inline int u256_num_limbs(const U256 *a){ if(a->limb[3]) return 4; if(a->limb[2]) return 3; if(a->limb[1]) return 2; if(a->limb[0]) return 1; return 0; }

static inline uint64_t div_1word(uint64_t *q, const uint64_t *u, int m, uint64_t v){
    __uint128_t rem=0;
    for(int i=m-1;i>=0;i--){
        __uint128_t cur = (rem<<64) | (__uint128_t)u[i];
        q[i] = (uint64_t)(cur / v);
        rem  = (uint64_t)(cur % v);
    }
    return (uint64_t)rem;
}

static inline uint64_t mul_sub_knuth(uint64_t *u, const uint64_t *v, int n, uint64_t qhat){
    __uint128_t carry=0;
    uint64_t borrow=0;
    for(int i=0;i<n;i++){
        __uint128_t prod = (__uint128_t)qhat * (__uint128_t)v[i] + carry;
        uint64_t pl = (uint64_t)prod;
        carry = (prod >> 64);
        __int128_t t = (__int128_t)u[i] - (__int128_t)pl - (__int128_t)borrow;
        u[i] = (uint64_t)t;
        borrow = (t < 0) ? 1 : 0;
    }
    __int128_t ttop = (__int128_t)u[n] - (__int128_t)carry - (__int128_t)borrow;
    u[n] = (uint64_t)ttop;
    return (ttop < 0) ? 1 : 0;
}

static inline uint64_t add_back_knuth(uint64_t *u, const uint64_t *v, int n){
    __uint128_t carry=0;
    for(int i=0;i<n;i++){
        __uint128_t t = (__uint128_t)u[i] + (__uint128_t)v[i] + carry;
        u[i] = (uint64_t)t;
        carry = (t >> 64);
    }
    __uint128_t ttop = (__uint128_t)u[n] + carry;
    u[n] = (uint64_t)ttop;
    return (uint64_t)(ttop >> 64);
}

int mp_divmod(U256 *q, U256 *r, const U256 *num, const U256 *den){
    if(!q||!r||!num||!den) return -1;
    if(mp_is_zero(den)) return -1;

    if(mp_cmp(num,den)<0){ mp_set_ui(q,0); *r=*num; return 0; }

    int m = u256_num_limbs(num);
    int n = u256_num_limbs(den);

    if(n==1){
        uint64_t v = den->limb[0];
        uint64_t u[4] = {num->limb[0],num->limb[1],num->limb[2],num->limb[3]};
        uint64_t qq[4] = {0,0,0,0};
        uint64_t rem = div_1word(qq,u,4,v);
        q->limb[0]=qq[0]; q->limb[1]=qq[1]; q->limb[2]=qq[2]; q->limb[3]=qq[3];
        r->limb[0]=rem; r->limb[1]=r->limb[2]=r->limb[3]=0;
        return 0;
    }

    uint64_t vnorm[4] = {0,0,0,0};
    uint64_t unorm[5] = {0,0,0,0,0};
    uint64_t vraw[4]  = {den->limb[0],den->limb[1],den->limb[2],den->limb[3]};
    uint64_t uraw[5]  = {num->limb[0],num->limb[1],num->limb[2],num->limb[3],0};

    unsigned s = (unsigned)u256_clz64(vraw[n-1]);

    if(s==0){
        for(int i=0;i<n;i++) vnorm[i]=vraw[i];
        for(int i=n;i<4;i++) vnorm[i]=0;
        for(int i=0;i<5;i++) unorm[i]=uraw[i];
    }else{
        uint64_t carry=0;
        for(int i=0;i<n;i++){
            uint64_t x=vraw[i];
            vnorm[i]=(x<<s)|carry;
            carry = x >> (64 - s);
        }
        for(int i=n;i<4;i++) vnorm[i]=0;

        carry=0;
        for(int i=0;i<5;i++){
            uint64_t x=uraw[i];
            unorm[i]=(x<<s)|carry;
            carry = x >> (64 - s);
        }
    }

    int qn = m - n + 1;
    uint64_t qlimb[4] = {0,0,0,0};

    const uint64_t v1 = vnorm[n-1];
    const uint64_t v2 = vnorm[n-2];

    for(int j=qn-1;j>=0;j--){
        __uint128_t uj2 = ((__uint128_t)unorm[j+n] << 64) | (__uint128_t)unorm[j+n-1];
        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        for(;;){
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 64) | (__uint128_t)unorm[j+n-2];
            if(left <= right) break;
            qhat--;
            rhat += v1;
            if(rhat < v1) break;
        }

        uint64_t *u_seg = &unorm[j];
        uint64_t borrow_out = mul_sub_knuth(u_seg, vnorm, n, qhat);
        if(borrow_out){
            (void)add_back_knuth(u_seg, vnorm, n);
            qhat--;
        }
        qlimb[j]=qhat;
    }

    uint64_t rlimb[4] = {0,0,0,0};
    if(s==0){
        for(int i=0;i<n;i++) rlimb[i]=unorm[i];
    }else{
        uint64_t carry=0;
        for(int i=n-1;i>=0;i--){
            uint64_t x=unorm[i];
            rlimb[i]=(x>>s)|carry;
            carry = x << (64 - s);
        }
    }

    q->limb[0]=qlimb[0]; q->limb[1]=qlimb[1]; q->limb[2]=qlimb[2]; q->limb[3]=qlimb[3];
    r->limb[0]=rlimb[0]; r->limb[1]=rlimb[1]; r->limb[2]=rlimb[2]; r->limb[3]=rlimb[3];
    return 0;
}

void mp_set_sint_mod(U256 *r, int64_t x, const U256 *mod) {
    if (x >= 0) {
        mp_set_ui(r, (uint64_t)x);
        if (mp_cmp(r, mod) >= 0) {
            U256 q, rem;
            if (mp_divmod(&q, &rem, r, mod) == 0) *r = rem;
        }
        return;
    }

    uint64_t absv = (uint64_t)(0 - (uint64_t)x);
    U256 a; mp_set_ui(&a, absv);

    if (mp_is_zero(&a)) {
        mp_set_ui(r, 0);
        return;
    }

    if (mp_cmp(&a, mod) >= 0) {
        U256 q, rem;
        if (mp_divmod(&q, &rem, &a, mod) == 0) a = rem;
    }

    if (mp_is_zero(&a)) {
        mp_set_ui(r, 0);
        return;
    }

    (void)mp_sub(r, mod, &a);
}

void mp_shl_2exp(U256 *r, const U256 *a, uint32_t k) { mpn_lshift(r->limb, a->limb, 4, k); }


