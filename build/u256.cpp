#include "u256.hpp"
#include <cstring>
#include <cstdlib>
#include <string>
#include <climits>

void u256_set_ui(U256* r, uint64_t x) {
    r->limb[0] = x;
    r->limb[1] = 0;
    r->limb[2] = 0;
    r->limb[3] = 0;
}

void u256_copy(U256* r, const U256* a) {
    std::memcpy(r->limb, a->limb, sizeof(r->limb));
}

int u256_is_zero(const U256* a) {
    return (a->limb[0] | a->limb[1] | a->limb[2] | a->limb[3]) == 0;
}

int u256_cmp(const U256* a, const U256* b) {
    for (int i = 3; i >= 0; --i) {
        if (a->limb[i] < b->limb[i]) return -1;
        if (a->limb[i] > b->limb[i]) return  1;
    }
    return 0;
}

int u256_cmp_ui(const U256* a, uint64_t x) {
    if (a->limb[3] | a->limb[2] | a->limb[1]) return 1;
    if (a->limb[0] < x) return -1;
    if (a->limb[0] > x) return 1;
    return 0;
}

static inline uint64_t add_carry_u64(uint64_t a, uint64_t b, uint64_t c, uint64_t* out) {
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

static inline uint64_t sub_borrow_u64(uint64_t a, uint64_t b, uint64_t c, uint64_t* out) {
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

uint64_t u256_add(U256* r, const U256* a, const U256* b) {
    uint64_t c = 0;
    c = add_carry_u64(a->limb[0], b->limb[0], c, &r->limb[0]);
    c = add_carry_u64(a->limb[1], b->limb[1], c, &r->limb[1]);
    c = add_carry_u64(a->limb[2], b->limb[2], c, &r->limb[2]);
    c = add_carry_u64(a->limb[3], b->limb[3], c, &r->limb[3]);
    return c;
}

uint64_t u256_sub(U256* r, const U256* a, const U256* b) {
    uint64_t c = 0;
    c = sub_borrow_u64(a->limb[0], b->limb[0], c, &r->limb[0]);
    c = sub_borrow_u64(a->limb[1], b->limb[1], c, &r->limb[1]);
    c = sub_borrow_u64(a->limb[2], b->limb[2], c, &r->limb[2]);
    c = sub_borrow_u64(a->limb[3], b->limb[3], c, &r->limb[3]);
    return c;
}

uint64_t u256_add_ui(U256* r, const U256* a, uint64_t b) {
    uint64_t c = 0;
    c = add_carry_u64(a->limb[0], b, c, &r->limb[0]);
    c = add_carry_u64(a->limb[1], 0, c, &r->limb[1]);
    c = add_carry_u64(a->limb[2], 0, c, &r->limb[2]);
    c = add_carry_u64(a->limb[3], 0, c, &r->limb[3]);
    return c;
}

uint64_t u256_sub_ui(U256* r, const U256* a, uint64_t b) {
    uint64_t c = 0;
    c = sub_borrow_u64(a->limb[0], b, c, &r->limb[0]);
    c = sub_borrow_u64(a->limb[1], 0, c, &r->limb[1]);
    c = sub_borrow_u64(a->limb[2], 0, c, &r->limb[2]);
    c = sub_borrow_u64(a->limb[3], 0, c, &r->limb[3]);
    return c;
}

int u256_tstbit(const U256* a, uint32_t bit) {
    if (bit >= 256) return 0;
    uint32_t w = bit >> 6;
    uint32_t s = bit & 63;
    return (int)((a->limb[w] >> s) & 1ULL);
}

void u256_fdiv_q_2exp(U256* r, const U256* a, uint32_t k) {
    if (k >= 256) {
        r->limb[0] = r->limb[1] = r->limb[2] = r->limb[3] = 0;
        return;
    }

    U256 tmp;
    const U256* src = a;
    if (r == a) {
        tmp = *a;
        src = &tmp;
    }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63;

    if (bitShift == 0) {
        for (int i = 0; i < 4; i++) {
            int si = i + (int)wordShift;
            r->limb[i] = (si < 4) ? src->limb[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64 - bitShift;
    for (int i = 0; i < 4; i++) {
        int si0 = i + (int)wordShift;
        int si1 = si0 + 1;

        uint64_t lo = (si0 < 4) ? src->limb[si0] : 0;
        uint64_t hi = (si1 < 4) ? src->limb[si1] : 0;

        r->limb[i] = (lo >> bitShift) | (hi << inv);
    }
}

void u256_export(uint8_t out[32], const U256* a) {
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

static inline uint64_t u256_mul_small_dec(U256* x, uint32_t m) {
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

static inline uint64_t u256_add_small_dec(U256* x, uint32_t add) {
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

int u256_set_str(U256* r, const char* str, int base) {
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

int u256_fits_sint(const U256* a) {
    if (!a) return 0;
    if (a->limb[1] || a->limb[2] || a->limb[3]) return 0;
    return a->limb[0] <= (uint64_t)INT_MAX;
}

void u256_add_mod(U256* r, const U256* a, const U256* b, const U256* mod) {
    uint64_t carry = u256_add(r, a, b);
    if (carry || u256_cmp(r, mod) >= 0) {
        (void)u256_sub(r, r, mod);
    }
}

void u256_add_ui_mod(U256* r, const U256* a, uint64_t b, const U256* mod) {
    uint64_t carry = u256_add_ui(r, a, b);
    if (carry || u256_cmp(r, mod) >= 0) {
        (void)u256_sub(r, r, mod);
    }
}

void u256_mul_small_mod(U256* r, const U256* a, uint32_t m, const U256* mod) {
    U256 res; u256_set_ui(&res, 0);
    U256 cur; u256_copy(&cur, a);

    uint32_t k = m;
    while (k) {
        if (k & 1u) {
            U256 tmp;
            u256_add_mod(&tmp, &res, &cur, mod);
            res = tmp;
        }
        k >>= 1u;
        if (k) {
            U256 tmp;
            u256_add_mod(&tmp, &cur, &cur, mod);
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

int u256_set_str_mod(U256* r, const char* str, int base, const U256* mod) {
    if (!r || !str || !mod) return -1;
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return -1;

    while (*str && is_space(*str)) str++;

    bool neg = false;
    if (*str == '+') str++;
    else if (*str == '-') { neg = true; str++; }

    U256 acc; u256_set_ui(&acc, 0);
    bool any = false;

    for (; *str; str++) {
        if (is_space(*str)) break;

        int d = digit_val(*str);
        if (d < 0 || d >= base) return -1;

        U256 t1;
        u256_mul_small_mod(&t1, &acc, (uint32_t)base, mod);

        U256 t2;
        u256_add_ui_mod(&t2, &t1, (uint64_t)d, mod);

        acc = t2;
        any = true;
    }

    if (!any) return -1;

    while (*str && is_space(*str)) str++;
    if (*str != '\0') return -1;

    if (neg && !u256_is_zero(&acc)) {
        U256 tmp;
        (void)u256_sub(&tmp, mod, &acc);
        acc = tmp;
    }

    *r = acc;
    return 0;
}

#if defined(__SIZEOF_INT128__)
static inline uint32_t u256_div_ui(U256* q, const U256* a, uint32_t base) {
    uint64_t rem = 0;
    for (int i = 3; i >= 0; i--) {
        __uint128_t cur = (((__uint128_t)rem) << 64) | (__uint128_t)a->limb[i];
        q->limb[i] = (uint64_t)(cur / base);
        rem        = (uint64_t)(cur % base);
    }
    return (uint32_t)rem;
}
#endif

char* u256_get_str_alloc(const U256* a, int base) {
    if (!a) return nullptr;
    if (base == 0) base = 10;
    if (base < 2 || base > 16) return nullptr;

    if (u256_is_zero(a)) {
        char* p = (char*)std::malloc(2);
        if (!p) return nullptr;
        p[0] = '0'; p[1] = 0;
        return p;
    }

#if !defined(__SIZEOF_INT128__)
    return nullptr;
#else
    U256 v; u256_copy(&v, a);
    std::string out;

    while (!u256_is_zero(&v)) {
        U256 q;
        uint32_t rem = u256_div_ui(&q, &v, (uint32_t)base);

        char digit = (rem < 10) ? (char)('0' + rem) : (char)('a' + (rem - 10));
        out.push_back(digit);
        v = q;
    }

    for (size_t i = 0, j = out.size() - 1; i < j; i++, j--) {
        char t = out[i]; out[i] = out[j]; out[j] = t;
    }

    char* p = (char*)std::malloc(out.size() + 1);
    if (!p) return nullptr;
    std::memcpy(p, out.c_str(), out.size() + 1);
    return p;
#endif
}

void u256_export_be(uint8_t out[32], const U256* a) {
    for (int i = 0; i < 32; i++) {
        uint8_t byte = (uint8_t)((a->limb[i / 8] >> (8 * (i % 8))) & 0xFF);
        out[31 - i] = byte;
    }
}

void u256_import_be(U256* r, const uint8_t in[32]) {
    r->limb[0] = r->limb[1] = r->limb[2] = r->limb[3] = 0;
    for (int i = 0; i < 32; i++) {
        uint8_t byte = in[31 - i];
        r->limb[i / 8] |= (uint64_t)byte << (8 * (i % 8));
    }
}

static inline int u256_clz64(uint64_t x){ return x?__builtin_clzll(x):64; }
static inline int u256_num_limbs(const U256* a){ if(a->limb[3]) return 4; if(a->limb[2]) return 3; if(a->limb[1]) return 2; if(a->limb[0]) return 1; return 0; }

static inline uint64_t div_1word(uint64_t* q, const uint64_t* u, int m, uint64_t v){
    __uint128_t rem=0;
    for(int i=m-1;i>=0;i--){
        __uint128_t cur = (rem<<64) | (__uint128_t)u[i];
        q[i] = (uint64_t)(cur / v);
        rem  = (uint64_t)(cur % v);
    }
    return (uint64_t)rem;
}

static inline uint64_t mul_sub_knuth(uint64_t* u, const uint64_t* v, int n, uint64_t qhat){
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

static inline uint64_t add_back_knuth(uint64_t* u, const uint64_t* v, int n){
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

int u256_divmod(U256* q, U256* r, const U256* num, const U256* den){
    if(!q||!r||!num||!den) return -1;
    if(u256_is_zero(den)) return -1;

    if(u256_cmp(num,den)<0){ u256_set_ui(q,0); *r=*num; return 0; }

    int m = u256_num_limbs(num); // 1..4
    int n = u256_num_limbs(den); // 1..4

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

    unsigned s = (unsigned)u256_clz64(vraw[n-1]); // 0..63

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

    int qn = m - n + 1; // 1..4
    uint64_t qlimb[4] = {0,0,0,0};

    const uint64_t v1 = vnorm[n-1];
    const uint64_t v2 = vnorm[n-2];

    for(int j=qn-1;j>=0;j--){
        __uint128_t uj2 = ((__uint128_t)unorm[j+n] << 64) | (__uint128_t)unorm[j+n-1];
        uint64_t qhat = (uint64_t)(uj2 / v1);
        uint64_t rhat = (uint64_t)(uj2 % v1);

        // adjust qhat (Knuth D3)
        for(;;){
            __uint128_t left  = (__uint128_t)qhat * (__uint128_t)v2;
            __uint128_t right = ((__uint128_t)rhat << 64) | (__uint128_t)unorm[j+n-2];
            if(left <= right) break;
            qhat--;
            rhat += v1;
            if(rhat < v1) break; // overflow => stop
        }

        uint64_t* u_seg = &unorm[j];
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

void u256_set_sint_mod(U256* r, int64_t x, const U256* mod) {
    if (x >= 0) {
        u256_set_ui(r, (uint64_t)x);
        if (u256_cmp(r, mod) >= 0) {
            U256 q, rem;
            if (u256_divmod(&q, &rem, r, mod) == 0) *r = rem;
        }
        return;
    }

    uint64_t absv = (uint64_t)(-x);
    U256 a; u256_set_ui(&a, absv);

    if (u256_is_zero(&a)) {
        u256_set_ui(r, 0);
        return;
    }

    if (u256_cmp(&a, mod) >= 0) {
        U256 q, rem;
        if (u256_divmod(&q, &rem, &a, mod) == 0) a = rem;
    }

    if (u256_is_zero(&a)) {
        u256_set_ui(r, 0);
        return;
    }

    (void)u256_sub(r, mod, &a);
}

void u256_shl_2exp(U256* r, const U256* a, uint32_t k) {
    if (k >= 256) {
        r->limb[0] = r->limb[1] = r->limb[2] = r->limb[3] = 0;
        return;
    }

    U256 tmp;
    const U256* src = a;
    if (r == a) {
        tmp = *a;
        src = &tmp;
    }

    const uint32_t wordShift = k >> 6;
    const uint32_t bitShift  = k & 63;

    if (bitShift == 0) {
        for (int i = 3; i >= 0; --i) {
            int si = i - (int)wordShift;
            r->limb[i] = (si >= 0) ? src->limb[si] : 0;
        }
        return;
    }

    const uint32_t inv = 64 - bitShift;
    for (int i = 3; i >= 0; --i) {
        int si0 = i - (int)wordShift;
        int si1 = si0 - 1;

        uint64_t lo = (si0 >= 0) ? src->limb[si0] : 0;
        uint64_t hi = (si1 >= 0) ? src->limb[si1] : 0;

        r->limb[i] = (lo << bitShift) | (hi >> inv);
    }
}

#if !defined(__SIZEOF_INT128__)
#error "u256_mont_* requires __int128"
#endif

static inline uint64_t limbs_mul_1(uint64_t* r, const uint64_t* a, size_t n, uint64_t b) {
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t t = (__uint128_t)a[i] * b + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline uint64_t limbs_addmul_1(uint64_t* r, const uint64_t* a, size_t n, uint64_t b) {
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t t = (__uint128_t)a[i] * b + r[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    return (uint64_t)carry;
}

static inline void limbs_add_inplace(uint64_t* r, size_t rlen, const uint64_t* a, size_t alen) {
    __uint128_t carry = 0;
    size_t i = 0;

    for (; i < alen; i++) {
        __uint128_t t = (__uint128_t)r[i] + a[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
    }
    while (carry && i < rlen) {
        __uint128_t t = (__uint128_t)r[i] + carry;
        r[i] = (uint64_t)t;
        carry = t >> 64;
        i++;
    }
}

static inline int limbs_cmp(const uint64_t* a, const uint64_t* b, size_t n) {
    for (size_t i = n; i-- > 0;) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

static inline uint64_t limbs_sub_n(uint64_t* r, const uint64_t* a, const uint64_t* b, size_t n) {
    uint64_t borrow = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t ai = a[i];
        uint64_t bi = b[i];
        uint64_t ri = ai - bi - borrow;
        borrow = (borrow ? (ai <= bi) : (ai < bi));
        r[i] = ri;
    }
    return borrow;
}

void u256_mont_mul_4(uint64_t out[4], const uint64_t a[4], const uint64_t b[4],
                     const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[4] = limbs_mul_1(p0, b, 4, a[0]);

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);

    p1[4] = limbs_addmul_1(p1, b, 4, a[1]);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);

    p2[4] = limbs_addmul_1(p2, b, 4, a[2]);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);

    p3[4] = limbs_addmul_1(p3, b, 4, a[3]);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

void u256_mont_mul_4_u64(uint64_t out[4], const uint64_t a[4], uint64_t b,
                         const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[4] = limbs_mul_1(p0, a, 4, b);

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

void u256_mont_reduce_4(uint64_t out[4], const uint64_t in_mont[4],
                        const uint64_t mod5[5], uint64_t np) {
    constexpr size_t N = 5;

    uint64_t p0[N] = {0};
    uint64_t p1[N] = {0};
    uint64_t p2[N] = {0};
    uint64_t p3[N] = {0};

    p0[0] = in_mont[0];
    p0[1] = in_mont[1];
    p0[2] = in_mont[2];
    p0[3] = in_mont[3];
    p0[4] = 0;

    uint64_t m = np * p0[0];
    p1[1] = limbs_addmul_1(p0, mod5, N, m);
    limbs_add_inplace(p1, N, p0 + 1, N - 1);

    m = np * p1[0];
    p2[1] = limbs_addmul_1(p1, mod5, N, m);
    limbs_add_inplace(p2, N, p1 + 1, N - 1);

    m = np * p2[0];
    p3[1] = limbs_addmul_1(p2, mod5, N, m);
    limbs_add_inplace(p3, N, p2 + 1, N - 1);

    m = np * p3[0];
    (void)limbs_addmul_1(p3, mod5, N, m);

    out[0] = p3[1];
    out[1] = p3[2];
    out[2] = p3[3];
    out[3] = p3[4];

    if (limbs_cmp(out, mod5, 4) >= 0) {
        (void)limbs_sub_n(out, out, mod5, 4);
    }
}

void u256_scalar32_from_u64(uint8_t out[32], uint64_t v) {
    U256 x;
    u256_set_ui(&x, v);
    u256_export(out, &x);
}

int u256_scalar32_from_dec(uint8_t out[32], const char* dec) {
    U256 x;
    int rc = u256_set_str(&x, dec, 10);
    if (rc != 0) return rc;
    u256_export(out, &x);
    return 0;
}
