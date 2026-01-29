#include "fq.hpp"
#include "u256.hpp"
#include <cstring>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <climits>

static bool initialized = false;

static bool Fq_init() {
    if (initialized) return false;
    initialized = true;
    return true;
}

void Fq_getModulusU256(U256 *q) {
    q->limb[0] = (uint64_t)Fq_q.longVal[0];
    q->limb[1] = (uint64_t)Fq_q.longVal[1];
    q->limb[2] = (uint64_t)Fq_q.longVal[2];
    q->limb[3] = (uint64_t)Fq_q.longVal[3];
}

 void Fq_fromU256(PFqElement pE, const U256 *v) {
    if (mp_fits_sint(v)) {
        pE->type = Fq_SHORT;
        pE->shortVal = (int32_t)v->limb[0];
        return;
    }
    pE->type = Fq_LONG;
    pE->longVal[0] = v->limb[0];
    pE->longVal[1] = v->limb[1];
    pE->longVal[2] = v->limb[2];
    pE->longVal[3] = v->limb[3];
}

void Fq_toU256(U256 *out, PFqElement pE) {
    FqElement tmp;
    Fq_toNormal(&tmp, pE);

    if (!(tmp.type & Fq_LONG)) {
        U256 mod; Fq_getModulusU256(&mod);
        mp_set_sint_mod(out, (int64_t)tmp.shortVal, &mod);
        return;
    }

    out->limb[0] = (uint64_t)tmp.longVal[0];
    out->limb[1] = (uint64_t)tmp.longVal[1];
    out->limb[2] = (uint64_t)tmp.longVal[2];
    out->limb[3] = (uint64_t)tmp.longVal[3];
}

static inline void fq_q_minus_2(uint8_t out_le[32]) {
    mp_limb_t t[4] = { (uint64_t)Fq_q.longVal[0], (uint64_t)Fq_q.longVal[1],
                       (uint64_t)Fq_q.longVal[2], (uint64_t)Fq_q.longVal[3] };

    (void)mp_sub_ui(t, t, 2u);
    mp_export(out_le, t);
}

static inline void Fq_toRawNormal(FqRawElement out, PFqElement a) {
    FqElement tmp;
    Fq_toNormal(&tmp, a);

    if (tmp.type & Fq_LONG) {
        out[0] = (uint64_t)tmp.longVal[0];
        out[1] = (uint64_t)tmp.longVal[1];
        out[2] = (uint64_t)tmp.longVal[2];
        out[3] = (uint64_t)tmp.longVal[3];
        return;
    }

    U256 mod;
    Fq_getModulusU256(&mod);

    U256 v;
    mp_set_sint_mod(&v, (int64_t)tmp.shortVal, &mod);

    out[0] = v.limb[0];
    out[1] = v.limb[1];
    out[2] = v.limb[2];
    out[3] = v.limb[3];
}

static inline void Fq_fromRawNormal(PFqElement out, const FqRawElement in) {
    if (in[1] == 0 && in[2] == 0 && in[3] == 0 && in[0] <= (uint64_t)INT_MAX) {
        out->type = Fq_SHORT;
        out->shortVal = (int32_t)in[0];
        return;
    }
    out->type = Fq_LONG;
    out->longVal[0] = (uint64_t)in[0];
    out->longVal[1] = (uint64_t)in[1];
    out->longVal[2] = (uint64_t)in[2];
    out->longVal[3] = (uint64_t)in[3];
}

static inline int bit_is_set_le(const uint8_t *s, int bit) {
    return (s[bit >> 3] & (uint8_t)(1u << (bit & 7))) != 0;
}

static void Fq_rawExpMont(FqRawElement out_mont, const FqRawElement base_mont, const uint8_t *exp_le, unsigned exp_size) {
    FqRawElement one_norm = {1u, 0u, 0u, 0u};
    FqRawElement one_mont;
    Fq_rawToMontgomery(one_mont, one_norm);

    bool oneFound = false;
    FqRawElement acc;
    FqRawElement copyBase;
    Fq_rawCopy(copyBase, base_mont);

    for (int i = (int)exp_size * 8 - 1; i >= 0; i--) {
        if (!oneFound) {
            if (!bit_is_set_le(exp_le, i)) continue;
            Fq_rawCopy(acc, copyBase);
            oneFound = true;
            continue;
        }
        Fq_rawMSquare(acc, acc);
        if (bit_is_set_le(exp_le, i)) {
            Fq_rawMMul(acc, acc, copyBase);
        }
    }

    if (!oneFound) {
        Fq_rawCopy(out_mont, one_mont);
        return;
    }
    Fq_rawCopy(out_mont, acc);
}

static char *mp_strdup_malloc(const std::string &s) {
    char *p = (char*)std::malloc(s.size() + 1);
    if (!p) return nullptr;
    std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

void Fq_str2element(PFqElement pE, char const *s, uint base) {
    U256 mod; Fq_getModulusU256(&mod);

    U256 v;
    if (mp_set_str_mod(&v, s, (int)base, &mod) != 0) {
        mp_set_ui(&v, 0);
    }
    Fq_fromU256(pE, &v);
}

char *Fq_element2str(PFqElement pE) {
    U256 v;
    Fq_toU256(&v, pE);
    return mp_strdup_malloc(mp_get_str(&v, 10));
}

void Fq_idiv(PFqElement r, PFqElement a, PFqElement b) {
    U256 ma, mb;
    Fq_toU256(&ma, a);
    Fq_toU256(&mb, b);

    U256 q, rem;
    if (mp_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fq_fromU256(r, &q);
}

void Fq_mod(PFqElement r, PFqElement a, PFqElement b) {
    U256 ma, mb;
    Fq_toU256(&ma, a);
    Fq_toU256(&mb, b);

    U256 q, rem;
    if (mp_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fq_fromU256(r, &rem);
}

void Fq_pow(PFqElement r, PFqElement a, PFqElement b) {
    U256 mb;
    Fq_toU256(&mb, b);

    uint8_t exp_le[32];
    mp_export(exp_le, &mb);

    FqRawElement base_norm;
    Fq_toRawNormal(base_norm, a);

    FqRawElement base_mont;
    Fq_rawToMontgomery(base_mont, base_norm);

    FqRawElement res_mont;
    Fq_rawExpMont(res_mont, base_mont, exp_le, (unsigned)sizeof(exp_le));

    FqRawElement res_norm;
    Fq_rawFromMontgomery(res_norm, res_mont);

    Fq_fromRawNormal(r, res_norm);
}

void Fq_inv(PFqElement r, PFqElement a) {
    uint8_t exp_le[32];
    fq_q_minus_2(exp_le);

    FqRawElement base_norm;
    Fq_toRawNormal(base_norm, a);

    FqRawElement base_mont;
    Fq_rawToMontgomery(base_mont, base_norm);

    FqRawElement res_mont;
    Fq_rawExpMont(res_mont, base_mont, exp_le, (unsigned)sizeof(exp_le));

    FqRawElement res_norm;
    Fq_rawFromMontgomery(res_norm, res_mont);

    Fq_fromRawNormal(r, res_norm);
}

void Fq_div(PFqElement r, PFqElement a, PFqElement b) {
    FqElement tmp;
    Fq_inv(&tmp, b);
    Fq_mul(r, a, &tmp);
}

void Fq_fail() {
    throw std::runtime_error("Fq error");
}

void Fq_longErr() {
    Fq_fail();
}

RawFq::RawFq() {
    Fq_init();
    set(fZero, 0);
    set(fOne, 1);
    neg(fNegOne, fOne);
}

RawFq::~RawFq() {}

void RawFq::fromString(Element& r, const std::string& s, uint32_t radix) {
    U256 mod; Fq_getModulusU256(&mod);
    U256 v;
    if (mp_set_str_mod(&v, s.c_str(), (int)radix, &mod) != 0) {
        mp_set_ui(&v, 0);
    }
    r.v[0] = v.limb[0];
    r.v[1] = v.limb[1];
    r.v[2] = v.limb[2];
    r.v[3] = v.limb[3];
    Fq_rawToMontgomery(r.v, r.v);
}

void RawFq::fromUI(Element& r, unsigned long int v) {
    U256 x;
    mp_set_ui(&x, (uint64_t)v);
    r.v[0] = x.limb[0];
    r.v[1] = x.limb[1];
    r.v[2] = x.limb[2];
    r.v[3] = x.limb[3];
    Fq_rawToMontgomery(r.v, r.v);
}

RawFq::Element RawFq::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFq::set(Element& r, int value) {
    U256 mod; Fq_getModulusU256(&mod);
    U256 v;
    mp_set_sint_mod(&v, (int64_t)value, &mod);

    r.v[0] = v.limb[0];
    r.v[1] = v.limb[1];
    r.v[2] = v.limb[2];
    r.v[3] = v.limb[3];
    Fq_rawToMontgomery(r.v, r.v);
}

std::string RawFq::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fq_rawFromMontgomery(tmp.v, a.v);

    U256 v;
    v.limb[0] = tmp.v[0];
    v.limb[1] = tmp.v[1];
    v.limb[2] = tmp.v[2];
    v.limb[3] = tmp.v[3];

    return mp_get_str(&v, (int)radix);
}

void RawFq::inv(Element& r, const Element& a) {
    uint8_t exp_le[32];
    fq_q_minus_2(exp_le);
    exp(r, a, exp_le, (unsigned)sizeof(exp_le));
}

void RawFq::div(Element& r, const Element& a, const Element& b) {
    Element tmp;
    inv(tmp, b);
    mul(r, a, tmp);
}

#define BIT_IS_SET(s, p) (s[(p)>>3] & (1 << ((p) & 0x7)))
void RawFq::exp(Element& r, const Element& base, uint8_t* scalar, unsigned int scalarSize) {
    bool oneFound = false;
    Element copyBase;
    copy(copyBase, base);

    for (int i = (int)scalarSize * 8 - 1; i >= 0; i--) {
        if (!oneFound) {
            if (!BIT_IS_SET(scalar, i)) continue;
            copy(r, copyBase);
            oneFound = true;
            continue;
        }
        square(r, r);
        if (BIT_IS_SET(scalar, i)) {
            mul(r, r, copyBase);
        }
    }
    if (!oneFound) {
        copy(r, fOne);
    }
}

int RawFq::toRprBE(const Element& element, uint8_t* data, int bytes) {
    const int need = Fq_N64 * 8;
    if (bytes < need) return -need;

    Element tmp;
    Fq_rawFromMontgomery(tmp.v, element.v);

    U256 v;
    v.limb[0] = tmp.v[0];
    v.limb[1] = tmp.v[1];
    v.limb[2] = tmp.v[2];
    v.limb[3] = tmp.v[3];

    mp_export_be(data, &v);
    return need;
}

int RawFq::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fq_N64 * 8;
    if (bytes < need) return -need;

    U256 v;
    mp_import_be(&v, data);

    element.v[0] = v.limb[0];
    element.v[1] = v.limb[1];
    element.v[2] = v.limb[2];
    element.v[3] = v.limb[3];
    Fq_rawToMontgomery(element.v, element.v);
    return need;
}

static bool init = Fq_init();
RawFq RawFq::field;
