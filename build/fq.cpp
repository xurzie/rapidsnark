#include "fq.hpp"
#include "mp.hpp"
#include <cstring>
#include <string>
#include <stdexcept>
#include <climits>

void Fq_fromMP(PFqElement pE, const uint64_t *v) {
    if (mp_fits_int32(v)) {
        pE->type = Fq_SHORT;
        pE->shortVal = (int32_t)v[0];
        return;
    }
    pE->type = Fq_LONG;
    mp_copy(pE->longVal, v);
}

void Fq_toMP(uint64_t *out, PFqElement pE) {
    FqElement tmp;
    Fq_toNormal(&tmp, pE);
    if (!(tmp.type & Fq_LONG)) {
        mp_set_mod(out, (int64_t)tmp.shortVal, Fq_q.longVal);
        return;
    }
    mp_copy(out, tmp.longVal);
}

static inline void fq_q_minus_2(uint8_t out[MP_N]) {
    mp_uint_t t;
    mp_copy(t, Fq_q.longVal);
    mp_sub(t, t, 2u);
    mp_export(out, t);
}

static inline void Fq_toRawNormal(FqRawElement out, PFqElement a) {
    FqElement tmp;
    Fq_toNormal(&tmp, a);
    if (tmp.type & Fq_LONG) {
        mp_copy(out, tmp.longVal);
        return;
    }
    mp_uint_t v;
    mp_set_mod(v, (int64_t)tmp.shortVal, Fq_q.longVal);
    mp_copy(out, v);
}

static inline void Fq_fromRawNormal(PFqElement out, const FqRawElement in) {
    if (in[1] == 0 && in[2] == 0 && in[3] == 0 && in[0] <= (uint64_t)INT_MAX) {
        out->type = Fq_SHORT;
        out->shortVal = (int32_t)in[0];
        return;
    }
    out->type = Fq_LONG;
    mp_copy(out->longVal, in);
}

static inline int bit_is_set(const uint8_t *s, int bit) {
    return (s[bit >> 3] & (uint8_t)(1u << (bit & 7))) != 0;
}

static void Fq_rawExpMont(FqRawElement out_mont, const FqRawElement base_mont, const uint8_t *exp, unsigned exp_size) {
    FqRawElement one_norm;
    mp_set(one_norm, 1u);
    FqRawElement one_mont;
    Fq_rawToMontgomery(one_mont, one_norm);

    bool oneFound = false;
    FqRawElement acc;
    FqRawElement copyBase;
    Fq_rawCopy(copyBase, base_mont);

    for (int i = (int)exp_size * 8 - 1; i >= 0; i--) {
        if (!oneFound) {
            if (!bit_is_set(exp, i)) continue;
            Fq_rawCopy(acc, copyBase);
            oneFound = true;
            continue;
        }
        Fq_rawMSquare(acc, acc);
        if (bit_is_set(exp, i)) {
            Fq_rawMMul(acc, acc, copyBase);
        }
    }

    if (!oneFound) {
        Fq_rawCopy(out_mont, one_mont);
        return;
    }
    Fq_rawCopy(out_mont, acc);
}

void Fq_str2element(PFqElement pE, char const *s, uint base) {
    mp_uint_t v;
    if (!mp_set_mod(v, s, base, Fq_q.longVal)) {
        mp_set(v, 0);
    }
    Fq_fromMP(pE, v);
}

std::string Fq_element2str(PFqElement pE, uint32_t base) {
    mp_uint_t v;
    Fq_toMP(v, pE);
    return mp_get_str(v, base);
}

void Fq_idiv(PFqElement r, PFqElement a, PFqElement b) {
    mp_uint_t ma, mb;
    Fq_toMP(ma, a);
    Fq_toMP(mb, b);

    mp_uint_t q, rem;
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fq_fromMP(r, q);
}

void Fq_mod(PFqElement r, PFqElement a, PFqElement b) {
    mp_uint_t ma, mb;
    Fq_toMP(ma, a);
    Fq_toMP(mb, b);

    mp_uint_t q, rem;
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fq_fromMP(r, rem);
}

void Fq_pow(PFqElement r, PFqElement a, PFqElement b) {
    mp_uint_t mb;
    Fq_toMP(mb, b);

    uint8_t exp[MP_N];
    mp_export(exp, mb);

    FqRawElement base_norm;
    Fq_toRawNormal(base_norm, a);

    FqRawElement res_norm;
    mp_powm(res_norm, base_norm, exp, (unsigned)sizeof(exp), Fq_q.longVal);

    Fq_fromRawNormal(r, res_norm);
}

void Fq_inv(PFqElement r, PFqElement a) {
    FqRawElement base_norm;
    Fq_toRawNormal(base_norm, a);

    FqRawElement res_norm;
    mp_invert(res_norm, base_norm, Fq_q.longVal);

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
    set(fZero, 0);
    set(fOne, 1);
    neg(fNegOne, fOne);
}

RawFq::~RawFq() {}

void RawFq::fromString(Element& r, const std::string& s, uint32_t radix) {
    if (!mp_set_mod(r.v, s.c_str(), radix, Fq_q.longVal)) {
        mp_set(r.v, 0);
    }
    Fq_rawToMontgomery(r.v, r.v);
}

void RawFq::fromUI(Element& r, unsigned long int v) {
    mp_set(r.v, v);
    Fq_rawToMontgomery(r.v, r.v);
}

void RawFq::toMP(mp_uint_t r, const Element &a) {
    FqRawElement tmp;
    Fq_rawFromMontgomery(tmp, a.v);
    mp_copy(r, tmp);
}

void RawFq::fromMP(Element &a, const mp_uint_t r) {
    mp_copy(a.v, r);
    Fq_rawToMontgomery(a.v, a.v);
}

RawFq::Element RawFq::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFq::set(Element& r, int value) {
    mp_set_mod(r.v, (int64_t)value, Fq_q.longVal);
    Fq_rawToMontgomery(r.v, r.v);
}

std::string RawFq::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fq_rawFromMontgomery(tmp.v, a.v);
    mp_uint_t v;
    mp_copy(v, tmp.v);
    return mp_get_str(v, (int)radix);
}

void RawFq::inv(Element& r, const Element& a) {
    mp_uint_t an;
    toMP(an, a);

    mp_uint_t invn;
    mp_invert(invn, an, Fq_q.longVal);

    fromMP(r, invn);
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
    mp_uint_t v;
    mp_copy(v, tmp.v);
    mp_export_be(data, v);
    return need;
}

int RawFq::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fq_N64 * 8;
    if (bytes < need) return -need;
    mp_uint_t v;
    mp_import_be(v, data);
    mp_copy(element.v, v);
    Fq_rawToMontgomery(element.v, element.v);
    return need;
}

RawFq RawFq::field;
