#include "fr.hpp"
#include "mp.hpp"

#include <cstring>
#include <string>
#include <stdexcept>
#include <climits>

void Fr_fromMP(PFrElement pE, const uint64_t *v) {
    if (mp_fits_int32(v)) {
        pE->type = Fr_SHORT;
        pE->shortVal = (int32_t)v[0];
        return;
    }
    pE->type = Fr_LONG;
    mp_copy(pE->longVal, v);
}

void Fr_toMP(uint64_t *out, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);

    if (!(tmp.type & Fr_LONG)) {
        mp_set_mod(out, (int64_t)tmp.shortVal, Fr_q.longVal);
        return;
    }
    mp_copy(out, tmp.longVal);
}

static inline void fr_q_minus_2(uint8_t out_le[MP_N]) {
    mp_uint_t t;
    mp_copy(t, Fr_q.longVal);
    mp_sub(t, t, 2u);
    mp_export(out_le, t);
}

static inline void Fr_toRawNormal(FrRawElement out, PFrElement a) {
    FrElement tmp;
    Fr_toNormal(&tmp, a);
    if (tmp.type & Fr_LONG) {
        mp_copy(out, tmp.longVal);
        return;
    }
    mp_uint_t v;
    mp_set_mod(v, (int64_t)tmp.shortVal, Fr_q.longVal);
    mp_copy(out, v);
}

static inline void Fr_fromRawNormal(PFrElement out, const FrRawElement in) {
    if (in[1] == 0 && in[2] == 0 && in[3] == 0 && in[0] <= (uint64_t)INT_MAX) {
        out->type = Fr_SHORT;
        out->shortVal = (int32_t)in[0];
        return;
    }
    out->type = Fr_LONG;
    mp_copy(out->longVal, in);
}

static inline int bit_is_set_le(const uint8_t *s, int bit) {
    return (s[bit >> 3] & (uint8_t)(1u << (bit & 7))) != 0;
}

static void Fr_rawExpMont(FrRawElement out_mont, const FrRawElement base_mont, const uint8_t *exp, unsigned exp_size) {
    FrRawElement one_norm;
    mp_set(one_norm, 1u);
    FrRawElement one_mont;
    Fr_rawToMontgomery(one_mont, one_norm);

    bool oneFound = false;
    FrRawElement acc;
    FrRawElement copyBase;
    Fr_rawCopy(copyBase, base_mont);

    for (int i = (int)exp_size * 8 - 1; i >= 0; i--) {
        if (!oneFound) {
            if (!bit_is_set_le(exp, i)) continue;
            Fr_rawCopy(acc, copyBase);
            oneFound = true;
            continue;
        }
        Fr_rawMSquare(acc, acc);
        if (bit_is_set_le(exp, i)) {
            Fr_rawMMul(acc, acc, copyBase);
        }
    }

    if (!oneFound) {
        Fr_rawCopy(out_mont, one_mont);
        return;
    }
    Fr_rawCopy(out_mont, acc);
}

void Fr_str2element(PFrElement pE, char const* s, uint base) {
    mp_uint_t v;
    if (!mp_set_mod(v, s, base, Fr_q.longVal)) {
        mp_set(v, 0);
    }
    Fr_fromMP(pE, v);
}

std::string Fr_element2str(PFrElement pE, uint32_t base) {
    mp_uint_t v;
    Fr_toMP(v, pE);
    return mp_get_str(v, base);
}

void Fr_idiv(PFrElement r, PFrElement a, PFrElement b) {
    mp_uint_t ma, mb;
    Fr_toMP(ma, a);
    Fr_toMP(mb, b);

    mp_uint_t q, rem;
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromMP(r, q);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    mp_uint_t ma, mb;
    Fr_toMP(ma, a);
    Fr_toMP(mb, b);

    mp_uint_t q, rem;
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromMP(r, rem);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    mp_uint_t mb;
    Fr_toMP(mb, b);

    uint8_t exp[MP_N];
    mp_export(exp, mb);

    FrRawElement base_norm;
    Fr_toRawNormal(base_norm, a);

    FrRawElement res_norm;
    mp_powm(res_norm, base_norm, exp, (unsigned)sizeof(exp), Fr_q.longVal);

    Fr_fromRawNormal(r, res_norm);
}

void Fr_inv(PFrElement r, PFrElement a) {
    FrRawElement base_norm;
    Fr_toRawNormal(base_norm, a);

    FrRawElement res_norm;
    mp_invert(res_norm, base_norm, Fr_q.longVal);

    Fr_fromRawNormal(r, res_norm);
}

void Fr_div(PFrElement r, PFrElement a, PFrElement b) {
    FrElement invb;
    Fr_inv(&invb, b);
    Fr_mul(r, a, &invb);
}

void Fr_fail() {
    throw std::runtime_error("Fr error");
}

void Fr_longErr() {
    Fr_fail();
}

RawFr::RawFr() {
    set(fZero, 0);
    set(fOne, 1);
    neg(fNegOne, fOne);
}

RawFr::~RawFr() {}

void RawFr::fromString(Element& r, const std::string& s, uint32_t radix) {
    if (!mp_set_mod(r.v, s.c_str(), radix, Fr_q.longVal)) {
        mp_set(r.v, 0);
    }
    Fr_rawToMontgomery(r.v, r.v);
}

void RawFr::fromUI(Element& r, unsigned long int v) {
    mp_set(r.v, (uint64_t)v);
    Fr_rawToMontgomery(r.v, r.v);
}

void RawFr::toMP(mp_uint_t r, const Element &a) {
    FrRawElement tmp;
    Fr_rawFromMontgomery(tmp, a.v);
    mp_copy(r, tmp);
}

void RawFr::fromMP(Element &a, const mp_uint_t r) {
    mp_copy(a.v, r);
    Fr_rawToMontgomery(a.v, a.v);
}

RawFr::Element RawFr::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFr::set(Element& r, int value) {
    mp_set_mod(r.v, (int64_t)value, Fr_q.longVal);
    Fr_rawToMontgomery(r.v, r.v);
}

std::string RawFr::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, a.v);
    mp_uint_t v;
    mp_copy(v, tmp.v);
    return mp_get_str(v, (int)radix);
}

void RawFr::inv(Element& r, const Element& a) {
    mp_uint_t an;
    toMP(an, a);

    mp_uint_t invn;
    mp_invert(invn, an, Fr_q.longVal);

    fromMP(r, invn);
}

void RawFr::div(Element& r, const Element& a, const Element& b) {
    Element tmp;
    inv(tmp, b);
    mul(r, a, tmp);
}

#define BIT_IS_SET(s, p) (s[(p)>>3] & (1 << ((p) & 0x7)))
void RawFr::exp(Element& r, const Element& base, uint8_t* scalar, unsigned int scalarSize) {
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

int RawFr::toRprBE(const Element& element, uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, element.v);
    mp_uint_t v;
    mp_copy(v, tmp.v);
    mp_export_be(data, v);
    return need;
}

int RawFr::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;
    mp_uint_t v;
    mp_import_be(v, data);
    mp_copy(element.v, v);
    Fr_rawToMontgomery(element.v, element.v);
    return need;
}

RawFr RawFr::field;
