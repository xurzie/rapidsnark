#include "fr.hpp"
#include "mp.hpp"
#include <cstring>
#include <string>
#include <stdexcept>
#include <climits>

void Fr_toMP(mp_uint_t out, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);

    if (!(tmp.type & Fr_LONG)) {
        mp_set_mod(out, tmp.shortVal, Fr_q.longVal);
    } else {
        mp_copy(out, tmp.longVal);
    }
}

void Fr_fromMP(PFrElement pE, const mp_uint_t v) {
    if (mp_fits_int32(v)) {
        pE->type = Fr_SHORT;
        pE->shortVal = mp_get_int32(v);
    } else {
        pE->type = Fr_LONG;
        mp_copy(pE->longVal, v);
    }
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
    mp_div(q, rem, ma, mb);
    Fr_fromMP(r, q);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    mp_uint_t ma, mb;
    Fr_toMP(ma, a);
    Fr_toMP(mb, b);

    mp_uint_t q, rem;
    mp_div(q, rem, ma, mb);
    Fr_fromMP(r, rem);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    mp_uint_t mb;
    Fr_toMP(mb, b);

    mp_uint_t base;
    Fr_toMP(base, a);

    mp_uint_t res;
    mp_pow_mod(res, base, mb, Fr_q.longVal);

    Fr_fromMP(r, res);
}

void Fr_inv(PFrElement r, PFrElement a) {
    mp_uint_t base;
    Fr_toMP(base, a);

    mp_uint_t res;
    mp_inv_mod(res, base, Fr_q.longVal);

    Fr_fromMP(r, res);
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

RawFr::Element RawFr::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFr::set(Element& r, int value) {
    mp_set_mod(r.v, value, Fr_q.longVal);
    Fr_rawToMontgomery(r.v, r.v);
}

std::string RawFr::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, a.v);
    return mp_get_str(tmp.v, radix);
}

void RawFr::inv(Element& r, const Element& a) {
    Element t;
    Fr_rawFromMontgomery(t.v, a.v);
    mp_inv_mod(r.v, t.v, Fr_q.longVal);
    Fr_rawMMul(r.v, r.v, Fr_R2.longVal);
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

void RawFr::toMP(mp_uint_t r, const Element &a) {
    FrRawElement tmp;
    Fr_rawFromMontgomery(tmp, a.v);
    mp_copy(r, tmp);
}

void RawFr::fromMP(Element &a, const mp_uint_t r) {
    mp_copy(a.v, r);
    Fr_rawToMontgomery(a.v, a.v);
}

int RawFr::toRprBE(const Element& element, uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;

    mp_uint_t v;
    toMP(v, element);
    mp_export_be(data, v);

    return need;
}

int RawFr::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;

    mp_uint_t v;
    mp_import_be(v, data);
    fromMP(element, v);

    return need;
}

RawFr RawFr::field;
