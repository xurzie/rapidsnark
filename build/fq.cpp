#include "fq.hpp"
#include "mp.hpp"
#include <cstring>
#include <string>
#include <stdexcept>
#include <climits>

void Fq_toMP(mp_uint_t out, PFqElement pE) {
    FqElement tmp;
    Fq_toNormal(&tmp, pE);

    if (!(tmp.type & Fq_LONG)) {
        mp_set_mod(out, tmp.shortVal, Fq_q.longVal);
    } else {
        mp_copy(out, tmp.longVal);
    }
}

void Fq_fromMP(PFqElement pE, const mp_uint_t v) {
    if (mp_fits_int32(v)) {
        pE->type = Fq_SHORT;
        pE->shortVal = mp_get_int32(v);
    } else {
        pE->type = Fq_LONG;
        mp_copy(pE->longVal, v);
    }
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
    mp_div(q, rem, ma, mb);
    Fq_fromMP(r, q);
}

void Fq_mod(PFqElement r, PFqElement a, PFqElement b) {
    mp_uint_t ma, mb;
    Fq_toMP(ma, a);
    Fq_toMP(mb, b);

    mp_uint_t q, rem;
    mp_div(q, rem, ma, mb);
    Fq_fromMP(r, rem);
}

void Fq_pow(PFqElement r, PFqElement a, PFqElement b) {
    mp_uint_t mb;
    Fq_toMP(mb, b);

    mp_uint_t base;
    Fq_toMP(base, a);

    mp_uint_t res;
    mp_pow_mod(res, base, mb, Fq_q.longVal);

    Fq_fromMP(r, res);
}

void Fq_inv(PFqElement r, PFqElement a) {
    mp_uint_t base;
    Fq_toMP(base, a);

    mp_uint_t res;
    mp_inv_mod(res, base, Fq_q.longVal);

    Fq_fromMP(r, res);
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

RawFq::Element RawFq::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFq::set(Element& r, int value) {
    mp_set_mod(r.v, value, Fq_q.longVal);
    Fq_rawToMontgomery(r.v, r.v);
}

std::string RawFq::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fq_rawFromMontgomery(tmp.v, a.v);
    return mp_get_str(tmp.v, radix);
}

void RawFq::inv(Element& r, const Element& a) {
    Element t;
    Fq_rawFromMontgomery(t.v, a.v);
    mp_inv_mod(r.v, t.v, Fq_q.longVal);
    Fq_rawMMul(r.v, r.v, Fq_R2.longVal);
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

void RawFq::toMP(mp_uint_t r, const Element &a) {
    FqRawElement tmp;
    Fq_rawFromMontgomery(tmp, a.v);
    mp_copy(r, tmp);
}

void RawFq::fromMP(Element &a, const mp_uint_t r) {
    mp_copy(a.v, r);
    Fq_rawToMontgomery(a.v, a.v);
}

int RawFq::toRprBE(const Element& element, uint8_t* data, int bytes) {
    const int need = Fq_N64 * 8;
    if (bytes < need) return -need;

    mp_uint_t v;
    toMP(v, element);
    mp_export_be(data, v);

    return need;
}

int RawFq::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fq_N64 * 8;
    if (bytes < need) return -need;

    mp_uint_t v;
    mp_import_be(v, data);
    fromMP(element, v);

    return need;
}

RawFq RawFq::field;
