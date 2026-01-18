#include "fr.hpp"
#include "u256.hpp"

#include <cstring>
#include <cstdlib>
#include <string>
#include <stdexcept>

static bool initialized = false;

static bool Fr_init() {
    if (initialized) return false;
    initialized = true;
    return true;
}

static inline void fr_get_modulus(U256* q) {
    q->limb[0] = (uint64_t)Fr_q.longVal[0];
    q->limb[1] = (uint64_t)Fr_q.longVal[1];
    q->limb[2] = (uint64_t)Fr_q.longVal[2];
    q->limb[3] = (uint64_t)Fr_q.longVal[3];
}

static inline void fr_element_from_u256(PFrElement pE, const U256* v) {
    if (u256_fits_sint(v)) {
        pE->type = Fr_SHORT;
        pE->shortVal = (int32_t)v->limb[0];
        return;
    }
    pE->type = Fr_LONG;
    pE->longVal[0] = v->limb[0];
    pE->longVal[1] = v->limb[1];
    pE->longVal[2] = v->limb[2];
    pE->longVal[3] = v->limb[3];
}

static inline void fr_element_to_u256(U256* out, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);

    if (!(tmp.type & Fr_LONG)) {
        U256 mod; fr_get_modulus(&mod);
        u256_set_sint_mod(out, (int64_t)tmp.shortVal, &mod);
        return;
    }

    out->limb[0] = (uint64_t)tmp.longVal[0];
    out->limb[1] = (uint64_t)tmp.longVal[1];
    out->limb[2] = (uint64_t)tmp.longVal[2];
    out->limb[3] = (uint64_t)tmp.longVal[3];
}

static inline void u256_to_raw(FrRawElement out, const U256* a) {
    out[0] = a->limb[0];
    out[1] = a->limb[1];
    out[2] = a->limb[2];
    out[3] = a->limb[3];
}

static inline void u256_from_raw(U256* out, const FrRawElement in) {
    out->limb[0] = in[0];
    out->limb[1] = in[1];
    out->limb[2] = in[2];
    out->limb[3] = in[3];
}

void Fr_str2element(PFrElement pE, char const* s, uint base) {
    U256 mod; fr_get_modulus(&mod);

    U256 v;
    if (u256_set_str_mod(&v, s, (int)base, &mod) != 0) {
        u256_set_ui(&v, 0);
    }
    fr_element_from_u256(pE, &v);
}

char* Fr_element2str(PFrElement pE) {
    U256 v;
    fr_element_to_u256(&v, pE);
    return u256_get_str_alloc(&v, 10);
}

void Fr_idiv(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    fr_element_to_u256(&ma, a);
    fr_element_to_u256(&mb, b);

    U256 q, rem;
    if (u256_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    fr_element_from_u256(r, &q);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    fr_element_to_u256(&ma, a);
    fr_element_to_u256(&mb, b);

    U256 q, rem;
    if (u256_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    fr_element_from_u256(r, &rem);
}

static inline void fr_q_minus_2_le(uint8_t out_le[32]) {
    U256 mod; fr_get_modulus(&mod);
    U256 e;
    (void)u256_sub_ui(&e, &mod, 2);
    u256_export(out_le, &e);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    fr_element_to_u256(&ma, a);
    fr_element_to_u256(&mb, b);

    uint8_t exp_le[32];
    u256_export(exp_le, &mb);

    FrRawElement base_norm;
    u256_to_raw(base_norm, &ma);

    FrRawElement base_mont;
    Fr_rawCopy(base_mont, base_norm);
    Fr_rawToMontgomery(base_mont, base_mont);

    RawFr::Element B{}, R{};
    Fr_rawCopy(B.v, base_mont);

    RawFr::field.exp(R, B, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, R.v);

    U256 out;
    u256_from_raw(&out, res_norm);
    fr_element_from_u256(r, &out);
}

void Fr_inv(PFrElement r, PFrElement a) {
    U256 ma;
    fr_element_to_u256(&ma, a);

    uint8_t exp_le[32];
    fr_q_minus_2_le(exp_le);

    FrRawElement base_norm;
    u256_to_raw(base_norm, &ma);

    FrRawElement base_mont;
    Fr_rawCopy(base_mont, base_norm);
    Fr_rawToMontgomery(base_mont, base_mont);

    RawFr::Element B{}, R{};
    Fr_rawCopy(B.v, base_mont);

    RawFr::field.exp(R, B, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, R.v);

    U256 out;
    u256_from_raw(&out, res_norm);
    fr_element_from_u256(r, &out);
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
    Fr_init();
    set(fZero, 0);
    set(fOne, 1);
    neg(fNegOne, fOne);
}

RawFr::~RawFr() {}

void RawFr::fromString(Element& r, const std::string& s, uint32_t radix) {
    U256 mod; fr_get_modulus(&mod);
    U256 v;
    if (u256_set_str_mod(&v, s.c_str(), (int)radix, &mod) != 0) {
        u256_set_ui(&v, 0);
    }
    r.v[0] = v.limb[0];
    r.v[1] = v.limb[1];
    r.v[2] = v.limb[2];
    r.v[3] = v.limb[3];
    Fr_rawToMontgomery(r.v, r.v);
}

void RawFr::fromUI(Element& r, unsigned long int v) {
    U256 x;
    u256_set_ui(&x, (uint64_t)v);
    r.v[0] = x.limb[0];
    r.v[1] = x.limb[1];
    r.v[2] = x.limb[2];
    r.v[3] = x.limb[3];
    Fr_rawToMontgomery(r.v, r.v);
}

RawFr::Element RawFr::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFr::set(Element& r, int value) {
    U256 mod; fr_get_modulus(&mod);
    U256 v;
    u256_set_sint_mod(&v, (int64_t)value, &mod);

    r.v[0] = v.limb[0];
    r.v[1] = v.limb[1];
    r.v[2] = v.limb[2];
    r.v[3] = v.limb[3];
    Fr_rawToMontgomery(r.v, r.v);
}

std::string RawFr::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, a.v);

    U256 v;
    v.limb[0] = tmp.v[0];
    v.limb[1] = tmp.v[1];
    v.limb[2] = tmp.v[2];
    v.limb[3] = tmp.v[3];

    char* s = u256_get_str_alloc(&v, (int)radix);
    if (!s) return std::string();
    std::string out(s);
    std::free(s);
    return out;
}

void RawFr::inv(Element& r, const Element& a) {
    uint8_t exp_le[32];
    fr_q_minus_2_le(exp_le);
    exp(r, a, exp_le, (unsigned)sizeof(exp_le));
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

    U256 v;
    v.limb[0] = tmp.v[0];
    v.limb[1] = tmp.v[1];
    v.limb[2] = tmp.v[2];
    v.limb[3] = tmp.v[3];

    u256_export_be(data, &v);
    return need;
}

int RawFr::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;

    U256 v;
    u256_import_be(&v, data);

    element.v[0] = v.limb[0];
    element.v[1] = v.limb[1];
    element.v[2] = v.limb[2];
    element.v[3] = v.limb[3];
    Fr_rawToMontgomery(element.v, element.v);
    return need;
}

static bool init = Fr_init();
RawFr RawFr::field;
