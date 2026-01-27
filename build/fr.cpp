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

void Fr_getModulusU256(U256 *q) {
    q->limb[0] = (uint64_t)Fr_q.longVal[0];
    q->limb[1] = (uint64_t)Fr_q.longVal[1];
    q->limb[2] = (uint64_t)Fr_q.longVal[2];
    q->limb[3] = (uint64_t)Fr_q.longVal[3];
}

void Fr_fromU256(PFrElement pE, const U256 *v) {
    if (mp_fits_sint(v)) {
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

void Fr_toU256(U256 *out, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);

    if (!(tmp.type & Fr_LONG)) {
        U256 mod; Fr_getModulusU256(&mod);
        mp_set_sint_mod(out, (int64_t)tmp.shortVal, &mod);
        return;
    }

    out->limb[0] = (uint64_t)tmp.longVal[0];
    out->limb[1] = (uint64_t)tmp.longVal[1];
    out->limb[2] = (uint64_t)tmp.longVal[2];
    out->limb[3] = (uint64_t)tmp.longVal[3];
}

static inline void load_u256(U256 *out, const FrRawElement in) {
    out->limb[0] = in[0];
    out->limb[1] = in[1];
    out->limb[2] = in[2];
    out->limb[3] = in[3];
}

static inline void store_u256(FrRawElement out, const U256 *a) {
    out[0] = a->limb[0];
    out[1] = a->limb[1];
    out[2] = a->limb[2];
    out[3] = a->limb[3];
}

static char *mp_strdup_malloc(const std::string &s) {
    char *p = (char*)std::malloc(s.size() + 1);
    if (!p) return nullptr;
    std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

void Fr_str2element(PFrElement pE, char const* s, uint base) {
    U256 mod; Fr_getModulusU256(&mod);

    U256 v;
    if (mp_set_str_mod(&v, s, (int)base, &mod) != 0) {
        mp_set_ui(&v, 0);
    }
    Fr_fromU256(pE, &v);
}

char *Fr_element2str(PFrElement pE) {
    U256 v;
    Fr_toU256(&v, pE);
    return mp_strdup_malloc(mp_get_str(&v, 10));
}

void Fr_idiv(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    Fr_toU256(&ma, a);
    Fr_toU256(&mb, b);

    U256 q, rem;
    if (mp_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromU256(r, &q);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    Fr_toU256(&ma, a);
    Fr_toU256(&mb, b);

    U256 q, rem;
    if (mp_divmod(&q, &rem, &ma, &mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromU256(r, &rem);
}

static inline void fr_q_minus_2(uint8_t out[32]) {
    U256 mod; Fr_getModulusU256(&mod);
    U256 e;
    (void)mp_sub_ui(&e, &mod, 2);
    mp_export(out, &e);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    U256 ma, mb;
    Fr_toU256(&ma, a);
    Fr_toU256(&mb, b);

    uint8_t exp_le[32];
    mp_export(exp_le, &mb);

    FrRawElement base_norm;
    store_u256(base_norm, &ma);

    FrRawElement base_mont;
    Fr_rawCopy(base_mont, base_norm);
    Fr_rawToMontgomery(base_mont, base_mont);

    RawFr::Element B{}, R{};
    Fr_rawCopy(B.v, base_mont);

    RawFr::field.exp(R, B, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, R.v);

    U256 out;
    load_u256(&out, res_norm);
    Fr_fromU256(r, &out);
}

void Fr_inv(PFrElement r, PFrElement a) {
    U256 ma;
    Fr_toU256(&ma, a);

    uint8_t exp_le[32];
    fr_q_minus_2(exp_le);

    FrRawElement base_norm;
    store_u256(base_norm, &ma);

    FrRawElement base_mont;
    Fr_rawCopy(base_mont, base_norm);
    Fr_rawToMontgomery(base_mont, base_mont);

    RawFr::Element B{}, R{};
    Fr_rawCopy(B.v, base_mont);

    RawFr::field.exp(R, B, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, R.v);

    U256 out;
    load_u256(&out, res_norm);
    Fr_fromU256(r, &out);
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
    U256 mod; Fr_getModulusU256(&mod);
    U256 v;
    if (mp_set_str_mod(&v, s.c_str(), (int)radix, &mod) != 0) {
        mp_set_ui(&v, 0);
    }
    r.v[0] = v.limb[0];
    r.v[1] = v.limb[1];
    r.v[2] = v.limb[2];
    r.v[3] = v.limb[3];
    Fr_rawToMontgomery(r.v, r.v);
}

void RawFr::fromUI(Element& r, unsigned long int v) {
    U256 x;
    mp_set_ui(&x, (uint64_t)v);
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
    U256 mod; Fr_getModulusU256(&mod);
    U256 v;
    mp_set_sint_mod(&v, (int64_t)value, &mod);

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

    return mp_get_str(&v, (int)radix);
}

void RawFr::inv(Element& r, const Element& a) {
    uint8_t exp_le[32];
    fr_q_minus_2(exp_le);
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

    mp_export_be(data, &v);
    return need;
}

int RawFr::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;

    U256 v;
    mp_import_be(&v, data);

    element.v[0] = v.limb[0];
    element.v[1] = v.limb[1];
    element.v[2] = v.limb[2];
    element.v[3] = v.limb[3];
    Fr_rawToMontgomery(element.v, element.v);
    return need;
}

static bool init = Fr_init();
RawFr RawFr::field;
