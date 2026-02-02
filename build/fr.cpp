#include "fr.hpp"
#include "mp.hpp"

#include <cstring>
#include <string>
#include <stdexcept>
#include <climits>

static bool initialized = false;

static bool Fr_init() {
    if (initialized) return false;
    initialized = true;
    return true;
}

void Fr_getModulusU256(mp_limb_t *q) {
    q[0] = (uint64_t)Fr_q.longVal[0];
    q[1] = (uint64_t)Fr_q.longVal[1];
    q[2] = (uint64_t)Fr_q.longVal[2];
    q[3] = (uint64_t)Fr_q.longVal[3];
}

void Fr_fromU256(PFrElement pE, const mp_limb_t *v) {
    if (mp_fits_sint(v)) {
        pE->type = Fr_SHORT;
        pE->shortVal = (int32_t)v[0];
        return;
    }
    pE->type = Fr_LONG;
    pE->longVal[0] = v[0];
    pE->longVal[1] = v[1];
    pE->longVal[2] = v[2];
    pE->longVal[3] = v[3];
}

void Fr_toU256(mp_limb_t *out, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);

    if (!(tmp.type & Fr_LONG)) {
        mp_limb_t mod[4]; Fr_getModulusU256(mod);
        mp_set_sint_mod(out, (int64_t)tmp.shortVal, mod);
        return;
    }

    out[0] = (uint64_t)tmp.longVal[0];
    out[1] = (uint64_t)tmp.longVal[1];
    out[2] = (uint64_t)tmp.longVal[2];
    out[3] = (uint64_t)tmp.longVal[3];
}

static inline void fr_q_minus_2(uint8_t out_le[32]) {
    mp_limb_t t[4] = { (uint64_t)Fr_q.longVal[0], (uint64_t)Fr_q.longVal[1],
                       (uint64_t)Fr_q.longVal[2], (uint64_t)Fr_q.longVal[3] };

    mp_sub(t, t, 2u);
    mp_export(out_le, t);
}

static inline void Fr_toRawNormal(FrRawElement out, PFrElement a) {
    FrElement tmp;
    Fr_toNormal(&tmp, a);

    if (tmp.type & Fr_LONG) {
        out[0] = (uint64_t)tmp.longVal[0];
        out[1] = (uint64_t)tmp.longVal[1];
        out[2] = (uint64_t)tmp.longVal[2];
        out[3] = (uint64_t)tmp.longVal[3];
        return;
    }

    mp_limb_t mod[4];
    Fr_getModulusU256(mod);

    mp_limb_t v[4];
    mp_set_sint_mod(v, (int64_t)tmp.shortVal, mod);

    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
    out[3] = v[3];
}

static inline void Fr_fromRawNormal(PFrElement out, const FrRawElement in) {
    if (in[1] == 0 && in[2] == 0 && in[3] == 0 && in[0] <= (uint64_t)INT_MAX) {
        out->type = Fr_SHORT;
        out->shortVal = (int32_t)in[0];
        return;
    }
    out->type = Fr_LONG;
    out->longVal[0] = (uint64_t)in[0];
    out->longVal[1] = (uint64_t)in[1];
    out->longVal[2] = (uint64_t)in[2];
    out->longVal[3] = (uint64_t)in[3];
}

static inline int bit_is_set_le(const uint8_t *s, int bit) {
    return (s[bit >> 3] & (uint8_t)(1u << (bit & 7))) != 0;
}

static void Fr_rawExpMont(FrRawElement out_mont, const FrRawElement base_mont, const uint8_t *exp_le, unsigned exp_size) {
    FrRawElement one_norm = {1u, 0u, 0u, 0u};
    FrRawElement one_mont;
    Fr_rawToMontgomery(one_mont, one_norm);

    bool oneFound = false;
    FrRawElement acc;
    FrRawElement copyBase;
    Fr_rawCopy(copyBase, base_mont);

    for (int i = (int)exp_size * 8 - 1; i >= 0; i--) {
        if (!oneFound) {
            if (!bit_is_set_le(exp_le, i)) continue;
            Fr_rawCopy(acc, copyBase);
            oneFound = true;
            continue;
        }
        Fr_rawMSquare(acc, acc);
        if (bit_is_set_le(exp_le, i)) {
            Fr_rawMMul(acc, acc, copyBase);
        }
    }

    if (!oneFound) {
        Fr_rawCopy(out_mont, one_mont);
        return;
    }
    Fr_rawCopy(out_mont, acc);
}

static char *mp_strdup_malloc(const std::string &s) {
    char *p = (char*)std::malloc(s.size() + 1);
    if (!p) return nullptr;
    std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

void Fr_str2element(PFrElement pE, char const* s, uint base) {
    mp_limb_t mod[4]; Fr_getModulusU256(mod);

    mp_limb_t v[4];
    if (mp_set_str_mod(v, s, (int)base, mod) != 0) {
        mp_set_ui(v, 0);
    }
    Fr_fromU256(pE, v);
}

char *Fr_element2str(PFrElement pE) {
    mp_limb_t v[4];
    Fr_toU256(v, pE);
    return mp_strdup_malloc(mp_get_str(v, 10));
}

void Fr_idiv(PFrElement r, PFrElement a, PFrElement b) {
    mp_limb_t ma[4], mb[4];
    Fr_toU256(ma, a);
    Fr_toU256(mb, b);

    mp_limb_t q[4], rem[4];
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromU256(r, q);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    mp_limb_t ma[4], mb[4];
    Fr_toU256(ma, a);
    Fr_toU256(mb, b);

    mp_limb_t q[4], rem[4];
    if (mp_divmod(q, rem, ma, mb) != 0) {
        throw std::runtime_error("division by zero");
    }
    Fr_fromU256(r, rem);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    mp_limb_t mb[4];
    Fr_toU256(mb, b);

    uint8_t exp_le[32];
    mp_export(exp_le, mb);

    FrRawElement base_norm;
    Fr_toRawNormal(base_norm, a);

    FrRawElement base_mont;
    Fr_rawToMontgomery(base_mont, base_norm);

    FrRawElement res_mont;
    Fr_rawExpMont(res_mont, base_mont, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, res_mont);

    Fr_fromRawNormal(r, res_norm);
}

void Fr_inv(PFrElement r, PFrElement a) {
    uint8_t exp_le[32];
    fr_q_minus_2(exp_le);

    FrRawElement base_norm;
    Fr_toRawNormal(base_norm, a);

    FrRawElement base_mont;
    Fr_rawToMontgomery(base_mont, base_norm);

    FrRawElement res_mont;
    Fr_rawExpMont(res_mont, base_mont, exp_le, (unsigned)sizeof(exp_le));

    FrRawElement res_norm;
    Fr_rawFromMontgomery(res_norm, res_mont);

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
    Fr_init();
    set(fZero, 0);
    set(fOne, 1);
    neg(fNegOne, fOne);
}

RawFr::~RawFr() {}

void RawFr::fromString(Element& r, const std::string& s, uint32_t radix) {
    mp_limb_t mod[4]; Fr_getModulusU256(mod);
    mp_limb_t v[4];
    if (mp_set_str_mod(v, s.c_str(), (int)radix, mod) != 0) {
        mp_set_ui(v, 0);
    }
    r.v[0] = v[0];
    r.v[1] = v[1];
    r.v[2] = v[2];
    r.v[3] = v[3];
    Fr_rawToMontgomery(r.v, r.v);
}

void RawFr::fromUI(Element& r, unsigned long int v) {
    mp_limb_t x[4];
    mp_set_ui(x, (uint64_t)v);
    r.v[0] = x[0];
    r.v[1] = x[1];
    r.v[2] = x[2];
    r.v[3] = x[3];
    Fr_rawToMontgomery(r.v, r.v);
}

RawFr::Element RawFr::set(int value) {
    Element r;
    set(r, value);
    return r;
}

void RawFr::set(Element& r, int value) {
    mp_limb_t mod[4]; Fr_getModulusU256(mod);
    mp_limb_t v[4];
    mp_set_sint_mod(v, (int64_t)value, mod);

    r.v[0] = v[0];
    r.v[1] = v[1];
    r.v[2] = v[2];
    r.v[3] = v[3];
    Fr_rawToMontgomery(r.v, r.v);
}

std::string RawFr::toString(const Element& a, uint32_t radix) {
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, a.v);

    mp_limb_t v[4] = { tmp.v[0], tmp.v[1], tmp.v[2], tmp.v[3] };
    return mp_get_str(v, (int)radix);
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

    mp_limb_t v[4] = { tmp.v[0], tmp.v[1], tmp.v[2], tmp.v[3] };
    mp_export_be(data, v);
    return need;
}

int RawFr::fromRprBE(Element& element, const uint8_t* data, int bytes) {
    const int need = Fr_N64 * 8;
    if (bytes < need) return -need;

    mp_limb_t v[4];
    mp_import_be(v, data);

    element.v[0] = v[0];
    element.v[1] = v[1];
    element.v[2] = v[2];
    element.v[3] = v[3];
    Fr_rawToMontgomery(element.v, element.v);
    return need;
}

static bool init = Fr_init();
RawFr RawFr::field;
