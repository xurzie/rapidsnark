#ifndef MP_HPP
#define MP_HPP

#include <cstdint>
#include <cstddef>
#include <string>

#define MP_N64 4
#define MP_N   32

typedef uint64_t mp_uint_t[MP_N64];

void     mp_set(uint64_t *r, uint64_t a);
void     mp_copy(uint64_t *r, const uint64_t *a);

int      mp_cmp(const uint64_t *a, const uint64_t *b);
void     mp_zero(uint64_t *r);
bool     mp_is_zero(const uint64_t *a);

uint64_t mp_add(uint64_t *r, const uint64_t *a, const uint64_t *b);
uint64_t mp_add(uint64_t *r, const uint64_t *a, size_t an, const uint64_t *b, size_t bn);
uint64_t mp_sub(uint64_t *r, const uint64_t *a, const uint64_t *b);

uint64_t mp_add(uint64_t *r, const uint64_t *a, uint64_t b);
uint64_t mp_sub(uint64_t *r, const uint64_t *a, uint64_t b);

uint64_t mp_mul(uint64_t *r, const uint64_t *a, uint64_t b);
uint64_t mp_addmul(uint64_t *r, const uint64_t *a, size_t n, uint64_t b);

void     mp_and(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_or(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_xor(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_not(uint64_t *r, const uint64_t *a);

bool     mp_tstbit(const uint64_t *a, size_t bit);
void     mp_shl(uint64_t *r, const uint64_t *a, uint64_t k);
void     mp_shr(uint64_t *r, const uint64_t *a, uint64_t k);

int32_t     mp_get_int32(const mp_uint_t a);
bool        mp_fits_int32(const uint64_t *a);
bool        mp_set(uint64_t *r, const char *str, uint32_t base);
std::string mp_get_str(const uint64_t *a, uint32_t base);

void mp_set_mod(uint64_t *r, int64_t a, const uint64_t *mod);
bool mp_set_mod(uint64_t *r, const char *str, uint32_t base, const uint64_t *mod);

void mp_export_be(uint8_t *r, const uint64_t *a);
void mp_import_be(uint64_t *r, const uint8_t *a);

void mp_div(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den);
void mp_pow_mod(uint64_t *r, const uint64_t *base, const uint64_t *exp, const uint64_t *mod);
bool mp_inv_mod(uint64_t *r, const uint64_t *a, const uint64_t *mod);

#endif
