#ifndef MP_HPP
#define MP_HPP

#include <cstdint>
#include <cstddef>
#include <string>

#define MP_N64 4
#define MP_N   32

typedef uint64_t mp_uint_t[MP_N64];

void     mp_set(uint64_t *r, uint64_t x);
void     mp_copy(uint64_t *r, const uint64_t *a);

int      mp_cmp(const uint64_t *a, const uint64_t *b);
int      mp_is_zero(const uint64_t *a);

uint64_t mp_add(uint64_t *r, const uint64_t *a, const uint64_t *b);
uint64_t mp_add(uint64_t *rp, const uint64_t *ap, size_t an, const uint64_t *bp, size_t bn);
uint64_t mp_sub(uint64_t *r, const uint64_t *a, const uint64_t *b);

uint64_t mp_add(uint64_t *r, const uint64_t *a, uint64_t b);
uint64_t mp_sub(uint64_t *r, const uint64_t *a, uint64_t b);

uint64_t mp_add_1(uint64_t *r, const uint64_t *a, uint64_t b);
uint64_t mp_sub_1(uint64_t *r, const uint64_t *a, uint64_t b);

uint64_t mp_mul_1(uint64_t *r, const uint64_t *a, uint64_t m);
uint64_t mp_addmul_1(uint64_t *r, const uint64_t *a, uint64_t m);

uint64_t mp_lshift(uint64_t *r, const uint64_t *a, unsigned cnt);
uint64_t mp_rshift(uint64_t *r, const uint64_t *a, unsigned cnt);

void     mp_and(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_ior(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_xor(uint64_t *r, const uint64_t *a, const uint64_t *b);
void     mp_com(uint64_t *r, const uint64_t *a);

int      mp_tstbit(const uint64_t *a, size_t bit);
void     mp_shl(uint64_t *r, const uint64_t *a, uint32_t k);
void     mp_shr(uint64_t *r, const uint64_t *a, uint32_t k);

int         mp_set(uint64_t *r, const char *str, int base);
std::string mp_get_str(const uint64_t *a, int base);

int      mp_fits_int32(const uint64_t *a);

void mp_set_mod(uint64_t *r, int64_t x, const uint64_t *mod);
void mp_add_mod(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *mod);
void mp_mul_mod(uint64_t *r, const uint64_t *a, uint32_t m, const uint64_t *mod);
int  mp_set_mod(uint64_t *r, const char *str, int base, const uint64_t *mod);

void mp_export(uint8_t *out, const uint64_t *a);
void mp_export_be(uint8_t *out, const uint64_t *a);
void mp_import_be(uint64_t *r, const uint8_t *a);

int  mp_divmod(uint64_t *q, uint64_t *r, const uint64_t *num, const uint64_t *den);

#endif
