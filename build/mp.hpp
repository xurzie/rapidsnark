#ifndef MP_HPP
#define MP_HPP

#include <string>

using mp_limb_t = uint64_t;

#define MP_N64 4
#define MP_N 32

void mp_set_ui(mp_limb_t *r, uint64_t x);
void mp_copy(mp_limb_t *r, const mp_limb_t *a);

int  mp_cmp(const mp_limb_t *a, const mp_limb_t *b);
int  mp_is_zero(const mp_limb_t *a);

uint64_t mp_add(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
uint64_t mp_sub(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);

uint64_t mp_add(mp_limb_t *r, const mp_limb_t *a, uint64_t b);
uint64_t mp_sub(mp_limb_t *r, const mp_limb_t *a, uint64_t b);

int  mp_tstbit(const mp_limb_t *a, size_t bit);
void mp_fdiv_q_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);
void mp_shl_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);
void mp_shr_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);

int  mp_set_str(mp_limb_t *r, const char *str, int base);
std::string mp_get_str(const mp_limb_t *a, int base);

int  mp_fits_sint(const mp_limb_t *a);
void mp_set_sint_mod(mp_limb_t *r, int64_t x, const mp_limb_t *mod);
void mp_add_mod(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b, const mp_limb_t *mod);
void mp_mul_small_mod(mp_limb_t *r, const mp_limb_t *a, uint32_t m, const mp_limb_t *mod);
int  mp_set_str_mod(mp_limb_t *r, const char *str, int base, const mp_limb_t *mod);

void mp_export(uint8_t *out, const mp_limb_t *a);
void mp_export_be(uint8_t *out, const mp_limb_t *a);
void mp_import_be(mp_limb_t *r, const uint8_t in[32]);

int  mp_divmod(mp_limb_t *q, mp_limb_t *r, const mp_limb_t *num, const mp_limb_t *den);

#endif // MP_HPP