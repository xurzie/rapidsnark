#ifndef U256_HPP
#define U256_HPP

#include <cstdint>
#include <cstddef>
#include <string>

using mp_limb_t = uint64_t;

#ifndef MP_N64
#define MP_N64 4
#endif

using mp_uint_t = mp_limb_t[MP_N64];

struct U256 {
    mp_uint_t limb;
};

// U256 API

void mp_set_ui(U256 *r, uint64_t x);
void mp_copy(U256 *r, const U256 *a);

int  mp_cmp(const U256 *a, const U256 *b);
int  mp_cmp_ui(const U256 *a, uint64_t x);
int  mp_is_zero(const U256 *a);

uint64_t mp_add(U256 *r, const U256 *a, const U256 *b);
uint64_t mp_sub(U256 *r, const U256 *a, const U256 *b);
uint64_t mp_add_ui(U256 *r, const U256 *a, uint64_t b);
uint64_t mp_sub_ui(U256 *r, const U256 *a, uint64_t b);

int  mp_tstbit(const U256 *a, uint32_t bit);
void mp_fdiv_q_2exp(U256 *r, const U256 *a, uint32_t k);
void mp_shl_2exp(U256 *r, const U256 *a, uint32_t k);
void mp_shr_2exp(U256 *r, const U256 *a, uint32_t k);

int  mp_set_str(U256 *r, const char *str, int base);
void mp_export(uint8_t out[32], const U256 *a);

int  mp_fits_sint(const U256 *a);
void mp_set_sint_mod(U256 *r, int64_t x, const U256 *mod);
void mp_add_mod(U256 *r, const U256 *a, const U256 *b, const U256 *mod);
void mp_add_ui_mod(U256 *r, const U256 *a, uint64_t b, const U256 *mod);
void mp_mul_small_mod(U256 *r, const U256 *a, uint32_t m, const U256 *mod);
int  mp_set_str_mod(U256 *r, const char *str, int base, const U256 *mod);
std::string mp_get_str(const U256 *a, int base);

void mp_export_be(uint8_t out[32], const U256 *a);
void mp_import_be(U256 *r, const uint8_t in[32]);

int  mp_divmod(U256 *q, U256 *r, const U256 *num, const U256 *den);

// Limb API

void mp_set_ui(mp_limb_t *r, uint64_t x);
void mp_copy(mp_limb_t *r, const mp_limb_t *a);

int  mp_cmp(const mp_limb_t *a, const mp_limb_t *b);
int  mp_cmp_ui(const mp_limb_t *a, uint64_t x);
int  mp_is_zero(const mp_limb_t *a);

uint64_t mp_add(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
uint64_t mp_sub(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
uint64_t mp_add_ui(mp_limb_t *r, const mp_limb_t *a, uint64_t b);
uint64_t mp_sub_ui(mp_limb_t *r, const mp_limb_t *a, uint64_t b);

int  mp_tstbit(const mp_limb_t *a, size_t bit);
void mp_fdiv_q_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);
void mp_shl_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);
void mp_shr_2exp(mp_limb_t *r, const mp_limb_t *a, uint32_t k);

void mp_export(uint8_t out[32], const mp_limb_t *a);
void mp_export_be(uint8_t out[32], const mp_limb_t *a);
void mp_import_be(mp_limb_t *r, const uint8_t in[32]);

#endif // U256_HPP