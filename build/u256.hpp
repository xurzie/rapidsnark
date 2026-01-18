#pragma once
#include <cstdint>

struct U256 {
    uint64_t limb[4];
};

void u256_set_ui(U256* r, uint64_t x);
void u256_copy(U256* r, const U256* a);

int  u256_cmp(const U256* a, const U256* b);
int  u256_cmp_ui(const U256* a, uint64_t x);
int  u256_is_zero(const U256* a);

uint64_t u256_add(U256* r, const U256* a, const U256* b);
uint64_t u256_sub(U256* r, const U256* a, const U256* b);
uint64_t u256_add_ui(U256* r, const U256* a, uint64_t b);
uint64_t u256_sub_ui(U256* r, const U256* a, uint64_t b);

int  u256_tstbit(const U256* a, uint32_t bit);
void u256_fdiv_q_2exp(U256* r, const U256* a, uint32_t k);
int  u256_set_str(U256* r, const char* str, int base);
void u256_export(uint8_t out[32], const U256* a);

int  u256_fits_sint(const U256* a);
void u256_set_sint_mod(U256* r, int64_t x, const U256* mod);
void u256_add_mod(U256* r, const U256* a, const U256* b, const U256* mod);
void u256_add_ui_mod(U256* r, const U256* a, uint64_t b, const U256* mod);
void u256_mul_small_mod(U256* r, const U256* a, uint32_t m, const U256* mod);
int  u256_set_str_mod(U256* r, const char* str, int base, const U256* mod);
char* u256_get_str_alloc(const U256* a, int base);

void u256_export_be(uint8_t out[32], const U256* a);
void u256_import_be(U256* r, const uint8_t in[32]);

int  u256_divmod(U256* q, U256* r, const U256* num, const U256* den);

void u256_mont_mul_4(uint64_t out[4], const uint64_t a[4], const uint64_t b[4],
                     const uint64_t mod5[5], uint64_t np);
void u256_mont_mul_4_u64(uint64_t out[4], const uint64_t a[4], uint64_t b,
                         const uint64_t mod5[5], uint64_t np);
void u256_mont_reduce_4(uint64_t out[4], const uint64_t in_mont[4],
                        const uint64_t mod5[5], uint64_t np);
void u256_shl_2exp(U256* r, const U256* a, uint32_t k);

void u256_scalar32_from_u64(uint8_t out[32], uint64_t v);
int  u256_scalar32_from_dec(uint8_t out[32], const char* dec);

