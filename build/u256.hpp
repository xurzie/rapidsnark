#ifndef U256_HPP
#define U256_HPP

#include <cstdint>
#include <cstddef>
#include <string>

using mp_limb_t = uint64_t;
using mp_size_t = std::size_t;

void      mpn_zero(mp_limb_t *rp, mp_size_t n);
void      mpn_copyi(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n);
int       mpn_zero_p(const mp_limb_t *ap, mp_size_t n);
int       mpn_cmp(const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n);
int       mpn_cmp_1(const mp_limb_t *ap, mp_size_t n, mp_limb_t b);
mp_limb_t mpn_add_n(mp_limb_t *rp, const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n);
mp_limb_t mpn_sub_n(mp_limb_t *rp, const mp_limb_t *ap, const mp_limb_t *bp, mp_size_t n);
mp_limb_t mpn_add_1(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, mp_limb_t b);
mp_limb_t mpn_sub_1(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, mp_limb_t b);
int       mpn_tstbit(const mp_limb_t *ap, mp_size_t n, unsigned bit);
void      mpn_rshift(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, unsigned k);
void      mpn_lshift(mp_limb_t *rp, const mp_limb_t *ap, mp_size_t n, unsigned k);

struct U256 {
    mp_limb_t limb[4];
};

static inline mp_limb_t *mp_limbs(U256 *x) { return x->limb; }
static inline const mp_limb_t *mp_limbs(const U256 *x) { return x->limb; }
static inline constexpr mp_size_t mp_nlimbs(const U256 &) { return 4; }

template <std::size_t N>
static inline mp_limb_t *mp_limbs(mp_limb_t (&a)[N]) { return a; }

template <std::size_t N>
static inline const mp_limb_t *mp_limbs(const mp_limb_t (&a)[N]) { return a; }

template <std::size_t N>
static inline constexpr mp_size_t mp_nlimbs(const mp_limb_t (&)[N]) { return N; }

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

void mp_shl_2exp(U256 *r, const U256 *a, uint32_t k);

#endif // U256_HPP
