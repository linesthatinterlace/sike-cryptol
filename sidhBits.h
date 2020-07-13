#include <stdint.h>
#include <stdlib.h>
// Constants.

#define MODULUS 0x0cd07fff
#define E2 15
#define E3 8
#define MSBA 15
#define MSBB 12
#define XP21 0x0a7e21eb
#define XP22 0x09ab3afb
#define YP21 0x0366529e
#define YP22 0x00615a32
#define XQ21 0x09875b2b
#define XQ22 0x076a3372
#define YQ21 0x00216b6b
#define YQ22 0x0b6465ec
#define XP31 0x09aa9e25
#define XP32 0x00000000
#define YP31 0x040d8353
#define YP32 0x00000000
#define XQ31 0x0b96ade2
#define XQ32 0x00000000
#define YQ31 0x00000000
#define YQ32 0x082a1f97
#define CURVEA1 0x00000006
#define CURVEA2 0x00000000
#define CURVEB1 0x00000001
#define CURVEB2 0x00000000

// FP

typedef uint32_t fp;

void fp_Add(const fp *a, const fp *b, fp *c);

void fp_Constant(uint32_t a, fp *b);

void fp_Copy(const fp *src, fp *dst);

int fp_IsEqual(const fp *a, const fp *b);

int fp_IsEven(const fp *a);

void fp_Invert(const fp *a, fp *b);

int fp_IsBitSet(const fp *a, const int i);

int fp_IsConstant(const fp *a, const size_t constant);

void fp_Multiply(const fp *a, const fp *b, fp *c);

void fp_Negative(const fp *a, fp *b);

void fp_Pow(const fp *a, const fp *b, fp *c);

int fp_QuadNonRes(const fp *a);

void fp_Square(const fp *a, fp *b);

void fp_Sqrt(const fp *a, fp *b);

void fp_Subtract(const fp *a, const fp *b, fp *c);

void fp_Unity(fp *b);

void fp_Zero(fp *a);

// FP2

typedef struct {
  fp x0;
  fp x1;
} fp2;

void fp2_Add(const fp2 *a, const fp2 *b, fp2 *c);

void fp2_Copy(const fp2 *a, fp2 *b);

int fp2_IsEqual(const fp2 *a1, const fp2 *a2);

void fp2_Set(fp2 *fp2, uint32_t x0, uint32_t x1);

void fp2_Sub(const fp2 *a, const fp2 *b, fp2 *c);

void fp2_Multiply(const fp2 *a, const fp2 *b, fp2 *c);

void fp2_Square(const fp2 *a, fp2 *b);

void fp2_Invert(const fp2 *a, fp2 *b);

void fp2_Negative(const fp2 *a, fp2 *b);

int fp2_IsConst(const fp2 *a, uint32_t x0, uint32_t x1);

void fp2_Sqrt(const fp2 *a, fp2 *b);

// MONTGOMERY

typedef struct {
  fp2 x;
  fp2 y;
} mont_pt_t;

typedef struct {
  fp2 a;
  fp2 b;
  mont_pt_t P;
  mont_pt_t Q;
} mont_curve_int_t;

typedef struct {
  fp2 xP;
  fp2 xQ;
  fp2 xR;
} sike_public_key_t;

void mont_curve_copy(const mont_curve_int_t *curve,
                     mont_curve_int_t *curvecopy);

void mont_pt_copy(const mont_pt_t *src, mont_pt_t *dst);

void mont_double_and_add(const mont_curve_int_t *curve, const fp *k,
                         const mont_pt_t *P, mont_pt_t *Q, uint16_t msb);

void xDBL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R);

void xDBLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e,
           mont_pt_t *R);

void xADD(const mont_curve_int_t *curve, const mont_pt_t *P, const mont_pt_t *Q,
          mont_pt_t *R);

void xTPL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R);

void xTPLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e,
           mont_pt_t *R);

void j_inv(const mont_curve_int_t *E, fp2 *jinv);

void get_xR(const mont_curve_int_t *curve, const mont_pt_t *P,
            const mont_pt_t *Q, sike_public_key_t *pk);

void get_yP_yQ_A_B(const sike_public_key_t *pk, mont_pt_t *P, mont_pt_t *Q,
                   mont_curve_int_t *curve);

// PARAMS

typedef struct {
  mont_curve_int_t startingCurve;
  mont_pt_t param_P2;
  mont_pt_t param_Q2;
  mont_pt_t param_P3;
  mont_pt_t param_Q3;
} sike_params_t;

typedef struct {
  uint32_t param_p;
  uint32_t param_cA1;
  uint32_t param_cA2;
  uint32_t param_cB1;
  uint32_t param_cB2;
  uint32_t param_xp21;
  uint32_t param_xp22;
  uint32_t param_yp21;
  uint32_t param_yp22;
  uint32_t param_xq21;
  uint32_t param_xq22;
  uint32_t param_yq21;
  uint32_t param_yq22;
  uint32_t param_xp31;
  uint32_t param_xp32;
  uint32_t param_yp31;
  uint32_t param_yp32;
  uint32_t param_xq31;
  uint32_t param_xq32;
  uint32_t param_yq31;
  uint32_t param_yq32;
} sike_params_raw_t;

void sike_setup_params(const sike_params_raw_t *raw, sike_params_t *params);

extern const sike_params_raw_t sikeRawParams;

// ISOGENIES

void eval_2_iso(const mont_pt_t *P2, const mont_pt_t *P, mont_pt_t *isoP);

void curve_2_iso(const mont_pt_t *P2, const mont_curve_int_t *E,
                 mont_curve_int_t *isoE);

void eval_3_iso(const mont_pt_t *P3, const mont_pt_t *P, mont_pt_t *isoP);

void curve_3_iso(const mont_pt_t *P3, const mont_curve_int_t *E,
                 mont_curve_int_t *isoE);

void eval_4_iso(const mont_pt_t *P4, const mont_pt_t *P, mont_pt_t *isoP);

void curve_4_iso(const mont_pt_t *P4, const mont_curve_int_t *E,
                 mont_curve_int_t *isoE);

void iso_2_e(int e, const mont_curve_int_t *E, mont_pt_t *S,
             const mont_pt_t *P1, const mont_pt_t *P2, mont_curve_int_t *isoE,
             mont_pt_t *isoP1, mont_pt_t *isoP2);

void iso_3_e(int e, const mont_curve_int_t *E, mont_pt_t *S,
             const mont_pt_t *P1, const mont_pt_t *P2, mont_curve_int_t *isoE,
             mont_pt_t *isoP1, mont_pt_t *isoP2);

// SIDH

typedef uint16_t sike_private_key_2;
typedef uint16_t sike_private_key_3;

void sike_isogen_2(const sike_params_t *params, sike_public_key_t *pk,
                   const sike_private_key_2 *sk2);

void sike_isogen_3(const sike_params_t *params, sike_public_key_t *pk,
                   const sike_private_key_2 *sk2);

void sike_isoex_2(const sike_params_t *params, const sike_public_key_t *pkO,
                  const sike_private_key_2 *sk2I, fp2 *secret);

void sike_isoex_3(const sike_params_t *params, const sike_public_key_t *pkO,
                  const sike_private_key_3 *sk3I, fp2 *secret);