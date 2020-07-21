#include <stdint.h>
#include <stdlib.h>

// PARAMETER CONSTANTS

#define MODULUS 0x0cd07fff
#define PARE2 15
#define PARE3 8
#define MSBA 15
#define MSBB 12
#define NP   4
#define NSK2 2
#define NSK3 2
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

//
// TYPES
//

typedef uint32_t fp;

typedef struct {
  fp x0;
  fp x1;
} fp2;


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

typedef uint16_t sike_private_key_2;

typedef uint16_t sike_private_key_3;

typedef struct {
  fp cA1;
  fp cA2;
  fp cB1;
  fp cB2;
  fp xp21;
  fp xp22;
  fp yp21;
  fp yp22;
  fp xq21;
  fp xq22;
  fp yq21;
  fp yq22;
  fp xp31;
  fp xp32;
  fp yp31;
  fp yp32;
  fp xq31;
  fp xq32;
  fp yq31;
  fp yq32;
} sike_params_raw_t;

typedef struct {
  mont_curve_int_t EA;
  mont_curve_int_t EB;
} sike_params_t;

//
// PARAMETER SETUP
//

sike_params_t sike_setup_params_C(const sike_params_raw_t raw);

sike_params_raw_t sikeRawParams_C();

//
// ENCODING
//

void itoos_C (const uint32_t * to_enc, uint8_t * enc);

void ostoi_C (const uint8_t* to_dec, uint32_t * dec, size_t len);

void fptoos_C (const fp * to_enc, uint8_t * enc);

int ostofp_C (const uint8_t *to_dec, fp * dec);

void fp2toos_C(const fp2 * to_enc, uint8_t * enc);

int ostofp2_C(const uint8_t *to_dec, fp2 * dec);

void pktoos_C(const sike_public_key_t * to_enc, uint8_t * enc);

int ostopk_C(const uint8_t *to_dec, sike_public_key_t * dec);

//
// FP
//

fp fpAdd_C(const fp a, const fp b);

int fpCmp_C(const fp a, const fp b);

int fpIsBitSet_C(const fp a, const int i);

int fpIsEven_C(const fp a);

fp fpInvert_C(const fp a);

fp fpMultiply_C(const fp a, const fp b);

fp fpNegative_C(const fp a);

fp fpPow_C(const fp a, const fp b);

fp fpSquare_C(const fp a);

fp fpSubtract_C(const fp a, const fp b);

fp fpUnity_C();

fp fpZero_C();

//
// FP2
//

fp2 mkFP2_C(uint32_t x0, uint32_t x1);

fp2 fp2Add_C(const fp2 a, const fp2 b);

fp2 fp2Double_C(const fp2 a);

int fp2Cmp_C(const fp2 a1, const fp2 a2);

fp2 fp2Subtract_C(const fp2 a, const fp2 b);

fp2 fp2Multiply_C(const fp2 a, const fp2 b);

fp2 fp2Square_C(const fp2 a);

fp2 fp2Invert_C(const fp2 a);

fp2 fp2Negative_C(const fp2 a);

fp2 fp2Sqrt_C(const fp2 a);

fp2 fp2Unity_C();

fp2 fp2Zero_C();

//
// MONTGOMERY
//

mont_pt_t mont_double_and_add(const mont_curve_int_t curve, const fp k, const mont_pt_t P, uint16_t msb);
                         
mont_pt_t mont_inf_affine();

int mont_is_inf_affine(const mont_pt_t P);

mont_pt_t xDBL(const mont_curve_int_t curve, const mont_pt_t P);

mont_pt_t xDBLe(const mont_curve_int_t curve, const mont_pt_t P, int e);

mont_pt_t xADD(const mont_curve_int_t curve, const mont_pt_t P, const mont_pt_t Q);

mont_pt_t xTPL(const mont_curve_int_t curve, const mont_pt_t P);

mont_pt_t xTPLe(const mont_curve_int_t curve, const mont_pt_t P, int e);

mont_pt_t xNEGATE(const mont_pt_t P);

fp2 j_inv(const mont_curve_int_t E);

mont_pt_t mkPoint(const fp2 x, const fp2 y); 

mont_curve_int_t mkMC(const fp2 a, const fp2 b, const mont_pt_t P, const mont_pt_t Q); 


// ISOGENIES

mont_curve_int_t curve_2_iso(const mont_pt_t P2, const mont_curve_int_t E);
                 
mont_pt_t eval_2_iso(const mont_pt_t P2, const mont_pt_t P);

mont_curve_int_t curve_4_iso(const mont_pt_t P4, const mont_curve_int_t E);

mont_curve_int_t curve_3_iso(const mont_pt_t P3, const mont_curve_int_t E);

mont_pt_t eval_3_iso(const mont_pt_t P3, const mont_pt_t P);

mont_pt_t eval_4_iso(const mont_pt_t P4, const mont_pt_t P);

mont_curve_int_t iso_2_e(int e, const mont_curve_int_t E, mont_pt_t S);

mont_curve_int_t iso_3_e(int e, const mont_curve_int_t E, mont_pt_t S);

sike_public_key_t get_xR(const mont_curve_int_t curve);

mont_curve_int_t get_yP_yQ_A_B(const sike_public_key_t pk);

// SIDH


sike_public_key_t sike_isogen_2(const sike_params_t params, const sike_private_key_2 sk2);

sike_public_key_t sike_isogen_3(const sike_params_t params, const sike_private_key_3 sk3);

fp2 sike_isoex_2(const sike_public_key_t pkO, const sike_private_key_2 sk2I);

fp2 sike_isoex_3(const sike_public_key_t pkO, const sike_private_key_3 sk3I);

//
// CHECKS
//
fp2 gen3ex2(sike_private_key_2 a, sike_private_key_3 b);

fp2 gen2ex3(sike_private_key_2 a, sike_private_key_3 b);

int secretShared(sike_private_key_2 a, sike_private_key_2 b);

// PRINTING
void printFP2(const fp2 z);
void printPoint(const mont_pt_t P);
void printCurve(const mont_curve_int_t C);
void printPublicKey(const sike_public_key_t K);
void printFP(const fp z);
void printOctetString(const uint8_t* os, size_t len);
