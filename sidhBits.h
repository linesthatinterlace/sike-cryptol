//Constants.

#define MODULUS 0x0cd07fff
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

typedef uint32_t * mp;

void fp_Add(const mp a, const mp b, mp c);

void fp_Clear(mp* q);

void fp_Constant(unsigned long a, mp b);

void fp_Copy(mp dst, const mp src);

void fp_Init(mp* q);

int fp_IsEqual(const mp a, const mp b);

int fp_IsEven(const mp a);

void fp_Invert(const mp a, mp b);

int fp_IsBitSet(const mp a, const unsigned long i);

int fp_IsConstant(const mp a, const size_t constant);

void fp_Multiply(const mp a, const mp b, mp c);

void fp_Negative(const mp a, mp b);

void fp_Pow(const mp a, const mp b, mp c);

int fp_QuadNonRes(const mp a);

void fp_Square(const mp a, mp b);

void fp_Sqrt(const mp a, mp b);

void fp_Subtract(const mp a, const mp b, mp c);

void fp_Unity(mp b);

void fp_Zero(mp a);

// FP2

typedef struct {
  mp x0;
  mp x1;
} fp2;

void fp2_Add(const fp2* a, const fp2* b, fp2* c );

void fp2_Clear(fp2* fp2);

void fp2_Copy(const fp2* a, fp2* b );

void fp2_Init(fp2* fp2 );

void fp2_Init_set(fp2* fp2, unsigned long x0, unsigned long x1 );

int fp2_IsEqual(const fp2* a1, const fp2* a2 );

void fp2_Set(fp2* fp2, unsigned long x0, unsigned long x1 );

void fp2_Sub(const fp2* a, const fp2* b, fp2* c );

void fp2_Multiply(const fp2*  a, const fp2*  b, fp2*  c );

void fp2_Square(const fp2* a, fp2* b );

void fp2_Invert(const fp2* a, fp2* b );

void fp2_Negative(const fp2* a, fp2* b );

int fp2_IsConst(const fp2* a, unsigned long x0, unsigned long x1 );

void fp2_Sqrt(const fp2* a, fp2* b);


// MONTGOMERY

typedef struct {
  fp2 x;
  fp2 y;
} mont_pt_t;

typedef struct {
  fp2 a;
  fp2 b;
} mont_curve_int_t;

typedef struct {
  fp2 xP;
  fp2 xQ;
  fp2 xR;
} sike_public_key_t;

void mont_curve_init(mont_curve_int_t* curve);

void mont_curve_copy(const mont_curve_int_t* curve, mont_curve_int_t* curvecopy);

void mont_curve_clear(mont_curve_int_t* curve);

void mont_pt_init(mont_pt_t* pt);

void mont_pt_clear(mont_pt_t* pt);

void public_key_init(sike_public_key_t* pk);

void public_key_clear(sike_public_key_t* pk);

void mont_pt_copy(const mont_pt_t* src, mont_pt_t* dst);

void mont_double_and_add(const mont_curve_int_t *curve, const mp k, const mont_pt_t *P, mont_pt_t *Q, int msb);

void xDBL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R);

void xDBLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e, mont_pt_t *R);

void xADD(const mont_curve_int_t *curve, const mont_pt_t *P, const mont_pt_t *Q, mont_pt_t *R);

void xTPL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R);

void xTPLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e, mont_pt_t *R);

void j_inv(const mont_curve_int_t *E, fp2 *jinv);

void get_xR(const mont_curve_int_t *curve, const mont_pt_t *P, const mont_pt_t *Q, sike_public_key_t *pk);

void get_yP_yQ_A_B(const sike_public_key_t *pk, mont_pt_t *P, mont_pt_t *Q, mont_curve_int_t *curve);


// PARAMS

typedef struct {
    mont_curve_int_t startingCurve;
    mont_pt_t param_P2;
    mont_pt_t param_Q2;
    mont_pt_t param_P3;
    mont_pt_t param_Q3;

} sidh_params_t;

typedef struct {
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
} sidh_params_raw_t;

void sidh_setup_params(const sidh_params_raw_t *raw, sidh_params_t *params);
void sidh_teardown_params(sidh_params_t *params);

extern const sidh_params_raw_t sidhRawParams;