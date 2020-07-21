#include "sidh32.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// PARAMS

sike_params_raw_t sikeRawParams_C() {
sike_params_raw_t raw = {
  .cA1 = CURVEA1,
  .cA2 = CURVEA2,
  .cB1 = CURVEB1,
  .cB2 = CURVEB2,
  .xp21 = XP21,
  .xp22 = XP22,
  .yp21 = YP21,
  .yp22 = YP22,
  .xq21 = XQ21,
  .xq22 = XQ22,
  .yq21 = YQ21,
  .yq22 = YQ22,
  .xp31 = XP31,
  .xp32 = XP32,
  .yp31 = YP31,
  .yp32 = YP32,
  .xq31 = XQ31,
  .xq32 = XQ32,
  .yq31 = YQ31,
  .yq32 = YQ32};
  return raw;
}

sike_params_t sike_setup_params_C(const sike_params_raw_t raw) {
  fp2 A, B, xP2, yP2, xQ2, yQ2, xP3, yP3, xQ3, yQ3 = { 0 };
  mont_curve_int_t EA, EB = { 0 };
  mont_pt_t P2, Q2, P3, Q3 = { 0 };
  A = mkFP2_C(raw.cA1, raw.cA2);
  B = mkFP2_C(raw.cB1, raw.cB2);
  xP2 = mkFP2_C(raw.xp21, raw.xp22);
  yP2 = mkFP2_C(raw.yp21, raw.yp22);
  xQ2 = mkFP2_C(raw.xq21, raw.xq22);
  yQ2 = mkFP2_C(raw.yq21, raw.yq22);
  xP3 = mkFP2_C(raw.xp31, raw.xp32);
  yP3 = mkFP2_C(raw.yp31, raw.yp32);
  xQ3 = mkFP2_C(raw.xq31, raw.xq32);
  yQ3 = mkFP2_C(raw.yq31, raw.yq32);
  P2  = mkPoint(xP2, yP2);
  Q2  = mkPoint(xQ2, yQ2);
  P3  = mkPoint(xP3, yP3);
  Q3  = mkPoint(xQ3, yQ3);
  EA  = mkMC(A, B, P2, Q2);
  EB  = mkMC(A, B, P3, Q3);
  sike_params_t params = { .EA = EA, .EB = EB };
  return params;
}


// Encodings
void itoos_C (const uint32_t * to_enc, uint8_t* enc) {
  for (size_t i = 0; i < sizeof(*to_enc); i++) {
    enc[i] = (uint8_t)(*to_enc >> (8 * i));
  }
}

void ostoi_C (const uint8_t* to_dec, uint32_t *dec, size_t len) {
  uint32_t acc = 0;
  for (size_t i = 0; i < len; i++) {
    acc += to_dec[i] << (8 * i);
  }
  *dec = acc;
}

void fptoos_C (const fp * to_enc, uint8_t * enc) {
  fp a = *to_enc % MODULUS;
  itoos_C(&a, enc);
}

int ostofp_C (const uint8_t *to_dec, fp * dec) {
  ostoi_C(to_dec, dec, NP);
  return (*dec < MODULUS);
}

void fp2toos_C(const fp2* to_enc, uint8_t * enc) {
  fptoos_C(&to_enc->x0, enc );
  fptoos_C(&to_enc->x1, enc + NP);  
}


int ostofp_C2(const uint8_t *to_dec, fp2* dec) {
  int rc = 0;
  rc  = ostofp_C(to_dec, &dec->x0);
  rc |= ostofp_C(to_dec + NP, &dec->x1);
  return rc;
}

void pktoos_C(const sike_public_key_t* to_enc, uint8_t * enc) {
  fp2toos_C(&to_enc->xP, enc );
  fp2toos_C(&to_enc->xQ, enc + 2*NP);
  fp2toos_C(&to_enc->xR, enc + 4*NP);
}

int ostopk_C(const uint8_t *to_dec, sike_public_key_t* dec) {
  int rc = 0;
  rc  = ostofp_C2(to_dec, &dec->xP);
  rc |= ostofp_C2(to_dec + 2*NP, &dec->xQ);
  rc |= ostofp_C2(to_dec + 4*NP, &dec->xR);
  return rc;
}

// FP

fp fpAdd_C(const fp a, const fp b) {
  return (uint32_t)(((uint64_t)a + (uint64_t)b) % (uint64_t)MODULUS);
}

int fpCmp_C(const fp a, const fp b) {
  return (a % MODULUS) == (b % MODULUS);
}

int fpIsBitSet_C(const fp a, const int index) {
  return (a & (1 << index)) >> index;
}

int fpIsEven_C(const fp a) { return !fpIsBitSet_C(a, 0); }

fp fpMultiply_C(const fp a, const fp b) {
  return (uint32_t)(((uint64_t)a * (uint64_t)b) % (uint64_t)MODULUS);
}

fp fpNegative_C(const fp a) {
  fp r = a % MODULUS;
  return (r == 0) ? 0 : MODULUS - r;
}

fp fpSquare_C(const fp a) { 
  return fpMultiply_C(a, a);
}

fp fpPow_C(const fp a, const fp b) {
  fp t1;
  t1 = fpUnity_C();
  for (int i = 31; i >= 0; i--) {
    t1 = fpSquare_C(t1);
    if (fpIsBitSet_C(b, i)) {
      t1 = fpMultiply_C(t1, a);
    }
  }
  return t1;
}

fp fpInvert_C(const fp a) {
  fp t1;
  t1 = MODULUS - 2;
  return fpPow_C(a, t1);
}

fp fpSubtract_C(const fp a, const fp b) {
  fp t1;
  t1 = fpNegative_C(b);
  return fpAdd_C(a, t1);
}

fp fpUnity_C() { return 1; }

fp fpZero_C() { return 0; }


// FP2

fp2 mkFP2_C(const fp x0, const fp x1) {
  fp2 z = { .x0 = x0, .x1 = x1 };
  return z;
}

fp2 fp2Add_C(const fp2 a, const fp2 b) {
  fp2 z = { 0 };
  z.x0 = fpAdd_C(a.x0, b.x0);
  z.x1 = fpAdd_C(a.x1, b.x1);
  return z;
}

fp2 fp2Double_C(const fp2 a) {
  return fp2Add_C(a, a);
}

int fp2Cmp_C(const fp2 a1, const fp2 a2) {
  return fpCmp_C(a1.x0, a2.x0) && fpCmp_C(a1.x1, a2.x1);
}

fp2 fp2Subtract_C(const fp2 a, const fp2 b) {
  fp2 z = { 0 };
  z = fp2Add_C(a, fp2Negative_C(b));
  return z;
}

fp2 fp2Multiply_C(const fp2 a, const fp2 b) {
  fp mul0;
  fp mul1;
  fp adda;
  fp addb;
  fp2 z = { 0 };
  mul0 = fpMultiply_C(a.x0, b.x0);
  mul1 = fpMultiply_C(a.x1, b.x1);

  adda = fpAdd_C(a.x0, a.x1);
  addb = fpAdd_C(b.x0, b.x1);

  z.x0 = fpSubtract_C(mul0, mul1);

  mul0 = fpAdd_C(mul0, mul1);
  mul1 = fpMultiply_C(adda, addb);

  z.x1 = fpSubtract_C(mul1, mul0);
  
  return z;
}

fp2 fp2Square_C(const fp2 a) { return fp2Multiply_C(a, a); }

fp2 fp2Invert_C(const fp2 a) {
  fp mul0;
  fp mul1;
  fp2 z = { 0 };
  mul0 = fpSquare_C(a.x0);
  mul1 = fpSquare_C(a.x1);

  mul0 = fpAdd_C(mul0, mul1);
  mul0 = fpInvert_C(mul0);

  mul1 = fpNegative_C(a.x1);

  z.x0 = fpMultiply_C(a.x0, mul0);
  z.x1 = fpMultiply_C(mul1, mul0);
  return z;
}

fp2 fp2Negative_C(const fp2 a) {
  fp2 z = { 0 };
  z.x0 = fpNegative_C(a.x0);
  z.x1 = fpNegative_C(a.x1);
  return z;
}

fp2 fp2Sqrt_C(const fp2 a) {
  fp t0, t1, t2, t3, p14, x0, x1, inv2;
  fp2 z = { 0 };
  
  inv2 = 2;

  // (p + 1) / 4
  p14 = (MODULUS + 1) >> 2;

  t0 = fpSquare_C(a.x0);
  t1 = fpSquare_C(a.x1);
  t1 = fpAdd_C(t0, t1);
  t1 = fpPow_C(t1, p14);
  t0 = fpAdd_C(a.x0, t1);
  if ( fpCmp_C(t0, fpZero_C()) ) {
    t1 = fpAdd_C(t1, t1);
    t0 = fpSubtract_C(t0, t1); 
  }  
  inv2 = fpInvert_C(inv2);
  t0 = fpMultiply_C(t0, inv2); // u
  t1 = fpPow_C(t0, p14);
  t2 = fpInvert_C(t1);
  t2 = fpMultiply_C(t2, a.x1);
  t2 = fpMultiply_C(t2, inv2);
  t3 = fpSquare_C(t1);
  if ( !fpCmp_C(t3, t0) ) {
    t3 = t2;
    t2 = t1;
    t1 = t3;
  }
  if ( (!fpIsEven_C(t1)) || ( (fpCmp_C(t1, fpZero_C())) && (!fpIsEven_C(t2)) ) ) {
    x0 = fpNegative_C(t1);
    x1 = fpNegative_C(t2);
  } else {
    x0 = t1;
    x1 = t2;
  }
  z.x0 = x0;
  z.x1 = x1;
  return z;
}

fp2 fp2Unity_C() { return mkFP2_C(1, 0); }

fp2 fp2Zero_C() { return mkFP2_C(0, 0); }

//
// MONTGOMERY
//

mont_pt_t mkPoint(const fp2 x, const fp2 y) {
  mont_pt_t point = {.x = x, .y = y};
  return point;
} 

mont_curve_int_t mkMC(const fp2 a, const fp2 b, const mont_pt_t P, const mont_pt_t Q) {
  mont_curve_int_t curve = {.a = a, .b = b, .P = P, .Q = Q};
  return curve;
} 

// infinity is represented as a point with (0, 0) 
mont_pt_t mont_inf_affine() {
  mont_pt_t P = { 0 };
  P.x = fp2Zero_C();
  P.y = fp2Unity_C();
  return P;
}

// returns 1 for True, 0 for False 
int mont_is_inf_affine(const mont_pt_t P) {
  return (fp2Cmp_C(P.x, mont_inf_affine().x)) && (fp2Cmp_C(P.y, mont_inf_affine().y));
}

mont_pt_t mont_double_and_add(const mont_curve_int_t curve, const fp k, const mont_pt_t P, uint16_t msb) {

  int i;

  mont_pt_t kP = { 0 };

  kP = mont_inf_affine();
  for (i = msb - 1; i >= 0; i--) {
    kP = xDBL(curve, kP);
    if (fpIsBitSet_C(k, i)) {
      kP = xADD(curve, kP, P);
    }
  }

  return kP;
}

mont_pt_t xDBL(const mont_curve_int_t curve, const mont_pt_t P) {

  const fp2 a = curve.a;
  const fp2 b = curve.b;

  mont_pt_t R = { 0 };

  // x3 = b*(3*x1^2+2*a*x1+1)^2/(2*b*y1)^2-a-x1-x1
  // y3 =
  // (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)-b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3-y1

  fp2 t0 = { 0 }, t1 = { 0 }, t2 = { 0 }, xR = { 0 }, yR = { 0 };

  t0 = fp2Negative_C(P.y);

  if (mont_is_inf_affine(P)) {
    R = mont_inf_affine();
  } else if (fp2Cmp_C(P.y, t0)) {
    // P == -P 
    R = mont_inf_affine();
  } else {

    t2 = fp2Unity_C(); // t2 = 1

    t0 = fp2Square_C(P.x); // t0 = x1^2
    t1 = fp2Double_C(t0); // t1 = 2*x1^2
    t0 = fp2Add_C(t0, t1); // t0 = 3*x1^2

    t1 = fp2Multiply_C(a, P.x); // t1 = a*x1
    t1 = fp2Double_C(t1);      // t1 = 2*a*x1

    t0 = fp2Add_C(t0, t1); // t0 = 3*x1^2+2*a*x1
    t0 = fp2Add_C(t0, t2); // t0 = 3*x1^2+2*a*x1+1

    t1 = fp2Multiply_C(b, P.y); // t1 = b*y1
    t1 = fp2Double_C(t1);      // t1 = 2*b*y1
    t1 = fp2Invert_C(t1);        // t1 = 1 / (2*b*y1)

    t0 = fp2Multiply_C(t0, t1); // t0 = (3*x1^2+2*a*x1+1) / (2*b*y1)

    t1 = fp2Square_C(t0); // t1 = (3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2

    t2 = fp2Multiply_C(b, t1); // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2
    t2 = fp2Subtract_C(t2, a);      // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a
    t2 = fp2Subtract_C(t2, P.x);
    
    // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1
    t2 = fp2Subtract_C(t2, P.x);
    
    // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1 - x1

    t1 = fp2Multiply_C(t0, t1); // t1 = (3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3
    t1 = fp2Multiply_C(b, t1);   // t1 = b*(3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3
    t1 = fp2Add_C(t1, P.y);    // t1 = b*(3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3 + y1

    yR = fp2Double_C(P.x); // x3 = 2*x1
    yR = fp2Add_C(yR, P.x); // y3 = 2*x1+x1
    yR = fp2Add_C(yR, a);     // y3 = 2*x1+x1+a
    yR = fp2Multiply_C(yR, t0);        // y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)
    yR = fp2Subtract_C(yR, t1); // y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1) - (b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3 + y1)
    xR = t2; // x3 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1 - x1
    
    R = mkPoint(xR, yR);
  }
  return R;
}

mont_pt_t xDBLe(const mont_curve_int_t curve, const mont_pt_t P, int e) {
  mont_pt_t R = P;
  for (int j = 0; j < e; ++j) {
    R = xDBL(curve, R);
  }
  return R;
}

mont_pt_t xADD(const mont_curve_int_t curve, const mont_pt_t P, const mont_pt_t Q) {

  // x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2
  // y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1
  // y3 = ((2*x1)+x2+a) * ((y2-y1)/(x2-x1)) - b*((y2-y1)^3/(x2-x1)^3) - y1

  const fp2 a = curve.a;
  const fp2 b = curve.b;

  mont_pt_t R = { 0 };

  fp2 t0 = { 0 };

  t0 = fp2Negative_C(Q.y);

  if (mont_is_inf_affine(P)) {
    R = Q;
  } else if (mont_is_inf_affine(Q)) {
    R = P;
  } else if (fp2Cmp_C(P.x, Q.x) && fp2Cmp_C(P.y, Q.y)) {
    // P == Q 
    R = xDBL(curve, P);
  } else if (fp2Cmp_C(P.x, Q.x) && fp2Cmp_C(P.y, t0)) {
    // P == -Q 
    R = mont_inf_affine();
  } else {
    // P != Q or -Q  
    fp2 t1 = { 0 }, t2 = { 0 };
    t0 = fp2Subtract_C(Q.y, P.y);  // t0 = y2-y1
    t1 = fp2Subtract_C(Q.x, P.x);  // t1 = x2-x1
    t1 = fp2Invert_C(t1);        // t1 = 1/(x2-x1)
    t0 = fp2Multiply_C(t0, t1); // t0 = (y2-y1)/(x2-x1)

    t1 = fp2Square_C(t0); // t1 = (y2-y1)^2/(x2-x1)^2

    t2 = fp2Double_C(P.x);  // t2 = 2*x1
    t2 = fp2Add_C(t2, Q.x);    // t2 = 2*x1+x2
    t2 = fp2Add_C(t2, a);        // t2 = 2*x1+x2+a
    t2 = fp2Multiply_C(t2, t0); // t2 = (2*x1+x2+a)*(y2-y1)/(x2-x1)

    t0 = fp2Multiply_C(t0, t1); // t0 = (y2-y1)^3/(x2-x1)^3
    t0 = fp2Multiply_C(b, t0);   // t0 = b*(y2-y1)^3/(x2-x1)^3
    t0 = fp2Add_C(t0, P.y);    // t0 = b*(y2-y1)^3/(x2-x1)^3+y1

    t0 = fp2Subtract_C(t2, t0); // t0 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1

    t1 = fp2Multiply_C(b, t1); // t1 = b*(y2-y1)^2/(x2-x1)^2
    t1 = fp2Subtract_C(t1, a);      // t1 = b*(y2-y1)^2/(x2-x1)^2-a
    t1 = fp2Subtract_C(t1, P.x);  // t1 = b*(y2-y1)^2/(x2-x1)^2-a-x1

    R.x = fp2Subtract_C(t1, Q.x); // x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2

    R.y = t0; // y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-(b*(y2-y1)^3/(x2-x1)^3+y1)
  }
  return R;
}

mont_pt_t xTPL(const mont_curve_int_t curve, const mont_pt_t P) {
  mont_pt_t T = { 0 };

  T = xDBL(curve, P);
  T = xADD(curve, P, T);
  return T;
}

mont_pt_t xTPLe(const mont_curve_int_t curve, const mont_pt_t P, int e) {
  mont_pt_t R = P;
  for (int j = 0; j < e; ++j) {
    R = xTPL(curve, R);
  }
  return R;
}

mont_pt_t xNEGATE(const mont_pt_t P) {
  mont_pt_t R = P;
  R.y = fp2Negative_C(R.y);
  return R;
}

fp2 j_inv(const mont_curve_int_t E) {
  const fp2 a = E.a;

  fp2 t0 = { 0 }, t1 = { 0 }, jinv = { 0 };

  t0 = fp2Square_C(a);            // t0 = a^2
  jinv = mkFP2_C(3, 0);           // jinv = 3
  jinv = fp2Subtract_C(t0, jinv);      // jinv = a^2-3
  t1 = fp2Square_C(jinv);         // t1 = (a^2-3)^2
  jinv = fp2Multiply_C(jinv, t1); // jinv = (a^2-3)^3
  jinv = fp2Double_C(jinv);     // jinv = 2*(a^2-3)^3
  jinv = fp2Double_C(jinv);     // jinv = 4*(a^2-3)^3
  jinv = fp2Double_C(jinv);     // jinv = 8*(a^2-3)^3
  jinv = fp2Double_C(jinv);     // jinv = 16*(a^2-3)^3
  jinv = fp2Double_C(jinv);     // jinv = 32*(a^2-3)^3
  jinv = fp2Double_C(jinv);      // jinv = 64*(a^2-3)^3
  jinv = fp2Double_C(jinv);    // jinv = 128*(a^2-3)^3
  jinv = fp2Double_C(jinv);   // jinv = 256*(a^2-3)^3

  t1 = mkFP2_C(4, 0);            // t1 = 4
  t0 = fp2Subtract_C(t0, t1);        // t0 = a^2-4
  t0 = fp2Invert_C(t0);          // t0 = 1/(a^2-4)
  jinv = fp2Multiply_C(jinv, t0); // jinv = 256*(a^2-3)^3/(a^2-4)
  
  return jinv;
}


// ISOGENIES

mont_pt_t eval_2_iso(const mont_pt_t P2, const mont_pt_t P) {

  const fp2 x2 = P2.x;
  const fp2 x = P.x;
  const fp2 y = P.y;
  mont_pt_t isoP = { 0 };

  // xx:=(x^2*x2-x)/(x-x2);
  // yy:=y*(x^2*x2-2*x*x2^2+x2)/(x-x2)^2;

  fp2 t1 = {0}, t2 = {0}, t3 = {0};

  t1 = fp2Multiply_C(x, x2);   // t1:=x*x2;  // t1 = x*x2
  t2 = fp2Multiply_C(x, t1);  // t2:=x*t1;  // t2 = x^2*x2
  t3 = fp2Multiply_C(t1, x2); // t3:=t1*x2; // t3 = x*x2^2
  t3 = fp2Double_C(t3);     // t3:=t3+t3; // t3 = 2*x*x2^2
  t3 = fp2Subtract_C(t2, t3);     // t3:=t2-t3; // t3 = x^2*x2-2*x*x2^2
  t3 = fp2Add_C(t3, x2);      // t3:=t3+x2; // t3 = x^2*x2-2*x*x2^2+x2
  t3 = fp2Multiply_C(y, t3);  // t3:=y*t3;  // t3 = y*(x^2*x2-2*x*x2^2+x2)
  t2 = fp2Subtract_C(t2, x);       // t2:=t2-x;  // t2 = x^2*x2-x
  t1 = fp2Subtract_C(x, x2);        // t1:=x-x2;  // t1 = x-x2
  t1 = fp2Invert_C(t1);       // t1:=1/t1;  // t1 = 1/(x-x2)
  isoP.x = fp2Multiply_C(t2, t1); // xx:=t2*t1; // xx = (x^2*x2-x)/(x-x2)
  t1 = fp2Square_C(t1);       // t1:=t1^2;  // t1 = 1/(x-x2)^2
  isoP.y = fp2Multiply_C(t3, t1); // yy:=t3*t1; // yy = y*(x^2*x2-2*x*x2^2+x2)/(x-x2)^2
  
  return isoP;
}

mont_curve_int_t curve_2_iso(const mont_pt_t P2, const mont_curve_int_t E) {

  // aa:=2*(1-2*x2^2);
  // bb:=x2*b;

  const fp2 x2 = P2.x;
  const fp2 b = E.b;

  mont_curve_int_t isoE = { 0 };
  isoE.P = E.P;
  isoE.Q = E.Q;

  fp2 t1 = { 0 }, one = {0};
  one = fp2Unity_C();

  t1 = fp2Square_C(x2);     // t1 = x2^2
  t1 = fp2Double_C(t1);  // t1 = 2*x2^2
  t1 = fp2Subtract_C(one, t1); // t1 = 1-2*x2^2
  isoE.a = fp2Double_C(t1);   // aa = 2*(1-2*x2^2)
  isoE.b = fp2Multiply_C(x2, b); // bb = x2*b
  
  return isoE;
}


mont_pt_t eval_3_iso(const mont_pt_t P3, const mont_pt_t P) {

  const fp2 x3 = P3.x;
  const fp2 x = P.x;
  const fp2 y = P.y;
  mont_pt_t isoP = { 0 };
  
  // xx:=x*(x*x3-1)^2/(x-x3)^2;
  // yy:=y*(x*x3-1)*(x^2*x3-3*x*x3^2+x+x3)/(x-x3)^3;

  fp2 t1 = {0}, t2 = {0}, t3 = {0}, t4 = {0}, one = {0};

  one = fp2Unity_C();

  t1 = fp2Square_C(x);         // t1 = x^2
  t1 = fp2Multiply_C(t1, x3); // t1 = x^2*x3
  t2 = fp2Square_C(x3);        // t2 = x3^2
  t2 = fp2Multiply_C(x, t2);  // t2 = x*x3^2
  t3 = fp2Double_C(t2);     // t3 = 2*x*x3^2
  t2 = fp2Add_C(t2, t3);     // t2 = 3*x*x3^2
  t1 = fp2Subtract_C(t1, t2);     // t1 = x^2*x3-3*x*x3^2
  t1 = fp2Add_C(t1, x);       // t1 = x^2*x3-3*x*x3^2+x
  t1 = fp2Add_C(t1, x3);      // t1 = x^2*x3-3*x*x3^2+x+x3

  t2 = fp2Subtract_C(x, x3);         // t2 = x-x3
  t2 = fp2Invert_C(t2);        // t2 = 1/(x-x3)
  t3 = fp2Square_C(t2);        // t3 = 1/(x-x3)^2
  t2 = fp2Multiply_C(t2, t3); // t2 = 1/(x-x3)^3

  t4 = fp2Multiply_C(x, x3); // t4 = x*x3
  t4 = fp2Subtract_C(t4, one);  // t4 = x*x3-1

  t1 = fp2Multiply_C(t4, t1); // t1 = (x*x3-1)*(x^2*x3-3*x*x3^2+x+x3)
  t1 = fp2Multiply_C(t1, t2); // t1 = (x*x3-1)*(x^2*x3-3*x*x3^2+x+x3)/(x-x3)^3

  t2 = fp2Square_C(t4);        // t2 = (x*x3-1)^2
  t2 = fp2Multiply_C(t2, t3); // t2 = (x*x3-1)^2/(x-x3)^2

  isoP.x = fp2Multiply_C(x, t2); // xx = x*(x*x3-1)^2/(x-x3)^2
  isoP.y = fp2Multiply_C(y, t1); // yy = y*(x*x3-1)*(x^2*x3-3*x*x3^2+x+x3)/(x-x3)^3
  
  return isoP;
}


mont_curve_int_t curve_3_iso(const mont_pt_t P3, const mont_curve_int_t E) {

  // aa:=(a*x3-6*x3^2+6)*x3;
  // bb:=b*x3^2;

  const fp2 x3 = P3.x;
  const fp2 a = E.a;
  const fp2 b = E.b;
  
  mont_curve_int_t isoE = { 0 };
  isoE.P = E.P;
  isoE.Q = E.Q;
  
  fp2 t1 = {0}, t2 = {0};

  t1 = fp2Square_C(x3);      // t1 = x3^2
  isoE.b = fp2Multiply_C(b, t1); // bb = b*x3^2

  t1 = fp2Double_C(t1);    // t1 = 2*x3^2
  t2 = fp2Double_C(t1);    // t2 = 4*x3^2
  t1 = fp2Add_C(t1, t2);    // t1 = 6*x3^2
  t2 = mkFP2_C(6, 0);        // t2 = 6
  t1 = fp2Subtract_C(t1, t2);    // t1 = 6*x3^2-6
  t2 = fp2Multiply_C(a, x3);  // t2 = a*x3
  t1 = fp2Subtract_C(t2, t1);    // t1 = a*x3-6*x3^2+6
  isoE.a = fp2Multiply_C(t1, x3); // aa = (a*x3-6*x3^2+6)*x3
  
  return isoE;
}


mont_pt_t eval_4_iso(const mont_pt_t P4, const mont_pt_t P) {

  const fp2 x4 = P4.x;
  const fp2 x = P.x;
  const fp2 y = P.y;
  mont_pt_t isoP = { 0 };

  fp2 t1 = {0}, t2 = {0}, t3 = {0}, t4 = {0}, t5 = {0}, one = {0};
  one = fp2Unity_C();

  // xx:=-(x*x4^2+x-2*x4)*x*(x*x4-1)^2/((x-x4)^2*(2*x*x4-x4^2-1));
  // yy:=-2*y*x4^2*(x*x4-1)*(x^4*x4^2-4*x^3*x4^3+2*x^2*x4^4+x^4-4*x^3*x4+10*x^2*x4^2-4*x*x4^3-4*x*x4+x4^2+1)/((x-x4)^3*(2*x*x4-x4^2-1)^2);

  t1 = fp2Square_C(x);          // t1 = x^2
  t2 = fp2Square_C(t1);        // t2 = x^4
  t3 = fp2Square_C(x4);         // t3 = x4^2
  t4 = fp2Multiply_C(t2, t3); // t4 = x^4*x4^2
  t2 = fp2Add_C(t2, t4);      // t2 = x^4+x^4*x4^2
  t4 = fp2Multiply_C(t1, t3); // t4 = x^2*x4^2
  t4 = fp2Double_C(t4);      // t4 = 2*x^2*x4^2
  t5 = fp2Double_C(t4);      // t5 = 4*x^2*x4^2
  t5 = fp2Double_C(t5);      // t5 = 8*x^2*x4^2
  t4 = fp2Add_C(t4, t5);      // t4 = 10*x^2*x4^2
  t2 = fp2Add_C(t2, t4);      // t2 = x^4+x^4*x4^2+10*x^2*x4^2
  t4 = fp2Multiply_C(t3, t3); // t4 = x4^4
  t5 = fp2Multiply_C(t1, t4); // t5 = x^2*x4^4
  t5 = fp2Double_C(t5);      // t5 = 2*x^2*x4^4
  t2 = fp2Add_C(t2, t5);      // t2 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4
  t1 = fp2Multiply_C(t1, x);   // t1 = x^3
  t4 = fp2Multiply_C(x4, t3);  // t4 = x4^3
  t5 = fp2Multiply_C(t1, t4); // t5 = x^3*x4^3
  t5 = fp2Double_C(t5);      // t5 = 2*x^3*x4^3
  t5 = fp2Double_C(t5);      // t5 = 4*x^3*x4^3
  t2 = fp2Subtract_C(t2, t5); // t2 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3
  t1 = fp2Multiply_C(t1, x4); // t1 = x^3*x4
  t1 = fp2Double_C(t1);     // t1 = 2*x^3*x4
  t1 = fp2Double_C(t1);     // t1 = 4*x^3*x4
  t1 = fp2Subtract_C(t2, t1); // t1 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4
  t2 = fp2Multiply_C(x, t4); // t2 = x*x4^3
  t2 = fp2Double_C(t2);    // t2 = 2*x*x4^3
  t2 = fp2Double_C(t2);    // t2 = 4*x*x4^3
  t1 = fp2Subtract_C( t1, t2); // t1 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3
  t1 = fp2Add_C(t1, t3);  // t1 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2
  t1 = fp2Add_C(t1, one); // t1 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1
  t2 = fp2Multiply_C(x, x4); // t2 = x*x4
  t4 = fp2Subtract_C(t2, one);  // t4 = x*x4-1
  t2 = fp2Double_C(t2);   // t2 = 2*x*x4
  t5 = fp2Double_C(t2);   // t5 = 4*x*x4
  t1 = fp2Subtract_C(t1, t5); // t1 = x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4
  t1 = fp2Multiply_C(t4, t1); // t1 = (x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)
  t1 = fp2Multiply_C( t3, t1); // t1 = x4^2*(x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)
  t1 = fp2Multiply_C(y, t1); // t1 =x4^2*(x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)
  t1 = fp2Double_C(t1); // t1 = // 2*x4^2*(x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)
  isoP.y = fp2Negative_C(t1); // yy = -2*x4^2*(x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)
  t2 = fp2Subtract_C(t2, t3);      // t2 = 2*x*x4-x4^2
  t1 = fp2Subtract_C(t2, one);     // t1 = 2*x*x4-x4^2-1
  t2 = fp2Subtract_C(x, x4);         // t2 = x-x4
  t1 = fp2Multiply_C(t2, t1); // t1 = (x-x4)*(2*x*x4-x4^2-1)
  t5 = fp2Square_C(t1);        // t5 = (x-x4)^2*(2*x*x4-x4^2-1)^2
  t5 = fp2Multiply_C(t5, t2); // t5 = (x-x4)^3*(2*x*x4-x4^2-1)^2
  t5 = fp2Invert_C(t5);        // t5 = 1/((x-x4)^3*(2*x*x4-x4^2-1)^2)
  isoP.y = fp2Multiply_C(isoP.y, t5);
    // yy = -2*x4^2*(x*x4-1)*(x^4+x^4*x4^2+10*x^2*x4^2+2*x^2*x4^4-4*x^3*x4^3-4*x^3*x4-4*x*x4^3+x4^2+1-4*x*x4)/((x-x4)^3*(2*x*x4-x4^2-1)^2)

  t1 = fp2Multiply_C(t1, t2); // t1 = (x-x4)^2*(2*x*x4-x4^2-1)
  t1 = fp2Invert_C(t1);        // t1 = 1/((x-x4)^2*(2*x*x4-x4^2-1))
  t4 = fp2Square_C(t4);        // t4 = (x*x4-1)^2
  t1 = fp2Multiply_C(t1, t4); // t1 = (x*x4-1)^2/((x-x4)^2*(2*x*x4-x4^2-1))
  t1 = fp2Multiply_C(x, t1);   // t1 = x*(x*x4-1)^2/((x-x4)^2*(2*x*x4-x4^2-1))
  t2 = fp2Multiply_C(x, t3);   // t2 = x*x4^2
  t2 = fp2Add_C(t2, x);        // t2 = x*x4^2+x
  t3 = fp2Double_C(x4);        // t3 = 2*x4
  t2 = fp2Subtract_C(t2, t3);      // t2 = x*x4^2+x-2*x4
  t2 = fp2Negative_C(t2);      // t2 = -(x*x4^2+x-2*x4)
  isoP.x = fp2Multiply_C(t1, t2);
    // xx = -(x*x4^2+x-2*x4)*x*(x*x4-1)^2/((x-x4)^2*(2*x*x4-x4^2-1))
    
  return isoP;
}

mont_curve_int_t curve_4_iso(const mont_pt_t P4, const mont_curve_int_t E) {

  const fp2 x4 = P4.x;
  const fp2 b = E.b;
  
  mont_curve_int_t isoE = { 0 };
  isoE.P = E.P;
  isoE.Q = E.Q;

  fp2 t1 = {0}, t2 = {0};

  // aa:=4*x4^4-2;
  // bb:=-(1/2)*x4*(x4^2+1)*b = -(1/2)*(x4^3+x4)*b;

  t1 = fp2Square_C(x4);  // t1 = x4^2
  isoE.a = fp2Square_C(t1);  // aa = x4^4
  isoE.a = fp2Double_C(isoE.a);  // aa = 2*x4^4
  isoE.a = fp2Double_C(isoE.a);  // aa = 4*x4^4
  t2 = mkFP2_C(2, 0);   // t2 = 2
  isoE.a = fp2Subtract_C(isoE.a, t2); // aa = 4*x4^4-2

  t1 = fp2Multiply_C(x4, t1); // t1 = x4^3
  t1 = fp2Add_C(t1, x4);      // t1 = x4^3+x4
  t1 = fp2Multiply_C(t1, b);  // t1 = (x4^3+x4)*b
  t2 = fp2Invert_C(t2);       // t2 = 1/2 -> precofpute
  t2 = fp2Negative_C(t2);     // t2 = -(1/2)
  isoE.b = fp2Multiply_C(t2, t1);
  
  return isoE;
}


mont_curve_int_t iso_2_e(int e, const mont_curve_int_t E, mont_pt_t S) {
  
  const mont_pt_t P1 = E.P;
  const mont_pt_t P2 = E.Q;
  
  int p1Inf = mont_is_inf_affine(P1), p2Inf = mont_is_inf_affine(P2);

  mont_pt_t T = {0};
  mont_curve_int_t isoE = E;

  if (e % 2) {
    T = xDBLe(isoE, S, e - 1);
    isoE = curve_2_iso(T, isoE);
    S = eval_2_iso(T, S);
    if (!p1Inf) isoE.P = eval_2_iso(T, isoE.P);
    if (!p2Inf) isoE.Q = eval_2_iso(T, isoE.Q);
    e--;
  }
  
  for (int i = e - 2; i >= 0; i -= 2) {
    T = xDBLe(isoE, S, i);
    isoE = curve_4_iso(T, isoE);

    S = eval_4_iso(T, S);
    
    if (!p1Inf) isoE.P = eval_4_iso(T, isoE.P);
    if (!p2Inf)isoE.Q = eval_4_iso(T, isoE.Q);

  }
  return isoE;
}


mont_curve_int_t iso_3_e(int e, const mont_curve_int_t E, mont_pt_t S) {

  const mont_pt_t P1 = E.P;
  const mont_pt_t P2 = E.Q;

  int p1Inf = mont_is_inf_affine(P1), p2Inf = mont_is_inf_affine(P2);

  mont_pt_t T = {0};
  mont_curve_int_t isoE = E;


  for (int i = e - 1; i >= 0; --i) {
    T = xTPLe(isoE, S, i);
    isoE = curve_3_iso(T, isoE);

    S = eval_3_iso(T, S);

    if (!p1Inf) isoE.P = eval_3_iso(T, isoE.P);
    if (!p2Inf) isoE.Q = eval_3_iso(T, isoE.Q);
    
  }
  
  return isoE;
}


sike_public_key_t mkPublicKey(const fp2 xP, const fp2 xQ, const fp2 xR) {
  sike_public_key_t pk = {.xP = xP, .xQ = xQ, .xR = xR};
  return pk;
} 

sike_public_key_t get_xR(const mont_curve_int_t curve) {

  mont_pt_t R = {0};
  sike_public_key_t pk = { 0 };
  const fp2 xP = curve.P.x;
  const fp2 xQ = curve.Q.x;
  
  
  R = xNEGATE(curve.Q);
  R = xADD(curve, curve.P, R);
  pk = mkPublicKey(xP, xQ, R.x);

  return pk;
}


fp2 get_A(const sike_public_key_t pk) {
  fp2 t0 = {0}, t1 = {0}, one = {0};
  one = fp2Unity_C();
  fp2 a = {0};
  
  const fp2 xP = pk.xP;
  const fp2 xQ = pk.xQ;
  const fp2 xR = pk.xR;
  
  t1 = fp2Add_C(xP, xQ);
  t0 = fp2Multiply_C(xP, xQ);
  a = fp2Multiply_C(xR, t1);
  a = fp2Add_C(t0, a);
  t0 = fp2Multiply_C(xR, t0);
  a = fp2Subtract_C(a, one);
  t0 = fp2Double_C(t0);
  t1 = fp2Add_C(t1, xR);
  t0 = fp2Double_C(t0);
  a = fp2Square_C(a);
  t0 = fp2Invert_C(t0);
  a = fp2Multiply_C(a, t0);
  a = fp2Subtract_C(a, t1);
  
  return a;
}


mont_curve_int_t get_yP_yQ_A_B(const sike_public_key_t pk) {

  mont_curve_int_t curve = {0};
  const fp2 xP = pk.xP;
  const fp2 xQ = pk.xQ;
  const fp2 xR = pk.xR;

  mont_pt_t T = {0};
  mont_pt_t P = {0};
  mont_pt_t Q = {0};

  fp2 t1 = { 0 }, t2 = { 0 }, a = { 0 }, b = { 0 };
  
  a = get_A(pk);
  
  t1 = fp2Square_C(xP);       // t1 = xP^2
  t2 = fp2Multiply_C(xP, t1); // t2 = xP^3
  t1 = fp2Multiply_C(a, t1);  // t1 = a*xP^2
  t1 = fp2Add_C(t2, t1);      // t1 = xP^3+a*xP^2
  t1 = fp2Add_C(t1, xP);      // t1 = xP^3+a*xP^2+xP
  P.y = fp2Sqrt_C(t1);      // yP = sqrt(xP^3+a*xP^2+xP)
  t1 = fp2Square_C(xQ);       // t1 = xQ^2
  t2 = fp2Multiply_C(xQ, t1); // t2 = xQ^3
  t1 = fp2Multiply_C(a, t1);  // t1 = a*xQ^2
  t1 = fp2Add_C(t2, t1);      // t1 = xQ^3+a*xQ^2
  t1 = fp2Add_C(t1, xQ);      // t1 = xQ^3+a*xQ^2+xQ
  Q.y = fp2Sqrt_C(t1);      // yQ = sqrt(xQ^3+a*xQ^2+xQ)

  P.x = xP;
  Q.x = xQ;
  
  b = fp2Unity_C();
  
  curve = mkMC(a, b, P, Q);

  T = xNEGATE(Q);
  T = xADD(curve, P, T);

  if (!fp2Cmp_C(T.x, xR)) {    
    curve.Q.y = fp2Negative_C(curve.Q.y);
  }
  
  return curve;
}


//
// SIDH
//

sike_public_key_t sike_isogen_2(const sike_params_t params, const sike_private_key_2 sk2) {

  sike_public_key_t pk = { 0 };
                     
  const mont_curve_int_t C = params.EA;
  const mont_curve_int_t D = params.EB;

  mont_curve_int_t pkInt = D;

  uint16_t e = PARE2;
  uint16_t msb = MSBA;

  mont_pt_t S = {0};

  fp sk = sk2;
  
  // Generate kernel
  // S:=P2+SK_2*Q2;
  S = mont_double_and_add(C, sk, C.Q, msb);
  S = xADD(C, C.P, S);
  pkInt = iso_2_e((int)e, pkInt, S);
  pk = get_xR(pkInt);
  
  return pk;
}


sike_public_key_t sike_isogen_3(const sike_params_t params, const sike_private_key_3 sk3) {
  
  sike_public_key_t pk = { 0 };
                     
  const mont_curve_int_t C = params.EB;
  const mont_curve_int_t D = params.EA;

  mont_curve_int_t pkInt = D;

  uint16_t e = PARE3;
  uint16_t msb = MSBB;

  mont_pt_t S = {0};

  fp sk = sk3;
  
  // Generate kernel
  // S:=P2+SK_2*Q2;
  S = mont_double_and_add(C, sk, C.Q, msb);
  S = xADD(C, C.P, S);
  pkInt = iso_3_e((int)e, pkInt, S);       
  pk = get_xR(pkInt);
  
  return pk;
}

fp2 sike_isoex_2(const sike_public_key_t pkO, const sike_private_key_2 sk2I) {
  
  fp2 secret = { 0 };
  mont_curve_int_t E = {0};
  
  uint16_t e = PARE2;
  uint16_t msb = MSBA;

  mont_pt_t S = {0};
  E = get_yP_yQ_A_B(pkO);

  fp skI = sk2I;

  S = mont_double_and_add(E, skI, E.Q, msb);
  S = xADD(E, E.P, S);
  E.P = mont_inf_affine();
  E.Q = mont_inf_affine();
  E = iso_2_e((int)e, E, S);
  secret = j_inv(E);
  
  return secret;
}


fp2 sike_isoex_3(const sike_public_key_t pkO, const sike_private_key_3 sk3I) {
  
  fp2 secret = { 0 };
  mont_curve_int_t E = {0};
  
  uint16_t e = PARE3;
  uint16_t msb = MSBB;

  mont_pt_t S = {0};
  E = get_yP_yQ_A_B(pkO);

  fp skI = sk3I;

  S = mont_double_and_add(E, skI, E.Q, msb);
  S = xADD(E, E.P, S);
  E.P = mont_inf_affine();
  E.Q = mont_inf_affine();
  E = iso_3_e((int)e, E, S);
  secret = j_inv(E);
  
  return secret;
}

// Pretty printing
void printFP(const fp z) {
  printf("0x%08x\n", z);
}

void printFP2(const fp2 z) {
  printf("{x0 = 0x%08x, x1 = 0x%08x}\n", z.x0, z.x1);
}

void printPoint(const mont_pt_t P) {
  printf("{x = {x0 = 0x%08x, x1 = 0x%08x}, y = {x0 = 0x%08x, x1 = 0x%08x}}\n",
         P.x.x0, P.x.x1, P.y.x0, P.y.x1);
}

void printCurve(const mont_curve_int_t C) {
  printf("{a = {x0 = 0x%08x, x1 = 0x%08x}, b = {x0 = 0x%08x, x1 = 0x%08x}}\n",
         C.a.x0, C.a.x1, C.b.x0, C.b.x1);
}

void printPublicKey(const sike_public_key_t K) {
  printf("%08x %08x, %08x %08x, %08x %08x\n", K.xP.x0, K.xP.x1, K.xQ.x0,
         K.xQ.x1, K.xR.x0, K.xR.x1);
}

void printOctetString(const uint8_t* os, size_t len) {
  printf("[");
  size_t i;
  for(i = 0; i < len - 1; i++) {
    printf("0x%02x, ", os[i]);
  }
  printf("0x%02x", os[i]);
  printf("]\n");
}


// Check functions
fp2 gen3ex2(sike_private_key_2 a, sike_private_key_3 b) {
  sike_params_t params = { 0 };
  params = sike_setup_params_C(sikeRawParams_C());
  sike_public_key_t pkB = { 0 };
  pkB = sike_isogen_3(params, b);
  return sike_isoex_2(pkB, a);
}

fp2 gen2ex3(sike_private_key_2 a, sike_private_key_3 b) {
  sike_params_t params = { 0 };
  params = sike_setup_params_C(sikeRawParams_C());
  sike_public_key_t pkA = { 0 };
  pkA = sike_isogen_2(params, a);
  return sike_isoex_3(pkA, b);
}

int secretShared(sike_private_key_2 a, sike_private_key_2 b) {
  fp2 c1 = { 0 }, c2 = { 0 };
  c1 = gen2ex3(a, b);
  c2 = gen3ex2(a, b);
  return fp2Cmp_C(c1, c2);
}

int main(int argc, char *argv[]) {
  sike_private_key_2 skA;
  sike_private_key_3 skB;
  time_t t;
  srand((unsigned) time(&t));
  if (argc == 1) {
    skA = rand();
    skB = rand();
  }
  else {
    skA = (uint16_t)strtol(argv[1], NULL, 10); // Range here is 0 to 32767.
    if (argc == 2) {
      skB = rand();
    }
    else {
      skB = (uint16_t)strtol(argv[2], NULL, 10);   // Range here is 0 to 6560, although in practice it's 0 to 4095.
    }
  }
  printf("skA: %d, skB: %d\n", skA, skB);
  fp2 c = gen2ex3(skA, skB);
  printf("gen2ex3: ");
  printFP2(c);
  printf("gen3ex2: ");
  fp2 d = gen3ex2(skA, skB);
  printFP2(d);
  return !fp2Cmp_C(c, d);
}
