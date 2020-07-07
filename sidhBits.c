#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MODULUS 0x0cd07fff


typedef uint32_t * mp;

void fp_Init (const mp p, mp* q) {
    *q = calloc(1, sizeof(*p));
}

void fp_Clear(const mp p, mp* q) {
    free(*q);
}

void fp_Add(const mp p, const mp a, const mp b, mp c) {
    *c = (uint32_t)( ((uint64_t)*a + (uint64_t)*b ) % (uint64_t)*p );
}

void fp_Constant(const mp p, unsigned long a, mp b) {
    *b = a % *p;
}

void fp_Copy(const mp p, mp dst, const mp src) {
    memcpy(dst, src, sizeof(*p));
}

int fp_IsEqual(const mp p, const mp a, const mp b) {
    return (*a % *p) == (*b % *p);
}

void fp_Unity(const mp p, mp b) {
    fp_Constant(p, 1, b);
}

void fp_Multiply(const mp p, const mp a, const mp b, mp c) {
    *c = (uint32_t)( ((uint64_t)*a * (uint64_t)*b ) % (uint64_t)*p);
}

void fp_Square(const mp p, const mp a, mp b) {
    fp_Multiply(p, a, a, b);
}

void fp_Negative(const mp p, const mp a, mp b) {
    mp t1, t2;
    fp_Init(p, &t1);
    fp_Constant(p, (*a % *p), t1);
    *b = (*t1 == 0) ? 0 : *p - (*a % *p);
    fp_Clear(p, &t1);
}

void fp_Subtract(const mp p, const mp a, const mp b, mp c) {
    mp t1;
    fp_Init(p, &t1);
    fp_Negative(p, b, t1);
    fp_Add(p, a, t1, c);
    fp_Clear(p, &t1);
}

int fp_IsBitSet(const mp p, const mp a, const unsigned long index) {
    return (*a & ( 1 << index )) >> index;
}

void fp_Pow(const mp p, const mp a, const mp b, mp c) {
    mp t1;
    fp_Init(p, &t1);
    fp_Unity(p, t1);
    for (int i = (sizeof(p)*4 - 1); i >= 0; i--)
    {                 
        fp_Square(p, t1, t1);

        if (fp_IsBitSet(p, b, i)) {
            fp_Multiply(p, t1, a, t1);
        }
    }
    fp_Copy(p, c, t1);
    fp_Clear(p, &t1);
}

void fp_Invert(const mp p, const mp a, mp b) {
    mp t1;
    fp_Init(p, &t1);
    fp_Constant(p, *p - 2, t1);
    fp_Pow(p, a, t1, b);
    fp_Clear(p, &t1);
}


int fp_IsConstant(const mp p, const mp a, const size_t constant) {
    return (*a % *p) == (constant % *p);
}

void fp_Zero(const mp p, mp b) {
    fp_Constant(p, 0, b);
}

int fp_IsEven(const mp p, const mp a) {
    return !fp_IsBitSet(p, a, 0);
}

void fp_Sqrt(const mp p, const mp a, mp b) {
    mp t1, p14;
    fp_Init(p, &p14);
    fp_Init(p, &t1);
    
    fp_Constant(p, (*p + 1) >> 2, p14);
    fp_Pow(p, a, p14, t1);
    if (fp_IsEven(p, t1)) {
        fp_Copy(p, b, t1);
    } else {
        fp_Negative(p, t1, t1);
        fp_Copy(p, b, t1);
    }
    fp_Clear(p, &p14);
    fp_Clear(p, &t1);

}

int fp_QuadNonRes(const mp p, const mp a) {
    mp t1, t2;
    fp_Init(p, &t1);
    fp_Init(p, &t2);
    fp_Sqrt(p, a, t1);
    fp_Square(p, t1, t1);
    if ( fp_IsEqual(p, a, t1) ) return 0 ; else return 1;
}





typedef struct {
  mp x0;
  mp x1;
} fp2;

void
fp2_Add( const mp p, const fp2* a, const fp2* b, fp2* c )
{
  fp_Add(p, a->x0, b->x0, c->x0);
  fp_Add(p, a->x1, b->x1, c->x1);
}

void fp2_Clear( const mp p, fp2* fp2)
{
  fp_Clear(p, &(fp2->x0));
  fp_Clear(p, &(fp2->x1));
}


void
fp2_Copy( const mp p, const fp2* a, fp2* b )
{
  fp_Copy( p, b->x0, a->x0 );
  fp_Copy( p, b->x1, a->x1 );
}


void
fp2_Init( const mp p, fp2* fp2 )
{
  fp_Init(p, &(fp2->x0));
  fp_Init(p, &(fp2->x1));
}

int
fp2_IsEqual( const mp p, const fp2* a1, const fp2* a2 )
{
  return fp_IsEqual(p, a1->x0, a2->x0) && fp_IsEqual(p, a1->x1, a2->x1);
}


void
fp2_Init_set( const mp p, fp2* fp2, unsigned long x0, unsigned long x1 )
{
  fp2_Init(p, fp2);
  fp_Constant(p, x0, fp2->x0);
  fp_Constant(p, x1, fp2->x1);
}


void
fp2_Set( const mp p, fp2* fp2, unsigned long x0, unsigned long x1 )
{
  fp_Constant(p, x0, fp2->x0);
  fp_Constant(p, x1, fp2->x1);
}


void
fp2_Invert( const mp p, const fp2* a, fp2* b )
{
  mp mul0;
  mp mul1;
  fp_Init(p, &mul0);
  fp_Init(p, &mul1);

  fp_Square( p, a->x0, mul0 );
  fp_Square( p, a->x1, mul1 );

  fp_Add( p, mul0, mul1, mul0 );
  fp_Invert( p, mul0, mul0 );

  fp_Negative( p, a->x1, mul1 );

  fp_Multiply( p, a->x0, mul0, b->x0 );
  fp_Multiply( p, mul1, mul0, b->x1 );

  fp_Clear(p, &mul0);
  fp_Clear(p, &mul1);
}


void
fp2_Multiply( const mp p, const fp2*  a, const fp2*  b, fp2*  c )
{
  mp mul0;
  mp mul1;
  mp adda;
  mp addb;

  fp_Init(p, &mul0);
  fp_Init(p, &mul1);
  fp_Init(p, &adda);
  fp_Init(p, &addb);

  fp_Multiply( p, a->x0, b->x0, mul0 );
  fp_Multiply( p, a->x1, b->x1, mul1 );

  fp_Add( p, a->x0, a->x1, adda );
  fp_Add( p, b->x0, b->x1, addb );

  fp_Subtract( p, mul0, mul1, c->x0 );

  fp_Add( p, mul0, mul1, mul0 );
  fp_Multiply( p, adda, addb, mul1 );

  fp_Subtract( p, mul1, mul0, c->x1 );

  fp_Clear(p, &mul0);
  fp_Clear(p, &mul1);
  fp_Clear(p, &adda);
  fp_Clear(p, &addb);
}

void
fp2_Negative( const mp p, const fp2* a, fp2* b )
{
  fp_Negative( p, a->x0, b->x0 );
  fp_Negative( p, a->x1, b->x1 );
}

void
fp2_Square( const mp p, const fp2* a, fp2* b )
{
  fp2_Multiply( p, a, a, b );
}

void
fp2_Sub( const mp p, const fp2* a, const fp2* b, fp2* c )
{
  fp_Subtract( p, a->x0, b->x0, c->x0 );
  fp_Subtract( p, a->x1, b->x1, c->x1 );
}


int
fp2_IsConst( const mp p, const fp2* a, unsigned long x0, unsigned long x1 ) {
  return fp_IsConstant(p, a->x0, x0) && fp_IsConstant(p, a->x1, x1);
}


void fp2_Sqrt(const mp p, const fp2* a, fp2* b)
{
    mp t0, t1;
    fp_Init(p, &t0);
    fp_Init(p, &t1);
    if (( fp_IsConstant(p, a->x1, 0) ) && (fp_QuadNonRes(p, a->x0))) {
            fp_Zero(p, t0);
            fp_Sqrt(p, a-> x0, t1);
            fp_Copy(p, b->x0, t0);
            fp_Copy(p, b->x1, t1);
    }
    else {
        mp t2, t3, p14, p34, inv2;
        fp_Init(p, &t2);
        fp_Init(p, &t3);
        fp_Init(p, &p14);
        fp_Init(p, &p34);
        fp_Init(p, &inv2);
        fp_Constant(p, 2, inv2);

        // (p + 1) / 4
        fp_Constant(p,(*p + 1) >> 2, p14);

        // (p - 3) / 4
        fp_Constant(p, (*p - 3) >> 2, p34);

        fp_Invert(p, inv2, inv2);
        fp_Square(p, a->x0, t0); 
        fp_Square(p, a->x1, t1); 
        fp_Add(p, t0, t1, t0); 
        fp_Pow(p, t0, p14, t1); 
        fp_Add(p, a->x0, t1, t0); 
        fp_Multiply(p, t0, inv2, t0);

        //p->half(p, t0);
        fp_Pow(p, t0, p34, t2); 
        //fp_Multiply(p, t0, t2, t1);

        fp_Pow(p, t0, p14, t1);
        fp_Multiply(p, t2, a->x1, t2); 
        fp_Multiply(p, t2, inv2, t2); 

        fp_Square(p, t1, t3);
        if (!fp_IsEqual(p, t3, t0)) {
            fp_Copy(p, t0, t1);
            fp_Copy(p, t1, t2);
            fp_Negative(p, t0, t2);
        }
        if (fp_IsEven(p, t1)) {
            fp_Copy(p, b->x0, t1);
            fp_Copy(p, b->x1, t2);
        } else {
            fp_Negative(p, t1, b->x0);
            fp_Negative(p, t2, b->x1);
        }
        fp_Clear(p, &t2);
        fp_Clear(p, &t3);
        fp_Clear(p, &p14);
        fp_Clear(p, &p34);
        fp_Clear(p, &inv2);
    }
    fp_Clear(p, &t0);
    fp_Clear(p, &t1);
  }



int main () {
    mp p = calloc(1, sizeof(uint32_t));
    *p =  MODULUS;
    /*mp i;
    fp_Init(p, &i);
    fp_Constant(p, 0x00000003, i);
    fp_Invert(p, i, i);*/
    
    fp2 i;
    fp2 j;
    fp2 k;
    fp2_Init(p, &i);
    fp2_Init(p, &j);
    fp2_Init(p, &k);
    fp2_Init_set(p, &i, 0x00000009, 0x00000000);//[0x083ac890, 0x0cbdc022]
    fp2_Sqrt(p, &i, &i);
    printf("[0x%08x, 0x%08x]\n", *i.x0, *i.x1);
    
   // printf("0x%08x\n", *i);
    return 0;
}    
