/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) chain_nm_5_external_cost_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const casadi_int casadi_s0[28] = {24, 1, 0, 24, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
static const casadi_int casadi_s1[7] = {3, 1, 0, 3, 0, 1, 2};
static const casadi_int casadi_s2[31] = {27, 1, 0, 27, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
static const casadi_int casadi_s3[57] = {27, 27, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};

/* chain_nm_5_external_cost:(i0[24],i1[3])->(o0[27],o1[27x27,27nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a2;
  a0=5.0000000000000000e-01;
  a1=arg[1] ? arg[1][0] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][0]=a1;
  a1=arg[1] ? arg[1][1] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][1]=a1;
  a1=arg[1] ? arg[1][2] : 0;
  a1=(a1+a1);
  a0=(a0*a1);
  if (res[0]!=0) res[0][2]=a0;
  a0=5.0000000000000001e-03;
  a1=arg[0] ? arg[0][0] : 0;
  a2=2.4547359999999999e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][3]=a1;
  a1=arg[0] ? arg[0][1] : 0;
  a2=2.8850290000000002e-18;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][4]=a1;
  a1=arg[0] ? arg[0][2] : 0;
  a2=-3.6941920000000000e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][5]=a1;
  a1=arg[0] ? arg[0][3] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][6]=a1;
  a1=arg[0] ? arg[0][4] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][7]=a1;
  a1=arg[0] ? arg[0][5] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][8]=a1;
  a1=arg[0] ? arg[0][6] : 0;
  a2=5.0498149999999997e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][9]=a1;
  a1=arg[0] ? arg[0][7] : 0;
  a2=5.9350009999999999e-18;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][10]=a1;
  a1=arg[0] ? arg[0][8] : 0;
  a2=-4.2382419999999998e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][11]=a1;
  a1=arg[0] ? arg[0][9] : 0;
  a2=-6.8034679999999999e-31;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][12]=a1;
  a1=arg[0] ? arg[0][10] : 0;
  a2=1.2876850000000001e-31;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][13]=a1;
  a1=arg[0] ? arg[0][11] : 0;
  a2=7.5782790000000005e-32;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][14]=a1;
  a1=arg[0] ? arg[0][12] : 0;
  a2=7.5454920000000003e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][15]=a1;
  a1=arg[0] ? arg[0][13] : 0;
  a2=8.8681460000000005e-18;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][16]=a1;
  a1=arg[0] ? arg[0][14] : 0;
  a2=-1.5288599999999999e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][17]=a1;
  a1=arg[0] ? arg[0][15] : 0;
  a2=-1.1907460000000000e-40;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][18]=a1;
  a1=arg[0] ? arg[0][16] : 0;
  a2=2.6957830000000002e-41;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][19]=a1;
  a1=arg[0] ? arg[0][17] : 0;
  a2=1.0175350000000000e-41;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][20]=a1;
  a1=arg[0] ? arg[0][18] : 0;
  a2=9.9453820000000004e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][21]=a1;
  a1=arg[0] ? arg[0][19] : 0;
  a2=1.1688710000000000e-17;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][22]=a1;
  a1=arg[0] ? arg[0][20] : 0;
  a2=4.1850540000000003e-01;
  a1=(a1-a2);
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][23]=a1;
  a1=arg[0] ? arg[0][21] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][24]=a1;
  a1=arg[0] ? arg[0][22] : 0;
  a1=(a1+a1);
  a1=(a0*a1);
  if (res[0]!=0) res[0][25]=a1;
  a1=arg[0] ? arg[0][23] : 0;
  a1=(a1+a1);
  a0=(a0*a1);
  if (res[0]!=0) res[0][26]=a0;
  a0=1.;
  if (res[1]!=0) res[1][0]=a0;
  if (res[1]!=0) res[1][1]=a0;
  if (res[1]!=0) res[1][2]=a0;
  a0=1.0000000000000000e-02;
  if (res[1]!=0) res[1][3]=a0;
  if (res[1]!=0) res[1][4]=a0;
  if (res[1]!=0) res[1][5]=a0;
  if (res[1]!=0) res[1][6]=a0;
  if (res[1]!=0) res[1][7]=a0;
  if (res[1]!=0) res[1][8]=a0;
  if (res[1]!=0) res[1][9]=a0;
  if (res[1]!=0) res[1][10]=a0;
  if (res[1]!=0) res[1][11]=a0;
  if (res[1]!=0) res[1][12]=a0;
  if (res[1]!=0) res[1][13]=a0;
  if (res[1]!=0) res[1][14]=a0;
  if (res[1]!=0) res[1][15]=a0;
  if (res[1]!=0) res[1][16]=a0;
  if (res[1]!=0) res[1][17]=a0;
  if (res[1]!=0) res[1][18]=a0;
  if (res[1]!=0) res[1][19]=a0;
  if (res[1]!=0) res[1][20]=a0;
  if (res[1]!=0) res[1][21]=a0;
  if (res[1]!=0) res[1][22]=a0;
  if (res[1]!=0) res[1][23]=a0;
  if (res[1]!=0) res[1][24]=a0;
  if (res[1]!=0) res[1][25]=a0;
  if (res[1]!=0) res[1][26]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int chain_nm_5_external_cost(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void chain_nm_5_external_cost_incref(void) {
}

CASADI_SYMBOL_EXPORT void chain_nm_5_external_cost_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int chain_nm_5_external_cost_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int chain_nm_5_external_cost_n_out(void) { return 2;}

CASADI_SYMBOL_EXPORT const char* chain_nm_5_external_cost_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* chain_nm_5_external_cost_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* chain_nm_5_external_cost_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* chain_nm_5_external_cost_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s3;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int chain_nm_5_external_cost_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
