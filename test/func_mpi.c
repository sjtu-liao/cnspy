#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"
#include <sys/time.h>
#include "mpi.h"
#include "mpf2mpfr.h"
#include "mpi_gmp.h"
#include "mpi_mpf.h"


#define tsm_init();\
    double st, et;\
    int _id, _size;\
    MPI_Init(&argc, &argv);\
    MPI_Comm_size(MPI_COMM_WORLD, &_size);\
    MPI_Comm_rank(MPI_COMM_WORLD, &_id);\
    mpfr_set_default_prec(_prec);\
    commit_mpf(&(MPI_MPF),_prec,MPI_COMM_WORLD);\
    create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);\
    mpfr_t _son, _father,_zero;\
    mpfr_inits2(_prec, _son, _father,_zero, (mpfr_ptr) 0);\
    void *packed_1, *packed_rec1,*packed_2, *packed_rec2; \
    packed_1 = allocbuf_mpf(_prec, 1);\
    packed_rec1 = allocbuf_mpf(_prec, 1);\
    packed_2 = allocbuf_mpf(_prec, 1);\
    packed_rec2 = allocbuf_mpf(_prec, 1);

#define tsm_fina();\
    mpfr_clears(_son, _father,_zero, NULL);\
    free_mpf_op(&(MPI_MPF_SUM));\
    free_mpf(&(MPI_MPF));\
    free(packed_1);\
    free(packed_rec1);\
    free(packed_2);\
    free(packed_rec2);\
    MPI_Finalize();

#define tsm_clock_start();\
    st=MPI_Wtime();

#define tsm_clock_end();\
    if(_id==0)\
    {\
        et=MPI_Wtime();\
        printf("CPU's time is %fs\n", et-st);\
    }

#define tsm_block_start(x);\
    if(_id==x){

#define tsm_block_end();\
    }

/**
 * @brief base function C = ln B
 * @param B (mpfr_t*) input series
 * @param C (mpfr_t*) output series of ln B
 * @param k (int) order of the present iteration
 */

#define tsm_log(B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_log(C[0], B[0], GMP_RNDN);\
    } \
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\        
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father, B[k], _father, GMP_RNDN);\
        mpfr_div(C[k], _father, B[0], GMP_RNDN);\
    }

/**
 * @brief base function C = exp(B)
 * @param B (mpfr_t*) input series
 * @param C (mpfr_t*) output series of exp(B)
 * @param k (int) order of the present iteration
 */

#define tsm_exp(B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_exp(C[0], B[0], GMP_RNDN);\
    } \
    else\
    {\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(C[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = tan(A) and  help function 1/cos(A)**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of 1/cos(A)**2
 * @param C (mpfr_t*) output series of tan(A)
 * @param k (int) order of the present iteration
 */

#define tsm_tan(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_tan(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], C[0], GMP_RNDN);\
        mpfr_add_si(B[0], B[0], 1, GMP_RNDN);\
    } \
    else\
    {\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(C[k],_father,GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, C[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(_father, _father, 2, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = tanh(A) and  help function 1/cosh(A)**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of 1/cosh(A)**2
 * @param C (mpfr_t*) output series of tanh(A)
 * @param k (int) order of the present iteration
 */

#define tsm_tanh(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_tanh(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], C[0], GMP_RNDN);\
        mpfr_si_sub(B[0], 1, B[0], GMP_RNDN);\
    } \
    else\
    {\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(C[k],_father,GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, C[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(_father, _father, -2, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = cot(A) and  help function -1/sin(A)**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -1/sin(A)**2
 * @param C (mpfr_t*) output series of cot(A)
 * @param k (int) order of the present iteration
 */

#define tsm_cot(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_cot(C[0], A[0], GMP_RNDN);\
        mpfr_csc(B[0],A[0],GMP_RNDN);\
        mpfr_sqr(B[0], B[0], GMP_RNDN);\
        mpfr_neg(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(C[k],_father,GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, C[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(_father, _father, -2, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = coth(A) and  help function -1/sinh(A)**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -1/sinh(A)**2
 * @param C (mpfr_t*) output series of coth(A)
 * @param k (int) order of the present iteration
 */

#define tsm_coth(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_coth(C[0], A[0], GMP_RNDN);\
        mpfr_csch(B[0],A[0],GMP_RNDN);\
        mpfr_sqr(B[0], B[0], GMP_RNDN);\
        mpfr_neg(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(C[k],_father,GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, C[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(_father, _father, -2, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = arcsin(A) and  help function sqrt(1-A**2)
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of sqrt(1-A**2)
 * @param C (mpfr_t*) output series of arcsin(A)
 * @param k (int) order of the present iteration
 */

#define tsm_asin(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_asin(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_si_sub(B[0],1,B[0],GMP_RNDN);\
        mpfr_sqrt(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=size)\
        {\
            mpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_neg(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = arcsinh(A) and  help function sqrt(1+A**2)
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of sqrt(1+A**2)
 * @param C (mpfr_t*) output series of arcsinh(A)
 * @param k (int) order of the present iteration
 */

#define tsm_asinh(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_asinh(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_add_si(B[0],B[0],1,GMP_RNDN);\
        mpfr_sqrt(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = arccos(A) and  help function -sqrt(1-A**2)
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -sqrt(1-A**2)
 * @param C (mpfr_t*) output series of arccos(A)
 * @param k (int) order of the present iteration
 */

#define tsm_acos(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_acos(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_si_sub(B[0],1,B[0],GMP_RNDN);\
        mpfr_sqrt(B[0],B[0],GMP_RNDN);\
        mpfr_neg(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_neg(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = arccosh(A) and  help function sqrt(A**2-1)
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of sqrt(A**2-1)
 * @param C (mpfr_t*) output series of arccosh(A)
 * @param k (int) order of the present iteration
 */

#define tsm_acosh(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_acosh(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_sub_si(B[0],B[0],1,GMP_RNDN);\
        mpfr_sqrt(B[0],B[0],GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_set(B[k],_father,GMP_RNDN);\
    }

/**
 * @brief base function C = arctan(A) and  help function A**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of A**2
 * @param C (mpfr_t*) output series of arctan(A)
 * @param k (int) order of the present iteration
 */

#define tsm_atan(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_atan(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_add_si(B[0], B[0], 1, GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(B[k],_father,2,GMP_RNDN);\
    }

/**
 * @brief base function C = arctanh(A) and  help function -A**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -A**2
 * @param C (mpfr_t*) output series of arctanh(A)
 * @param k (int) order of the present iteration
 */

#define tsm_atanh(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_atanh(C[0], A[0], GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_si_sub(B[0], 1, B[0], GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(B[k],_father,-2,GMP_RNDN);\
    }

/**
 * @brief base function C = arccot(A) and  help function -A**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -A**2
 * @param C (mpfr_t*) output series of arccot(A)
 * @param k (int) order of the present iteration
 */

#define tsm_acot(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_si_div(_son,1,A[0],GMP_RNDN);\
        mpfr_atan(C[0], _son, GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_si_sub(B[0], -1, B[0], GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(B[k],_father,-2,GMP_RNDN);\
    }

/**
 * @brief base function C = arccoth(A) and  help function -A**2
 * @param A (mpfr_t*) input series
 * @param B (mpfr_t*) output series of -A**2
 * @param C (mpfr_t*) output series of arccoth(A)
 * @param k (int) order of the present iteration
 */

#define tsm_acoth(A,B,C,k);\
    mpfr_set_str(_son, "0.0", 10, GMP_RNDN);\
    mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
    if (k == 0) {\
        mpfr_si_div(_son,1,A[0],GMP_RNDN);\
        mpfr_atanh(C[0], _son, GMP_RNDN);\
        mpfr_sqr(B[0], A[0], GMP_RNDN);\
        mpfr_si_sub(B[0], 1, B[0], GMP_RNDN);\
    }\
    else\
    {\
        for (int _i = _id+1; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, C[_i], _i, GMP_RNDN);\
            mpfr_mul(_son, B[k-_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_sub(_father,A[k],_father,GMP_RNDN);\
        mpfr_div(C[k],_father,B[0],GMP_RNDN);\
        mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
        for (int _i = _id; _i < k; _i+=_size)\
        {\
            mpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);\
            mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
            mpfr_add(_father,_father,_son,GMP_RNDN);\
        }\
        pack_mpf(_father, 1, packed_1);\
        MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
        MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
        unpack_mpf(packed_rec1, _father, 1);\     
        mpfr_div_si(_father, _father, k, GMP_RNDN);\
        mpfr_mul_si(B[k],_father,-2,GMP_RNDN);\
    }
