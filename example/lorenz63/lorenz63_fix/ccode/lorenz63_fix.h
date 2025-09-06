/* This file is part of the cnspy(1.0).
 * Copyright (C) 2023 Shanghai Jiao Tong University and Authors.
 * Contact: Bo Zhang <zilpher@sjtu.edu.cn>
 * Licensed under the <MPL-2.0>;
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at <https://www.mozilla.org/en-US/MPL/2.0/>.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "gmp.h"
#include "mpfr.h"
#include "mpi.h"
#include "mpf2mpfr.h"
#include "mpi_mpf.h"

#define tsm_init();\
	double st, et;\
	int _i, _j,_k;\
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

#define tsm_change_prec(PREC);\
	mpfr_set_default_prec(PREC);\
	commit_mpf(&(MPI_MPF), PREC, MPI_COMM_WORLD);\
	create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);\
	mpfr_prec_round(_son,PREC,GMP_RNDN);\
	mpfr_prec_round(_father,PREC,GMP_RNDN);\
	mpfr_prec_round(_zero,PREC,GMP_RNDN);\
	packed_1 = allocbuf_mpf(PREC, 1);\
	packed_rec1 = allocbuf_mpf(PREC, 1);\
	packed_2 = allocbuf_mpf(PREC, 1);\
	packed_rec2 = allocbuf_mpf(PREC, 1);

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
		printf("CPU's time is %fs", et-st);\
	}

#define tsm_block_start(x);\
	if(_id==x){
#define tsm_block_end();\
	}

#define tsm_mul_cc(a,b,c,k);\
	mpfr_mul(c,a,b,GMP_RNDN);

#define tsm_mul_cs(a,B,C,k);\
	mpfr_mul(C[k],a,B[k],GMP_RNDN);

#define tsm_mul_sc(A,b,C,k);\
	mpfr_mul(C[k],A[k],b,GMP_RNDN);

#define tsm_mul_ss(A,B,C,k);\
	mpfr_set_str(_son,"0.0",10,GMP_RNDN);\
	mpfr_set_str(_father,"0.0",10,GMP_RNDN);\
	for (_i = _id; _i < k+1; _i+=_size)\
	{\
		mpfr_mul(_son,A[_i],B[k-_i],GMP_RNDN);\
		mpfr_add(_father,_father,_son,GMP_RNDN);\
	}\
	pack_mpf(_father, 1, packed_1);\
	MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
	MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
	unpack_mpf(packed_rec1, _father, 1);\
	mpfr_set(C[k],_father,GMP_RNDN);

#define tsm_add_cc(a,b,c,k);\
	mpfr_add(c,a,b,GMP_RNDN);

#define tsm_add_cs(a,B,C,k);\
	if(k==0)\
	{\
		mpfr_add(C[0],a,B[0],GMP_RNDN);\
	}\
	else\
	{\
		mpfr_set(C[k],B[k],GMP_RNDN);\
	}

#define tsm_add_sc(A,b,C,k);\
	if(k==0)\
	{\
		mpfr_add(C[0],A[0],b,GMP_RNDN);\
	}\
	else\
	{\
		mpfr_set(C[k],A[k],GMP_RNDN);\
	}

#define tsm_add_ss(A,B,C,k);\
	mpfr_add(C[i],A[i],B[i],GMP_RNDN);

#define tsm_sum2(A,b,C,k);\
	mpfr_set_str(_father, "0.0", 10, GMP_RNDN);\
	for (_i = _id; _i < k+1; _i+=_size)\
	{\
		mpfr_set_str(_son,"1.0",10,GMP_RNDN);\
		for(_j=1;_j<_i+1;_j++)\
			mpfr_mul(_son,_son,b,GMP_RNDN);\
		mpfr_mul(_son, A[_i], _son, GMP_RNDN);\
		mpfr_add(_father, _father, _son, GMP_RNDN);\
	}\
	pack_mpf(_father, 1, packed_1);\
	MPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);\
	MPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);\
	unpack_mpf(packed_rec1, _father, 1);\
	mpfr_set(C,_father,GMP_RNDN);
