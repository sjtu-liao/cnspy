/* This file is part of the cnspy(1.0).
 * Copyright (C) 2023 Shanghai Jiao Tong University and Authors.
 * Contact: Bo Zhang <zilpher@sjtu.edu.cn>
 * Licensed under the <MPL-2.0>;
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at <https://www.mozilla.org/en-US/MPL/2.0/>.
 */

#include "lorenz63_fix.h"
#define _prec 1366
#define _order 617
int main(int argc, char *argv[]){
	tsm_init();
	tsm_clock_start();
	FILE *fp=fopen("..\ans\result_lorenz63_fix.txt","w");
	int i, j, order;
	order = _order;
	mpfr_t x_0[order+1], x_1[order+1], x_2[order+1];
	mpfr_t saver_0[order+1], saver_1[order+1], saver_2[order+1], saver_3[order+1], saver_4[order+1], saver_5[order+1], saver_6[order+1], saver_7[order+1], saver_8[order+1], saver_9[order+1], saver_10[order+1], saver_11[order+1];
	mpfr_t constant_0, constant_1, constant_2, constant_3, constant_4;
	mpfr_t t[order+1];
	mpfr_t inc_time, final_time;
	mpfr_t x_0_p, x_1_p, x_2_p;
	mpfr_inits2(_prec, inc_time, final_time, (mpfr_ptr) 0);
	mpfr_inits2(_prec, x_0_p, x_1_p, x_2_p, constant_0, constant_1, constant_2, constant_3, constant_4, (mpfr_ptr) 0);
	for(i=0; i<order+1; i++){
		mpfr_inits2(_prec, t[i], x_0[i], x_1[i], x_2[i], saver_0[i], saver_1[i], saver_2[i], saver_3[i], saver_4[i], saver_5[i], saver_6[i], saver_7[i], saver_8[i], saver_9[i], saver_10[i], saver_11[i], (mpfr_ptr) 0);
	}
	mpfr_set_str(final_time, "1000", 10, GMP_RNDN);
	mpfr_set_str(inc_time, "0.01", 10, GMP_RNDN);
	for(i=0; i<order+1; i++){
		mpfr_set_str(t[i], "0.0", 10, GMP_RNDN);
	}
	mpfr_set_str(constant_0, "-2.666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667", 10, GMP_RNDN);
	mpfr_set_str(constant_1, "28.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10, GMP_RNDN);
	mpfr_set_str(constant_2, "-10.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10, GMP_RNDN);
	mpfr_set_str(constant_3, "10.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10, GMP_RNDN);
	mpfr_set_str(constant_4, "-1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10, GMP_RNDN);
	mpfr_set_str(t[0], "0", 10, GMP_RNDN);
	mpfr_set_str(t[1], "1.0", 10, GMP_RNDN);
	mpfr_set_str(x_0_p, "-15.8", 10, GMP_RNDN);
	mpfr_set_str(x_1_p, "-17.48", 10, GMP_RNDN);
	mpfr_set_str(x_2_p, "35.64", 10, GMP_RNDN);
	int fprint_flag=1;
	mpfr_t print_time, print_time_inc, fake_time;
	mpfr_inits2(_prec, print_time, print_time_inc, fake_time, (mpfr_ptr) 0);
	mpfr_set_str(print_time, "0.1", 10, GMP_RNDN);
	mpfr_set_str(print_time_inc, "0.1", 10, GMP_RNDN);
	mpfr_set_str(fake_time, "0.0", 10, GMP_RNDN);
	mpfr_t tol, max_Mm1, max_M, abs_tmp_1, abs_tmp_2, CM, VSk1, VSk2, VSY1, VSY2;
	mpfr_inits2(_prec, tol, max_Mm1, max_M, abs_tmp_1, abs_tmp_2, CM, VSk1, VSk2, VSY1, VSY2, (mpfr_ptr) 0);
	mpfr_set_str(tol,"10.0",10,GMP_RNDN);
	mpfr_pow_si(tol, tol, -_prec, GMP_RNDN);
	while (mpfr_less_p(t[0], final_time)){
		mpfr_set(x_0[0],x_0_p,GMP_RNDN);
		mpfr_set(x_1[0],x_1_p,GMP_RNDN);
		mpfr_set(x_2[0],x_2_p,GMP_RNDN);
		/*------------Taylor Series Method------------*/
		for (i = 0; i < order; i++){
			tsm_mul_cs(constant_3,x_1,saver_0,i);
			tsm_mul_cs(constant_2,x_0,saver_1,i);
			tsm_add_ss(saver_0,saver_1,saver_2,i);
			tsm_mul_cs(constant_4,x_1,saver_3,i);
			tsm_mul_cs(constant_1,x_0,saver_4,i);
			tsm_mul_cs(constant_4,x_0,saver_5,i);
			tsm_mul_ss(saver_5,x_2,saver_6,i);
			tsm_add_ss(saver_3,saver_4,saver_7,i);
			tsm_add_ss(saver_7,saver_6,saver_8,i);
			tsm_mul_cs(constant_0,x_2,saver_9,i);
			tsm_mul_ss(x_0,x_1,saver_10,i);
			tsm_add_ss(saver_9,saver_10,saver_11,i);
			mpfr_div_ui(x_0[i+1],saver_2[i],(i+1),GMP_RNDN);
			mpfr_div_ui(x_1[i+1],saver_8[i],(i+1),GMP_RNDN);
			mpfr_div_ui(x_2[i+1],saver_11[i],(i+1),GMP_RNDN);
		}
		mpfr_abs(abs_tmp_1,x_0[order-1],GMP_RNDN);
		mpfr_abs(abs_tmp_2,x_1[order-1],GMP_RNDN);
		mpfr_max(max_Mm1, abs_tmp_1, abs_tmp_2, GMP_RNDN);
		mpfr_abs(abs_tmp_1,x_2[order-1],GMP_RNDN);
		mpfr_max(max_Mm1, max_Mm1, abs_tmp_1, GMP_RNDN);
		mpfr_abs(abs_tmp_1,x_0[order],GMP_RNDN);
		mpfr_abs(abs_tmp_2,x_1[order],GMP_RNDN);
		mpfr_max(max_M, abs_tmp_1, abs_tmp_2, GMP_RNDN);
		mpfr_abs(abs_tmp_1,x_2[order],GMP_RNDN);
		mpfr_max(max_M, max_M, abs_tmp_1, GMP_RNDN);
		mpfr_set_str(CM, "1.0", 10, GMP_RNDN);
		mpfr_div_si(CM, CM, order, GMP_RNDN);
		mpfr_pow(VSk1, tol, CM, GMP_RNDN);
		mpfr_set_str(CM, "-1.0", 10, GMP_RNDN);
		mpfr_div_si(VSY1, CM, order-1, GMP_RNDN);
		mpfr_set_str(CM, "1.0", 10, GMP_RNDN);
		mpfr_div_si(CM, CM, order+1, GMP_RNDN);
		mpfr_pow(VSk2, tol, CM, GMP_RNDN);
		mpfr_set_str(CM, "-1.0", 10, GMP_RNDN);
		mpfr_div_si(VSY2, CM, order, GMP_RNDN);
		mpfr_pow(max_Mm1, max_Mm1, VSY1, GMP_RNDN);
		mpfr_mul(max_Mm1, VSk1, max_Mm1, GMP_RNDN);
		mpfr_pow(max_M, max_M, VSY2, GMP_RNDN);
		mpfr_mul(max_M, VSk2, max_M, GMP_RNDN);
		mpfr_min(inc_time, max_Mm1, max_M, GMP_RNDN);
		fprint_flag = 0;
		mpfr_add(fake_time,t[0],inc_time,GMP_RNDN);
		if (mpfr_greater_p(fake_time, print_time)){
			mpfr_sub(inc_time, print_time, t[0], GMP_RNDN);
			mpfr_add(print_time, print_time, print_time_inc, GMP_RNDN);
			fprint_flag = 1;
		}
		tsm_sum2(x_0, inc_time, x_0_p, order);
		tsm_sum2(x_1, inc_time, x_1_p, order);
		tsm_sum2(x_2, inc_time, x_2_p, order);
		mpfr_add(t[0],t[0],inc_time,GMP_RNDN);
		if (fprint_flag){
		tsm_block_start(0);
		mpfr_fprintf(fp,"%.2Rf\t",t[0]);
		mpfr_fprintf(fp,"%1.500Re\t",x_0_p);
		mpfr_fprintf(fp,"%1.500Re\t",x_1_p);
		mpfr_fprintf(fp,"%1.500Re\n",x_2_p);
		tsm_block_end();
		}
	}
	tsm_clock_end();
	tsm_fina();
	return 0;
}
