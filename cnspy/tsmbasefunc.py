import sympy as sp
import numpy as np

class tsm_base_func():
    """
    Father class of basic functions of taylor series method.
    
    There are different types of each base function, constant 
    with constant(cc), constants with series(cs and sc), series 
    with series().
    
    The basic functions that can be used at present are: 
    `Add`, `Mul`, `Pow`, `Sin`, `Cos`, `Sinh`, `Cosh`.

    Yields
    ======
    Strings to .h file or .c file.
    
    """    
    def __init__(self):
        pass

    def add_h(self,mpi_flag):
        '''
        Basic function of `Add` to .h file.
        -----------------------------------
        if c(t) = a(t) + b(t):
        ss: c[k] = a[k] + b[k];
        cs: c[0] = a + b[0]; c[k] = b[k];
        sc: c[0] = a[0] + b; c[k] = a[k];
        cc: c = a + b;
        '''
        if mpi_flag:
            type_dict = {'cc': "\\\n".join(["#define tsm_add_cc(a,b,c,k);","\tmpfr_add(c,a,b,GMP_RNDN);"])+"\n",
                         'cs': "\\\n".join(["#define tsm_add_cs(a,B,C,k);","\tif(k==0)","\t{","\t\tmpfr_add(C[0],a,B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tmpfr_set(C[k],B[k],GMP_RNDN);","\t}"])+"\n",
                         'sc': "\\\n".join(["#define tsm_add_sc(A,b,C,k);","\tif(k==0)","\t{","\t\tmpfr_add(C[0],A[0],b,GMP_RNDN);","\t}","\telse","\t{","\t\tmpfr_set(C[k],A[k],GMP_RNDN);","\t}"])+"\n",
                         'ss': "\\\n".join(["#define tsm_add_ss(A,B,C,k);","\tmpfr_add(C[i],A[i],B[i],GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'cc': "\\\n".join(["#define tsm_add_cc(a,b,c,k);","\tmpfr_add(c,a,b,GMP_RNDN);"])+"\n",
                         'cs': "\\\n".join(["#define tsm_add_cs(a,B,C,k);","\tif(k==0)","\t{","\t\tmpfr_add(C[0],a,B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tmpfr_set(C[k],B[k],GMP_RNDN);","\t}"])+"\n",
                         'sc': "\\\n".join(["#define tsm_add_sc(A,b,C,k);","\tif(k==0)","\t{","\t\tmpfr_add(C[0],A[0],b,GMP_RNDN);","\t}","\telse","\t{","\t\tmpfr_set(C[k],A[k],GMP_RNDN);","\t}"])+"\n",
                         'ss': "\\\n".join(["#define tsm_add_ss(A,B,C,k);","\tmpfr_add(C[i],A[i],B[i],GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def mul_h(self,mpi_flag):
        '''
        Basic function of `Mul` to .h file.
        -----------------------------------
        if c(t) = a(t) * b(t):
        ss: c[k] = sum{j=0->k} a[j] * b[k-j];
        cs: c[k] = a * b[k];
        sc: c[k] = b * a[k];
        cc: c = a * b;
        '''
        if mpi_flag:
            type_dict = {'cc': "\\\n".join(["#define tsm_mul_cc(a,b,c,k);","\tmpfr_mul(c,a,b,GMP_RNDN);"])+"\n",
                         'cs': "\\\n".join(["#define tsm_mul_cs(a,B,C,k);","\tmpfr_mul(C[k],a,B[k],GMP_RNDN);"])+"\n",
                         'sc': "\\\n".join(["#define tsm_mul_sc(A,b,C,k);","\tmpfr_mul(C[k],A[k],b,GMP_RNDN);"])+"\n",
                         'ss': "\\\n".join(["#define tsm_mul_ss(A,B,C,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tfor (_i = _id; _i < k+1; _i+=_size)","\t{","\t\tmpfr_mul(_son,A[_i],B[k-_i],GMP_RNDN);","\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t}","\tpack_mpf(_father, 1, packed_1);","\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\tunpack_mpf(packed_rec1, _father, 1);","\tmpfr_set(C[k],_father,GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'cc': "\\\n".join(["#define tsm_mul_cc(a,b,c,k);","\tmpfr_mul(c,a,b,GMP_RNDN);"])+"\n",
                         'cs': "\\\n".join(["#define tsm_mul_cs(a,B,C,k);","\tmpfr_mul(C[k],a,B[k],GMP_RNDN);"])+"\n",
                         'sc': "\\\n".join(["#define tsm_mul_sc(A,b,C,k);","\tmpfr_mul(C[k],A[k],b,GMP_RNDN);"])+"\n",
                         'ss': "\\\n".join(["#define tsm_mul_ss(A,B,C,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tfor(_i = 0;_i < k+1;_i++)","\t{","\t\tmpfr_mul(_son,A[_i],B[k-_i],GMP_RNDN);","\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t}","\tmpfr_set(C[k],_father,GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def pow_h(self,mpi_flag):
        '''
        Basic function of `Pow` to .h file.
        -----------------------------------
        if c(t) = b(t)^a:
        s: c[k] = sum{j=0->k-1} (k*a-j(a+1))*c[j]*b[k-j] / k*b[0];
        c: c = b^a;
        '''
        if mpi_flag:
            type_dict = {'c': "\\\n".join(["#define tsm_pow_c(a,b,c,k);","\tmpfr_pow(c,a,b,GMP_RNDN)"])+"\n",
                         's': "\\\n".join(["#define tsm_pow(A,b,C,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_zero,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_pow(_father,A[0],b,GMP_RNDN);","\t}","\telse","\t{","\t\tif(!mpfr_equal_p(A[0],_zero))","\t\t{","\t\t\tif (_id==0)","\t\t\t{","\t\t\t\tmpfr_mul_si(_father,b,k,GMP_RNDN);","\t\t\t\tmpfr_mul(_father,_father,C[0],GMP_RNDN);","\t\t\t\tmpfr_mul(_father,_father,A[k],GMP_RNDN);","\t\t\t}","\t\t\tfor(_i=_id+1;_i<k;_i+=_size)","\t\t\t{","\t\t\t\tmpfr_mul_si(_son,b,k-_i,GMP_RNDN);","\t\t\t\tmpfr_sub_si(_son,_son,_i,GMP_RNDN);","\t\t\t\tmpfr_mul(_son,_son,C[_i],GMP_RNDN);","\t\t\t\tmpfr_mul(_son,_son,A[k-_i],GMP_RNDN);","\t\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t\t}","\t\t\tpack_mpf(_father, 1, packed_1);","\t\t\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\t\t\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\t\t\tunpack_mpf(packed_rec1, _father, 1);","\t\t\tmpfr_mul_si(_son,A[0],k,GMP_RNDN);","\t\t\tmpfr_div(_father,_father,_son,GMP_RNDN);","\t\t}","\t\telse","\t\t{","\t\t\tmpfr_set(_father,_zero,GMP_RNDN);","\t\t}","\t}","\tmpfr_set(C[k],_father,GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'c': "\\\n".join(["#define tsm_pow_c(a,b,c,k);","\tmpfr_pow(c,a,b,GMP_RNDN)"])+"\n",
                         's': "\\\n".join(["#define tsm_pow(A,b,C,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_zero,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_pow(_father,A[0],b,GMP_RNDN);","\t}","\telse","\t{","\t\tif(!mpfr_equal_p(A[0],_zero))","\t\t{","\t\t\tmpfr_mul_si(_father,b,k,GMP_RNDN);","\t\t\tmpfr_mul(_father,_father,C[0],GMP_RNDN);","\t\t\tmpfr_mul(_father,_father,A[k],GMP_RNDN);","\t\t\tfor(_i=1;_i<k;_i++)","\t\t\t{","\t\t\t\tmpfr_mul_si(_son,b,k-_i,GMP_RNDN);","\t\t\t\tmpfr_sub_si(_son,_son,_i,GMP_RNDN);","\t\t\t\tmpfr_mul(_son,_son,C[_i],GMP_RNDN);","\t\t\t\tmpfr_mul(_son,_son,A[k-_i],GMP_RNDN);","\t\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t\t}","\t\t\tmpfr_mul_si(_son,A[0],k,GMP_RNDN);","\t\t\tmpfr_div(_father,_father,_son,GMP_RNDN);","\t\t}","\t\telse","\t\t{","\t\t\tmpfr_set(_father,_zero,GMP_RNDN);","\t\t}","\t}","\tmpfr_set(C[k],_father,GMP_RNDN);"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def ln_h(self,mpi_flag):
        '''
        Basic function of `ln` to .h file.
        -----------------------------------
        if c(t) = ln b(t):
        s: c[k] = (k*b[k]-sum{j=1->k-1} ja[j]b[k-j]) / k*b[0];
        c: c = ln b;
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_ln_c(B,C,k);","\tmpfr_log(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_ln(B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,     GMP_RNDN);","\tif (k == 0) {","\t\t mpfr_log(C[0], B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father, 1, packed_1);","\t\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1, _father, 1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father, B[k], _father, GMP_RNDN);","\t\tmpfr_div(C[k], _father, B[0], GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_ln_c(B,C,k);","\tmpfr_log(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_ln(B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,     GMP_RNDN);","\tif (k == 0) {","\t\t mpfr_log(C[0], B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father, B[k], _father, GMP_RNDN);","\t\tmpfr_div(C[k], _father, B[0], GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def exp_h(self,mpi_flag):
        '''
        Basic function of `exp` to .h file.
        -----------------------------------
        if c(t) = exp(b(t)):
        s: c[k] = sum{j=0->k-1} (k-j)a[j]b[k-j] / k;
        c: c = exp(b);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_exp_c(B,C,k);","\tmpfr_exp(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_exp(B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,     GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_exp(C[0], B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\t\tpack_mpf(_father, 1, packed_1);","\t\t\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\t\t\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\t\t\tunpack_mpf(packed_rec1, _father, 1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_exp_c(B,C,k);","\tmpfr_exp(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_exp(B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,     GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_exp(C[0], B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def sin_h(self,mpi_flag):
        '''
        Basic function of `Sin` and `Cos` to .h file.
        -----------------------------------
        if c(t) = cos(a(t)) and s(t) = sin(a(t)):
        s: s[k] = sum{j=0->k-1} c[j]*(k-j)*a[k-j]/k;
           c[k] = sum{j=0->k-1} (-1)^k*s[j]*(k-j)*a[k-j]/k;
        c: c = cos(a);
           s = sin(a);
        '''
        if mpi_flag:
            type_dict = {'c': "\\\n".join(["#define tsm_sin_c(a,b,k);","\tmpfr_sin(b,a,GMP_RNDN);"])+"\n",
                         's': "\\\n".join(["#define tsm_sin_cos(x,sinx,cosx,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_sin(sinx[0],x[0],GMP_RNDN);","\t\tmpfr_cos(cosx[0],x[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (_i =_id; _i<k;_i+=_size)","\t\t{","\t\t\tmpfr_mul(_son,cosx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father, 1, packed_1);","\t\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1, _father, 1);","\t\tmpfr_div_ui(sinx[k], _father,k,GMP_RNDN);","\t\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\t\tfor (_i=_id;_i<k;_i+=_size)","\t\t{","\t\t\tmpfr_mul(_son,sinx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father, 1, packed_1);","\t\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1, _father, 1);","\t\tmpfr_div_ui(cosx[k],_father,k,GMP_RNDN);","\t\tmpfr_neg(cosx[k],cosx[k],GMP_RNDN);","\t}"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'c': "\\\n".join(["#define tsm_sin_c(a,b,k);","\tmpfr_sin(b,a,GMP_RNDN);"])+"\n",
                         's': "\\\n".join(["#define tsm_sin_cos(x,sinx,cosx,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_sin(sinx[0],x[0],GMP_RNDN);","\t\tmpfr_cos(cosx[0],x[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (_i=0;_i<k;_i++)","\t\t{","\t\t\tmpfr_mul(_son,cosx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_ui(sinx[k], _father,k,GMP_RNDN);","\t\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\t\tfor (_i=0;_i<k;_i++)","\t\t{","\t\t\tmpfr_mul(_son,sinx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_ui(cosx[k],_father,k,GMP_RNDN);","\t\tmpfr_neg(cosx[k],cosx[k],GMP_RNDN);","\t}"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        
    def cos_h(self,mpi_flag):
        '''
        This function only contains constant `Cos` function.
        Other functions are in the self.sin_h(self,mpi_flag).
        '''
        if mpi_flag:
            type_dict = {'c': "\\\n".join(["#define tsm_cos_c(a,b,k);","\tmpfr_cos(b,a,GMP_RNDN);"])+"\n",
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'c': "\\\n".join(["#define tsm_cos_c(a,b,k);","\tmpfr_cos(b,a,GMP_RNDN);"])+"\n",
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        
    def sinh_h(self,mpi_flag):
        '''
        Basic function of `Sinh` and `Cosh` to .h file.
        -----------------------------------
        if c(t) = cosh(a(t)) and s(t) = sinh(a(t)):
        s: s[k] = sum{j=0->k-1} c[j]*(k-j)*a[k-j]/k;
           c[k] = sum{j=0->k-1} s[j]*(k-j)*a[k-j]/k;
        c: c = cosh(a);
           s = sinh(a);
        '''
        if mpi_flag:
            type_dict = {'c': "\\\n".join(["#define tsm_sin_c(a,b,k);","\tmpfr_sin(b,a,GMP_RNDN);"])+"\n",
                         's': "\\\n".join(["#define tsm_sinh_cosh(x,sinhx,coshx,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_sinh(sinhx[0],x[0],GMP_RNDN);","\t\tmpfr_cosh(coshx[0],x[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (_i=_id;_i<k;_i+=_size)","\t\t{","\t\t\tmpfr_mul(_son,coshx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_ui(sinhx[k],_father,k,GMP_RNDN);","\t\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\t\tfor (_i=0;_i<k;_i++)","\t\t{","\t\t\tmpfr_mul(_son,sinhx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_ui(coshx[k],_father,k,GMP_RNDN);","\t}"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'c': "\\\n".join(["#define tsm_sin_c(a,b,k);","\tmpfr_sin(b,a,GMP_RNDN);"])+"\n",
                         's': "\\\n".join(["#define tsm_sinh_cosh(x,sinhx,coshx,k);","\tmpfr_set_str(_son,\"0.0\",10,GMP_RNDN);","\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\tif (k==0)","\t{","\t\tmpfr_sinh(sinhx[0],x[0],GMP_RNDN);","\t\tmpfr_cosh(coshx[0],x[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (_i=0;_i<k;_i++)","\t\t{","\t\t\tmpfr_mul(_son,coshx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_ui(sinhx[k],_father,k,GMP_RNDN);","\t\tmpfr_set_str(_father,\"0.0\",10,GMP_RNDN);","\t\tfor (_i=0;_i<k;_i++)","\t\t{","\t\t\tmpfr_mul(_son,sinhx[_i],x[k-_i],GMP_RNDN);","\t\t\tmpfr_mul_ui(_son,_son,(k-_i),GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_ui(coshx[k],_father,k,GMP_RNDN);","\t}"])+"\n"
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        
    def cosh_h(self,mpi_flag):
        '''
        This function only contains constant `Sin` function.
        Other functions are in the self.sinh_h(self,mpi_flag).
        '''
        if mpi_flag:
            type_dict = {'c': "\\\n".join(["#define tsm_cosh_c(a,b,k);","\tmpfr_cosh(b,a,GMP_RNDN);"])+"\n",
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'c': "\\\n".join(["#define tsm_cosh_c(a,b,k);","\tmpfr_cosh(b,a,GMP_RNDN);"])+"\n",
                        }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def tan_h(self,mpi_flag):
        '''
        Basic function of `tan` to .h file.
        -----------------------------------
        if c(t) = tan(a(t)) and  help function b(t)=1/cos(a(t))**2:
        s: c[k] = sum{j=0->k-1} b[j]*(k-j)*a[k-j]/k;
           b[t] = 2*sum{j=0->k-1} c[j]*(k-j)*c[k-j]/k;
        c: c = tan(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_tan_c(B,C,k);","\tmpfr_tan(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_tan(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_tan(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], C[0], GMP_RNDN);","\t\tmpfr_add_si(B[0], B[0], 1, GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, 2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_tan_c(B,C,k);","\tmpfr_tan(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_tan(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_tan(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], C[0], GMP_RNDN);","\t\tmpfr_add_si(B[0], B[0], 1, GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, 2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def tanh_h(self,mpi_flag):
        '''
        Basic function of `tanh` to .h file.
        -----------------------------------
        if c(t) = tanh(a(t)) and  help function b(t)=1/cosh(a(t))**2:
        s: c[k] = sum{j=0->k-1} b[j]*(k-j)*a[k-j]/k;
           b[t] = -2*sum{j=0->k-1} c[j]*(k-j)*c[k-j]/k;
        c: c = tanh(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_tanh_c(B,C,k);","\tmpfr_tanh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_tanh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_tanh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], C[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_tanh_c(B,C,k);","\tmpfr_tanh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_tanh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_tanh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], C[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
   
    def cot_h(self,mpi_flag):
        '''
        Basic function of `cot` to .h file.
        -----------------------------------
        if c(t) = cot(a(t)) and  help function b(t)=1/sin(a(t))**2:
        s: c[k] = sum{j=0->k-1} b[j]*(k-j)*a[k-j]/k;
           b[t] = -2*sum{j=0->k-1} c[j]*(k-j)*c[k-j]/k;
        c: c = cot(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_cot_c(B,C,k);","\tmpfr_cot(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_cot(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_cot(C[0], A[0], GMP_RNDN);","\t\tmpfr_csc(B[0],A[0],GMP_RNDN);","\t\tmpfr_sqr(B[0], B[0], GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_cot_c(B,C,k);","\tmpfr_cot(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_cot(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_cot(C[0], A[0], GMP_RNDN);","\t\tmpfr_csc(B[0],A[0],GMP_RNDN);","\t\tmpfr_sqr(B[0], B[0], GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"

    def coth_h(self,mpi_flag):        
        '''
        Basic function of `coth` to .h file.
        -----------------------------------
        if c(t) = cot(a(t)) and  help function b(t)=1/sinh(a(t))**2:
        s: c[k] = sum{j=0->k-1} b[j]*(k-j)*a[k-j]/k;
           b[t] = -2*sum{j=0->k-1} c[j]*(k-j)*c[k-j]/k;
        c: c = coth(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_coth_c(B,C,k);","\tmpfr_coth(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_coth(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_coth(C[0], A[0], GMP_RNDN);","\t\tmpfr_csch(B[0],A[0],GMP_RNDN);","\t\tmpfr_sqr(B[0], B[0], GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_coth_c(B,C,k);","\tmpfr_coth(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_coth(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_coth(C[0], A[0], GMP_RNDN);","\t\tmpfr_csch(B[0],A[0],GMP_RNDN);","\t\tmpfr_sqr(B[0], B[0], GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, B[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(C[k],_father,GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, C[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(_father, _father, -2, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def asin_h(self,mpi_flag):
        '''
        Basic function of `asin` to .h file.
        -----------------------------------
        if c(t) = asin(a(t)) and  help function b(t)=sqrt(1-a(t)**2):
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*b[0];
           b[t] = -sum{j=0->k-1} a[j]*(k-j)*c[k-j]/k;
        c: c = asin(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_asin_c(B,C,k);","\tmpfr_asin(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_asin(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_asin(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0],1,B[0],GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_neg(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_asin_c(B,C,k);","\tmpfr_asin(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_asin(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_asin(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0],1,B[0],GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_neg(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def asinh_h(self,mpi_flag):
        '''
        Basic function of `asinh` to .h file.
        -----------------------------------
        if c(t) = asinh(a(t)) and  help function b(t)=sqrt(1+a(t)**2):
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*b[0];
           b[t] = sum{j=0->k-1} a[j]*(k-j)*c[k-j]/k;
        c: c = asinh(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_asinh_c(B,C,k);","\tmpfr_asinh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_asinh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_asinh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_add_si(B[0],B[0],1,GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_asinh_c(B,C,k);","\tmpfr_asinh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_asinh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_asinh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_add_si(B[0],B[0],1,GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def acos_h(self,mpi_flag):        
        '''
        Basic function of `acos` to .h file.
        -----------------------------------
        if c(t) = acos(a(t)) and  help function b(t)=-sqrt(1-a(t)**2):
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*b[0];
           b[t] = -sum{j=0->k-1} a[j]*(k-j)*c[k-j]/k;
        c: c = acos(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_acos_c(B,C,k);","\tmpfr_acos(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acos(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_acos(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0],1,B[0],GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_neg(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_acos_c(B,C,k);","\tmpfr_acos(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acos(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_acos(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0],1,B[0],GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t\tmpfr_neg(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_neg(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def acosh_h(self,mpi_flag):
        '''
        Basic function of `acosh` to .h file.
        -----------------------------------
        if c(t) = acosh(a(t)) and  help function b(t)=sqrt(1-a(t)**2):
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*b[0];
           b[t] = sum{j=0->k-1} a[j]*(k-j)*c[k-j]/k;
        c: c = acosh(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_acosh_c(B,C,k);","\tmpfr_acosh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acosh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_acosh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_sub_si(B[0],B[0],1,GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_acosh_c(B,C,k);","\tmpfr_acosh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acosh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_acosh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_sub_si(B[0],B[0],1,GMP_RNDN);","\t\tmpfr_sqrt(B[0],B[0],GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_set(B[k],_father,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def atan_h(self,mpi_flag):
        '''
        Basic function of `atan` to .h file.
        -----------------------------------
        if c(t) = atan(a(t)) and  help function b(t)=a(t)**2:
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*(1+b[0]);
           b[t] = 2*sum{j=0->k-1} a[j]*(k-j)*a[k-j]/k;
        c: c = atan(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_atan_c(B,C,k);","\tmpfr_atan(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_atan(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_atan(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_add_si(B[0], B[0], 1, GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_atan_c(B,C,k);","\tmpfr_atan(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_atan(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_atan(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_add_si(B[0], B[0], 1, GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def atanh_h(self,mpi_flag):
        '''
        Basic function of `atanh` to .h file.
        -----------------------------------
        if c(t) = atanh(a(t)) and  help function b(t)=-a(t)**2:
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*(1+b[0]);
           b[t] = -2*sum{j=0->k-1} a[j]*(k-j)*a[k-j]/k;
        c: c = atanh(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_atanh_c(B,C,k);","\tmpfr_atanh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_atanh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_atanh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_atanh_c(B,C,k);","\tmpfr_atanh(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_atanh(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_atanh(C[0], A[0], GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","\t\tmpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def acot_h(self,mpi_flag):
        '''
        Basic function of `acot` to .h file.
        -----------------------------------
        if c(t) = acot(a(t)) and  help function b(t)=-a(t)**2:
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*(-1+b[0]);
           b[t] = -2*sum{j=0->k-1} a[j]*(k-j)*a[k-j]/k;
        c: c = acot(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_acot_c(B,C,k);","\tmpfr_acot(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acot(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_si_div(_son,1,A[0],GMP_RNDN);","\t\tmpfr_atan(C[0], _son, GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","mpfr_si_sub(B[0], -1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_acot_c(B,C,k);","\tmpfr_acot(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acot(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_si_div(_son,1,A[0],GMP_RNDN);","\t\tmpfr_atan(C[0], _son, GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","mpfr_si_sub(B[0], -1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
    
    def acoth_h(self,mpi_flag):
        '''
        Basic function of `acoth` to .h file.
        -----------------------------------
        if c(t) = acoth(a(t)) and  help function b(t)=-a(t)**2:
        s: c[k] = (k*b[k]-sum{j=1->k-1} jc[j]b[k-j]) / k*(1+b[0]);
           b[t] = -2*sum{j=0->k-1} a[j]*(k-j)*a[k-j]/k;
        c: c = acoth(a);
        '''
        if mpi_flag:
            type_dict = {'s': "\\\n".join(["#define tsm_acoth_c(B,C,k);","\tmpfr_acoth(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acoth(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_si_div(_son,1,A[0],GMP_RNDN);","\t\tmpfr_atanh(C[0], _son, GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","mpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = _id+1; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = _id; _i < k; _i+=_size){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tpack_mpf(_father,1,packed_1);","\t\tMPI_Reduce(packed_1,packed_rec1,1,MPI_MPF,MPI_MPF_SUM,0,MPI_COMM_WORLD);","\t\tMPI_Bcast(packed_rec1,1,MPI_MPF,0,MPI_COMM_WORLD);","\t\tunpack_mpf(packed_rec1,_father,1);","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        else:
            type_dict = {'s': "\\\n".join(["#define tsm_acoth_c(B,C,k);","\tmpfr_acoth(C,B,GMP_RNDN);"])+"\n",
                         'c': "\\\n".join(["#define tsm_acoth(A,B,C,k);","\tmpfr_set_str(_son, \"0.0\", 10, GMP_RNDN);","\tmpfr_set_str(_father, \"0.0\", 10,GMP_RNDN);","\tif (k == 0) {","\t\tmpfr_si_div(_son,1,A[0],GMP_RNDN);","\t\tmpfr_atanh(C[0], _son, GMP_RNDN);","\t\tmpfr_sqr(B[0], A[0], GMP_RNDN);","mpfr_si_sub(B[0], 1, B[0], GMP_RNDN);","\t}","\telse","\t{","\t\tfor (int _i = 1; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, C[_i], _i, GMP_RNDN);","\t\t\tmpfr_mul(_son, B[k-_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_sub(_father,A[k],_father,GMP_RNDN);","\t\tmpfr_div(C[k],_father,B[0],GMP_RNDN);","\t\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\t\tfor (int _i = 0; _i < k; _i++){","\t\t\tmpfr_mul_si(_son, A[k-_i], k-_i, GMP_RNDN);","\t\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\t\tmpfr_add(_father,_father,_son,GMP_RNDN);","\t\t}","\t\tmpfr_div_si(_father, _father, k, GMP_RNDN);","\t\tmpfr_mul_si(B[k],_father,-2,GMP_RNDN);","\t}"])+"\n"
                         }
            return "\n".join([type_dict[key] for key in type_dict])+"\n"
        
    def add_cc(self,node1,node2,ans):
        return f'tsm_add_cc({node1},{node2},{ans},i);\n'

    def add_cs(self,node1,node2,ans):
        return f'tsm_add_cs({node1},{node2},{ans},i);\n'

    def add_sc(self,node1,node2,ans):
        return f'tsm_add_sc({node1},{node2},{ans},i);\n'
    
    def add_ss(self,node1,node2,ans):
        return f'tsm_add_ss({node1},{node2},{ans},i);\n'
    
    def mul_cc(self,node1,node2,ans):
        return f'tsm_mul_cc({node1},{node2},{ans},i);\n'

    def mul_cs(self,node1,node2,ans):
        return f'tsm_mul_cs({node1},{node2},{ans},i);\n'

    def mul_sc(self,node1,node2,ans):
        return f'tsm_mul_sc({node1},{node2},{ans},i);\n'

    def mul_ss(self,node1,node2,ans):
        return f'tsm_mul_ss({node1},{node2},{ans},i);\n' 
    
    def pow(self,node1,node2,ans):
        return f'tsm_pow({node1},{node2},{ans},i);\n'  
    
    def pow_c(self,node1,node2,ans):
        return f'tsm_pow_c({node1},{node2},{ans},i);\n'  
    
    def ln(self,node,ans):
        return f'tsm_ln({node},{ans},i);\n'
    
    def exp(self,node,ans):
        return f'tsm_exp({node},{ans},i);\n'
    
    def sin_c(self,node,ans):
        return f'tsm_sin_c({node}, {ans}, i);\n'
    
    def cos_c(self,node,ans):
        return f'tsm_cos_c({node}, {ans}, i);\n'
    
    def sinh_c(self,node,ans):
        return f'tsm_sinh_c({node}, {ans}, i);\n'
    
    def cosh_c(self,node,ans):
        return f'tsm_cosh_c({node}, {ans}, i);\n'
    
    def equality(self,node,ans):
        return f'mpfr_set({ans}[i], {node}[i], GMP_RNDN);\n'

    
    #operations that in pair
    
    def sin_cos(self,node1,ans1,ans2):
        return f'tsm_sin_cos({node1},{ans1},{ans2},i);\n'
    
    def sinh_cosh(self,node1,ans1,ans2):
        return f'tsm_sinh_cosh({node1},{ans1},{ans2},i);\n'

    #operations that needs helping function
    
    def tan(self,node,help_node,ans):
        return f'tsm_tan({node},{help_node},{ans},i);\n'
    
    def tanh(self,node,help_node,ans):
        return f'tsm_tanh({node},{help_node},{ans},i);\n'
    
    def cot(self,node,help_node,ans):
        return f'tsm_cot({node},{help_node},{ans},i);\n'
    
    def coth(self,node,help_node,ans):
        return f'tsm_coth({node},{help_node},{ans},i);\n'
    
    def asin(self,node,help_node,ans):
        return f'tsm_asin({node},{help_node},{ans},i);\n'
    
    def asinh(self,node,help_node,ans):
        return f'tsm_asinh({node},{help_node},{ans},i);\n'
    
    def acos(self,node,help_node,ans):
        return f'tsm_acos({node},{help_node},{ans},i);\n'
    
    def acosh(self,node,help_node,ans):
        return f'tsm_acosh({node},{help_node},{ans},i);\n'
    
    def atan(self,node,help_node,ans):
        return f'tsm_atan({node},{help_node},{ans},i);\n'
    
    def atanh(self,node,help_node,ans):
        return f'tsm_atanh({node},{help_node},{ans},i);\n'
    
    def acot(self,node,help_node,ans):
        return f'tsm_acot({node},{help_node},{ans},i);\n'
    
    def acoth(self,node,help_node,ans):
        return f'tsm_acoth({node},{help_node},{ans},i);\n'