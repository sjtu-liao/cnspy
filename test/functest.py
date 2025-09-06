import shutil, os
import subprocess
import sympy as sp
from cnspy.tsmbasefunc import tsm_base_func
from cnspy.printer import Cprinter_helper
import unittest

t = sp.symbols(r't')
def taylorcoeff(f,t0,n):
    coeff = []
    coeff.append(f.subs(t,t0))
    for i in range(1,n):
        coeff.append(sp.diff(f,t,i).subs(t,t0) / sp.factorial(i))
    return coeff

def str_tsm_general_code_h(fun,mpi_flag=False):
    tbf = tsm_base_func()
    if fun == sp.core.add.Add:
        return tbf.add_h(mpi_flag)
    elif fun == sp.core.mul.Mul:
        return tbf.mul_h(mpi_flag)
    elif fun == sp.core.power.Pow:
        return tbf.pow_h(mpi_flag)
    elif fun == sp.sin:
        return tbf.sin_h(mpi_flag)
    elif fun == sp.cos:
        return tbf.cos_h(mpi_flag)
    elif fun == sp.sinh:
        return tbf.sinh_h(mpi_flag)
    elif fun == sp.cosh:
        return tbf.cosh_h(mpi_flag)
    elif fun == sp.ln:
        return tbf.ln_h(mpi_flag)
    elif fun == sp.exp:
        return tbf.exp_h(mpi_flag)
    #tan, tanh, cot, coth, asin, asinh, acos, acosh, atan, atanh, acot, acoth
    elif fun == sp.tan:
        return tbf.tan_h(mpi_flag)
    elif fun == sp.tanh:
        return tbf.tanh_h(mpi_flag)
    elif fun == sp.cot:
        return tbf.cot_h(mpi_flag)
    elif fun == sp.coth:
        return tbf.coth_h(mpi_flag)
    elif fun == sp.asin:
        return tbf.asin_h(mpi_flag)
    elif fun == sp.asinh:
        return tbf.asinh_h(mpi_flag)
    elif fun == sp.acos:
        return tbf.acos_h(mpi_flag)
    elif fun == sp.acosh:
        return tbf.acosh_h(mpi_flag)
    elif fun == sp.atan:
        return tbf.atan_h(mpi_flag)
    elif fun == sp.atanh:
        return tbf.atanh_h(mpi_flag)
    elif fun == sp.acot:
        return tbf.acot_h(mpi_flag)
    elif fun == sp.acoth:
        return tbf.acoth_h(mpi_flag)

def str_tsm_function_code(fun):
    tbf = tsm_base_func()
    if fun == sp.ln:
        return "tsm_ln(b,a,i);"
    elif fun == sp.exp:
        return "tsm_exp(b,a,i);"
    #tan, tanh, cot, coth, asin, asinh, acos, acosh, atan, atanh, acot, acoth
    elif fun == sp.tan:
        return "tsm_tan(b,h,a,i);"
    elif fun == sp.tanh:
        return "tsm_tanh(b,h,a,i);"
    elif fun == sp.cot:
        return "tsm_cot(b,h,a,i);"
    elif fun == sp.coth:
        return "tsm_coth(b,h,a,i);"
    elif fun == sp.asin:
        return "tsm_asin(b,h,a,i);"
    elif fun == sp.asinh:
        return "tsm_asinh(b,h,a,i);"
    elif fun == sp.acos:
        return "tsm_acos(b,h,a,i);"
    elif fun == sp.acosh:
        return "tsm_acosh(b,h,a,i);"
    elif fun == sp.atan:
        return "tsm_atan(b,h,a,i);"
    elif fun == sp.atanh:
        return "tsm_atanh(b,h,a,i);"
    elif fun == sp.acot:
        return "tsm_acot(b,h,a,i);"
    elif fun == sp.acoth:
        return "tsm_acoth(b,h,a,i);"

def cnspy_basicfunction(f,n,mpi_flag=False):
    cph = Cprinter_helper()
    #Copy necessary header files
    destination_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_"+str(f),'')
    os.makedirs(os.path.dirname(destination_path), exist_ok=True)
    shutil.copy2(os.path.join(os.getcwd(),'cnspy','c','mpi_gmp.h'), destination_path)
    shutil.copy2(os.path.join(os.getcwd(),'cnspy','c','mpi_mpf.h'), destination_path)
    #print .h file
    cfilepath=os.path.join(destination_path,str(f)+'.c')
    ofilepath=os.path.join(destination_path,"test_"+str(f))
    with open(cfilepath,'w',encoding='utf-8')as file:
        file.write(f"#include <stdio.h>\n"+
                    f"#include <string.h>\n"+
                    f"#include <stdlib.h>\n"+
                    f"#include <stdbool.h>\n"+
                    f"#include <math.h>\n\n")
        if mpi_flag:
            file.write(f"#include \"gmp.h\"\n"+
                        f"#include \"mpfr.h\"\n"+
                        f"#include \"mpi.h\"\n"+
                        f"#include \"mpf2mpfr.h\"\n"+
                        f"#include \"mpi_mpf.h\"\n\n")
            file.write(f"#define tsm_init();\\\n"+
                        f"{cph.indent(1)}double st, et;\\\n"+
                        f"{cph.indent(1)}int _i, _j,_k;\\\n"+
                        f"{cph.indent(1)}int _id, _size;\\\n"+
                        f"{cph.indent(1)}MPI_Init(&argc, &argv);\\\n"+
                        f"{cph.indent(1)}MPI_Comm_size(MPI_COMM_WORLD, &_size);\\\n"+
                        f"{cph.indent(1)}MPI_Comm_rank(MPI_COMM_WORLD, &_id);\\\n"+
                        f"{cph.indent(1)}mpfr_set_default_prec(_prec);\\\n"+
                        f"{cph.indent(1)}commit_mpf(&(MPI_MPF),_prec,MPI_COMM_WORLD);\\\n"+
                        f"{cph.indent(1)}create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);\\\n"+
                        f"{cph.indent(1)}mpfr_t _son, _father,_zero;\\\n"+
                        f"{cph.indent(1)}mpfr_inits2(_prec, _son, _father,_zero, (mpfr_ptr) 0);\\\n"+
                        f"{cph.indent(1)}void *packed_1, *packed_rec1,*packed_2, *packed_rec2; \\\n"+
                        f"{cph.indent(1)}packed_1 = allocbuf_mpf(_prec, 1);\\\n"+
                        f"{cph.indent(1)}packed_rec1 = allocbuf_mpf(_prec, 1);\\\n"+
                        f"{cph.indent(1)}packed_2 = allocbuf_mpf(_prec, 1);\\\n"+
                        f"{cph.indent(1)}packed_rec2 = allocbuf_mpf(_prec, 1);\n\n")
            file.write(f"#define tsm_fina();\\\n"+
                        f"{cph.indent(1)}mpfr_clears(_son, _father,_zero, NULL);\\\n"+
                        f"{cph.indent(1)}free_mpf_op(&(MPI_MPF_SUM));\\\n"+
                        f"{cph.indent(1)}free_mpf(&(MPI_MPF));\\\n"+
                        f"{cph.indent(1)}free(packed_1);\\\n"+
                        f"{cph.indent(1)}free(packed_rec1);\\\n"+
                        f"{cph.indent(1)}free(packed_2);\\\n"+
                        f"{cph.indent(1)}free(packed_rec2);\\\n"+
                        f"{cph.indent(1)}MPI_Finalize();\n\n")
            file.write(f"#define tsm_clock_start();\\\n"+
                        f"{cph.indent(1)}st=MPI_Wtime();\n\n")
            file.write(f"#define tsm_clock_end();\\\n"+
                        f"{cph.indent(1)}if(_id==0)\\\n"+
                        f"{cph.indent(1)}{{\\\n"+
                        f"{cph.indent(2)}et=MPI_Wtime();\\\n"+
                        f"{cph.indent(2)}printf(\"CPU's time is %fs\", et-st);\\\n"+
                        f"{cph.indent(1)}}}\n\n")
            file.write(f"#define tsm_block_start(x);\\\n"+
                        f"{cph.indent(1)}if(_id==x){{\n")
            file.write(f"#define tsm_block_end();\\\n"+
                        f"{cph.indent(1)}}}\n\n")
        else:
            file.write(f"#include \"gmp.h\"\n"+
                        f"#include \"mpfr.h\"\n"+
                        f"#include <sys/time.h>\n\n")
            file.write(f"double timer()\n"+
                        f"{cph.indent(1)}{{\n"+
                        f"{cph.indent(2)}double t;\n"+
                        f"{cph.indent(2)}struct timeval tv;\n"+
                        f"{cph.indent(2)}gettimeofday(&tv, NULL);\n"+
                        f"{cph.indent(2)}t = (double)tv.tv_sec*1.0+(double)tv.tv_usec*1e-6;\n"+
                        f"{cph.indent(2)}return t;\n"+
                        f"{cph.indent(1)}}}\n\n")
            file.write(f"#define tsm_init();\\\n"+
                        f"{cph.indent(1)}double st, et;\\\n"+
                        f"{cph.indent(1)}int _i, _j,_k;\\\n"+
                        f"{cph.indent(1)}mpfr_t _son, _father,_zero;\\\n"+
                        f"{cph.indent(1)}mpfr_inits2(_prec, _son, _father,_zero, (mpfr_ptr) 0);\n\n")
            file.write(f"#define tsm_fina();\\\n"+
                        f"{cph.indent(1)}mpfr_clears(_son, _father,_zero, NULL);\n\n")
            file.write(f"#define tsm_clock_start();\\\n"+
                        f"{cph.indent(1)}st=timer();\n\n")
            file.write(f"#define tsm_clock_end();\\\n"+
                        f" {cph.indent(1)}et=timer();\\\n"+
                        f" {cph.indent(1)}printf(\"CPU's time is %fs\", et-st);\n\n")
        file.write(str_tsm_general_code_h(f,mpi_flag))
        #file.write(f"#include \"../func.c\"")
        file.write(f"\n#define _prec 100\n")
        file.write(f"#define order {n-1}\n")
        file.write("int main(int argc, char *argv[]){\ntsm_init();\nint i,j;\nmpfr_t t;\nmpfr_t b[order+1],a[order+1],h[order+1];\n")
        file.write(f"{cph.indent(1)}for (int i=0;i<=order;i++){{\n{cph.indent(2)}mpfr_inits2(_prec, a[i], b[i], h[i], (mpfr_ptr) 0);\n{cph.indent(2)}}}\n{cph.indent(1)}mpfr_set_str(t, \"0.0\", 10, MPFR_RNDN);")
        if f == sp.asin or f == sp.acos or f==sp.atanh:
            file.write(f"{cph.indent(1)}for (i=0;i<=order;i++){{\n"+
                        f"{cph.indent(2)}mpfr_set_str(a[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(b[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(h[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}}}\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[1], \"1.0\", 10, MPFR_RNDN);\n")
        elif f == sp.ln:
            file.write(f"{cph.indent(1)}for (i=0;i<=order;i++){{\n"+
                        f"{cph.indent(2)}mpfr_set_str(a[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(b[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(h[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}}}\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[0], \"1.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[1], \"1.0\", 10, MPFR_RNDN);\n")
        elif f == sp.acosh:
            file.write(f"{cph.indent(1)}for (i=0;i<=order;i++){{\n"+
                        f"{cph.indent(2)}mpfr_set_str(a[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(b[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(h[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}}}\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[0], \"2.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[1], \"1.0\", 10, MPFR_RNDN);\n")
        elif f==sp.acoth:
            file.write(f"{cph.indent(1)}for (i=0;i<=order;i++){{\n"+
                        f"{cph.indent(2)}mpfr_set_str(a[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(b[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}mpfr_set_str(h[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(2)}}}\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[0], \"3.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[1], \"1.0\", 10, MPFR_RNDN);\n")            
        else:
            file.write(f"{cph.indent(1)}for (i=0;i<=order;i++){{\n"+
                        f"{cph.indent(1)}mpfr_set_str(a[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}mpfr_set_str(b[i], \"1.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}mpfr_set_str(h[i], \"0.0\", 10, MPFR_RNDN);\n"+
                        f"{cph.indent(1)}}}\n"+
                        f"{cph.indent(1)}for (i=1;i<=order;i+=1){{\n"+
                        f"{cph.indent(2)}for (j=i;j>=1;j-=1){{\n"+
                        f"{cph.indent(3)}mpfr_div_ui(b[i],b[i],j,MPFR_RNDN);\n"+
                        f"{cph.indent(2)}}}\n"+
                        f"{cph.indent(1)}}}\n")
        file.write(f"{cph.indent(1)}for (i=0;i<=order;i++) {{"+str_tsm_function_code(f)+f"}}\n")
        file.write(f"{cph.indent(1)}for (i=0;i<=order;i++) mpfr_printf(\"%1.15Rf\t\",a[i]);\n")
        file.write(f"{cph.indent(1)}tsm_fina();\n{cph.indent(1)}return 0;\n}}")

    if mpi_flag:    
        return "/usr/local/Cellar/open-mpi/4.1.1_2/bin/mpicc "+cfilepath+" -o "+ofilepath+" -I/usr/local/Cellar/gmp/6.2.0/include -I/usr/local/Cellar/mpfr/4.1.0/include -lm -L/usr/local/Cellar/mpfr/4.1.0/lib -lmpfr -L/usr/local/Cellar/gmp/6.2.0/lib -lgmp -lm"
    else:
        return "gcc "+cfilepath+" -o "+ofilepath+" -I/usr/local/Cellar/gmp/6.2.0/include -I/usr/local/Cellar/mpfr/4.1.0/include -lm -L/usr/local/Cellar/mpfr/4.1.0/lib -lmpfr -L/usr/local/Cellar/gmp/6.2.0/lib -lgmp -lm"
 
def base_test_frame(test_func,real_func,order):
    destination_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_"+str(test_func),'')
    ofilepath=os.path.join(destination_path,"test_"+str(test_func))
    tc = taylorcoeff(real_func,0,order)
    sympy_result = [ tc[i].evalf() for i in range(order)]
    #compile_command = cnspy_basicfunction(test_func,order,mpi_flag=False)
    compile_command = cnspy_basicfunction(test_func,order,mpi_flag=True)
    subprocess.run(compile_command, shell=True, check=True)
    #cnspy_result_string = subprocess.check_output("./"+"test_"+str(test_func), cwd=destination_path).decode('utf-8').split()
    #print(subprocess.check_output(["mpirun", "--version"]))
    cnspy_result_string = subprocess.check_output(["mpirun","-n","4",destination_path+"test_"+str(test_func)]).decode('utf-8').split()
    cnspy_result = [float(num) for num in cnspy_result_string]
    return sympy_result, cnspy_result

class TestAddNumbers(unittest.TestCase):

    def test_ln(self):
        f = sp.ln
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(t+1),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)

    def test_exp(self):
        f = sp.exp
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15) 
                   
    def test_tan(self):
        f = sp.tan
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)   

    def test_tanh(self):
        f = sp.tanh
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)   

    def test_cot(self):
        f = sp.cot
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_coth(self):
        f = sp.coth
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_asin(self):
        f = sp.asin
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(t),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_asinh(self):
        f = sp.asinh
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_acos(self):
        f = sp.acos
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(t),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_acosh(self):
        f = sp.acosh
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(t+2),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)

    def test_atan(self):
        f = sp.atan
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_atanh(self):
        f = sp.atanh
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(t),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_acot(self):
        f = sp.acot
        order = 21
        sympy_result, cnspy_result = base_test_frame(f,f(sp.exp(t)),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)  

    def test_acoth(self):
        f = sp.acoth
        order = 15
        sympy_result, cnspy_result = base_test_frame(f,f(t+3),order)
        for i in range(order):
            self.assertAlmostEqual(cnspy_result[i], sympy_result[i], delta=1e-15)

if __name__ == '__main__':
    unittest.main()