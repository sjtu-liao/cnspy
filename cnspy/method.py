import sympy as sp
import numpy as np
from .funcs import trees
from .tsmbasefunc import tsm_base_func
from .printer import Cprinter_helper

class method():
    """
    Father method class of the strategies of Clean Numerical Simulation.The methods 
    generally divided in two type, arbitrary order and fixed order. Arbitrary
    order means the function must smooth at any point. For example, if a function
    cotains abs() or tabular values, the function isn't smooth at zero point and 
    you should use the fixed order method. However, the fixed method is not 
    finished yet.
    
    Parameters
    ==========

    methodname : Name of method.
    
    order : order of this method.
    nosd : number of significant digits of this method~(10e-nosd).
    wordsize : wordsize of this method.
    """
    def __init__(self,methodname,**kwargs):
        
        self.methodname = methodname
        if kwargs.get('order') == None:
            self.order = 4
        else:
            self.order = kwargs.get('order')
        if kwargs.get('nosd') == None:
            if kwargs.get('wordsize') == None:
                self.wordsize = 64
                self.nosd = int(self.wordsize*np.log10(2))+1
            else:
                self.wordsize = kwargs.get('wordsize')
                self.nosd = int(self.wordsize*np.log10(2))+1
        else:
            self.nosd = kwargs.get('nosd')
            self.wordsize = int(self.nosd*np.log2(10))+1

class butcher_tableau():
    '''
    Coefficients of the explicit Runge-Kutta methods.
    '''
    def __init__(self, c, A, b):
        self.c = c  # List of time step coefficients
        self.A = A  # Matrix of coefficients for stages
        self.b = b  # List of coefficients for the final solution
        self.ki = len(self.c)  # Number of stages
        
class explicit_Runge_Kutta_frame():
    '''
    Frame of the explicit Runge-Kutta methods.
    '''
    def __init__(self, butcher_tableau):
        self.butcher = butcher_tableau
    
    def solve(self, f, u0, t0, h, num_steps, para):
        t_values = np.zeros(num_steps + 1)
        u_values = np.zeros((len(u0), num_steps + 1))
        
        t_values[0] = t0
        u_values[:,0] = u0

        for i in range(num_steps):
            k = np.zeros((self.butcher.ki, len(u0)))
            for j in range(self.butcher.ki):
                ti = t_values[i] + h * self.butcher.c[j]
                ui = u_values[:,i].copy()
                for l in range(j):
                    ui += h * self.butcher.A[j][l] * k[l]
                k[j] = f(ti, ui, para)

            u_new = u_values[:,i].copy()
            for j in range(self.butcher.ki):
                u_new += h * self.butcher.b[j] * k[j]

            t_values[i + 1] = t_values[i] + h
            u_values[:,i + 1] = u_new.copy()

        return t_values, u_values

class RK4(explicit_Runge_Kutta_frame):
    '''
    The classic Runge-Kutta method of order 4.
    '''
    def __init__(self,f, u0, tspan, para):
        c = [0, 1/2, 1/2, 1]
        A = [[0, 0, 0, 0],
             [1/2, 0, 0, 0],
             [0, 1/2, 0, 0],
             [0, 0, 1, 0]]
        b = [1/6, 1/3, 1/3, 1/6]
        butcher = butcher_tableau(c, A, b)
        explicit_Runge_Kutta_frame.__init__(self, butcher)
        self.t, self.y = self.solve(f, u0, tspan[0], ((tspan[1]-tspan[0])/tspan[2]), tspan[2],para)
    
class tsm(method,tsm_base_func,Cprinter_helper):
    """
    Father method class of the strategies of Clean Numerical Simulation.The methods 
    generally divided in two type, arbitrary order and fixed order. Arbitrary
    order means the function must smooth at any point. For example, if a function
    cotains abs() or tabular values, the function isn't smooth at zero point and 
    you should use the fixed order method. However, the fixed method is not 
    finished yet.
    
    Parameters
    ==========

    trees : A list of nodes of all odinary differential equations.
    
    order : order of this method.
    
    wordsize : wordsize of this method.
    """
    __tree: trees
    def __init__(self,trees,order=4,wordsize=16):
        method.__init__(self,"Taylor Series Method")
        tsm_base_func.__init__(self)
        Cprinter_helper.__init__(self)
        self.order = order
        self.wordsize = wordsize
        self.__tree = trees

    def __str__(self):
        return f'Method of {self.methodname} with order of {self.order} genrated.'
    
    def __repr__(self):
        return f'<Class: {self.methodname}; Order: {self.order}>'

    def __str_tsm_general_code_h(self,fun,mpi_flag):

        if fun == sp.core.add.Add:
            return self.add_h(mpi_flag)
        elif fun == sp.core.mul.Mul:
            return self.mul_h(mpi_flag)
        elif fun == sp.core.power.Pow:
            return self.pow_h(mpi_flag)
        elif fun == sp.sin:
            return self.sin_h(mpi_flag)
        elif fun == sp.cos:
            return self.cos_h(mpi_flag)
        elif fun == sp.sinh:
            return self.sinh_h(mpi_flag)
        elif fun == sp.cosh:
            return self.cosh_h(mpi_flag)
        elif fun == sp.ln:
            return self.ln_h(mpi_flag)
        elif fun == sp.exp:
            return self.exp_h(mpi_flag)
        #tan, tanh, cot, coth, asin, asinh, acos, acosh, atan, atanh, acot, acoth
        elif fun == sp.tan:
            return self.tan_h(mpi_flag)
        elif fun == sp.tanh:
            return self.tanh_h(mpi_flag)
        elif fun == sp.cot:
            return self.cot_h(mpi_flag)
        elif fun == sp.coth:
            return self.coth_h(mpi_flag)
        elif fun == sp.asin:
            return self.asin_h(mpi_flag)
        elif fun == sp.asinh:
            return self.asinh_h(mpi_flag)
        elif fun == sp.acos:
            return self.acos_h(mpi_flag)
        elif fun == sp.acosh:
            return self.acosh_h(mpi_flag)
        elif fun == sp.atan:
            return self.atan_h(mpi_flag)
        elif fun == sp.atanh:
            return self.atanh_h(mpi_flag)
        elif fun == sp.acot:
            return self.acot_h(mpi_flag)
        elif fun == sp.acoth:
            return self.acoth_h(mpi_flag)
        elif fun == sp.core.relational.Eq:
            return ''
            
        else:
            raise ValueError(f"<Undefined functions:{fun}>")
   
    def __str_tsm_general_code_c(self,i,level):
        if len(self.__tree.bnode[i]) == 3:
            op1 = self.__tree.bnode[i][0]
            node1 = self.__tree.bnode[i][1]
            ans1 = self.__tree.bnode[i][2]
            if not self.__tree.slist[i][1]:
                if op1 == sp.sin:
                    return self.indent(level)+self.sin_c(node1,ans1)
                elif op1 == sp.cos:
                    return self.indent(level)+self.cos_c(node1,ans1)
                elif op1 == sp.sinh:
                    return self.indent(level)+self.sinh_c(node1,ans1)
                elif op1 == sp.cosh:
                    return self.indent(level)+self.cosh_c(node1,ans1)
                elif op1 == sp.core.relational.Eq:
                    return self.indent(level)+self.equality(node1,ans1)
                else:
                    raise ValueError("Nodefine constant function of 3 inputs")
            else:
                if op1 == sp.sin:
                    op2 = self.__tree.bnode[i+1][0]
                    node2 = self.__tree.bnode[i+1][1]
                    ans2 = self.__tree.bnode[i+1][2]
                    if op2 == sp.cos and node2 == node1:
                        return self.indent(level)+self.sin_cos(node1,ans1,ans2)
                    else:
                        raise ValueError("TSM code generate failed. Maybe caused by wrong Reording(sin_cos)")
                elif op1 == sp.cos:
                    return ''
                elif op1 == sp.sinh:
                    op2 = self.__tree.bnode[i+1][0]
                    node2 = self.__tree.bnode[i+1][1]
                    ans2 = self.__tree.bnode[i+1][2]
                    if op2 == sp.cosh and node2 == node1:
                        return self.indent(level)+self.sinh_cosh(node1,ans1,ans2)
                    else:
                        raise ValueError("TSM code generate failed. Maybe caused by wrong Reording(sinh_cosh)")
                elif op1 == sp.cosh:
                    return ''
                elif op1 == sp.ln:
                    return self.indent(level)+self.ln(node1,ans1)
                elif op1 == sp.exp:
                    return self.indent(level)+self.exp(node1,ans1)
                elif op1 == sp.core.relational.Eq:
                    return self.indent(level)+self.equality(node1,ans1)
                else:
                    raise ValueError("Nodefine function of 3 inputs")
        
        elif len(self.__tree.bnode[i]) == 4:
            op = self.__tree.bnode[i][0]
            node_1 = self.__tree.bnode[i][1]
            node_2 = self.__tree.bnode[i][2]
            ans = self.__tree.bnode[i][3]
            #node_1 is constant and node_2 is series
            if not self.__tree.slist[i][1] and self.__tree.slist[i][2]:
                if op == sp.core.add.Add:
                   return self.indent(level)+self.add_cs(node_1,node_2,ans)
                elif op == sp.core.mul.Mul:
                   return self.indent(level)+self.mul_cs(node_1,node_2,ans)
                #power is special, the node_2 always not a series. so power_c means the inputs are both constant, but i need to reconsider about the non series operation, where should i put them? at the begining of the code will be a better choice.
                elif op == sp.core.power.Pow:
                   return self.indent(level)+self.pow_c(node_1,node_2,ans)
            #node_1 is series and node_2 is constant
            elif not self.__tree.slist[i][2] and self.__tree.slist[i][1]:
                if op == sp.core.add.Add:
                   return self.indent(level)+self.add_sc(node_1,node_2,ans)
                elif op == sp.core.mul.Mul:
                   return self.indent(level)+self.mul_sc(node_1,node_2,ans) 
                elif op == sp.core.power.Pow:
                   return self.indent(level)+self.pow(node_1,node_2,ans)
            # the ans is not a series, which means the inputs are both not series.
            elif not self.__tree.slist[i][3] and not self.__tree.slist[i][1] and not self.__tree.slist[i][2]: 
                if op == sp.core.add.Add:
                   return self.indent(level)+self.add_cc(node_1,node_2,ans)
                elif op == sp.core.mul.Mul:
                   return self.indent(level)+self.mul_cc(node_1,node_2,ans)
                elif op == sp.core.power.Pow:
                   return self.indent(level)+self.pow_c(node_1,node_2,ans)     
            #the last case, the inputs are both series, the operation that needs help function will be added here.
            #tan, tanh, cot, coth, asin, asinh, acos, acosh, atan, atanh, acot, acoth
            else:    
                if op == sp.core.add.Add:
                    return self.indent(level)+self.add_ss(node_1,node_2,ans)
                elif op == sp.core.mul.Mul:
                    return self.indent(level)+self.mul_ss(node_1,node_2,ans)
                elif op == sp.core.power.Pow:
                    return self.indent(level)+self.pow(node_1,node_2,ans)
                elif op == sp.tan:
                    return self.indent(level)+self.tan(node_1,node_2,ans)
                elif op == sp.tanh:
                    return self.indent(level)+self.tanh(node_1,node_2,ans)
                elif op == sp.cot:
                    return self.indent(level)+self.cot(node_1,node_2,ans)
                elif op == sp.coth:
                    return self.indent(level)+self.coth(node_1,node_2,ans)
                elif op == sp.asin:
                    return self.indent(level)+self.asin(node_1,node_2,ans)
                elif op == sp.asinh:
                    return self.indent(level)+self.asinh(node_1,node_2,ans)
                elif op == sp.acos:
                    return self.indent(level)+self.acos(node_1,node_2,ans)
                elif op == sp.acosh:
                    return self.indent(level)+self.acosh(node_1,node_2,ans)
                elif op == sp.atan:
                    return self.indent(level)+self.atan(node_1,node_2,ans)
                elif op == sp.atanh:
                    return self.indent(level)+self.atanh(node_1,node_2,ans)
                elif op == sp.acot:
                    return self.indent(level)+self.acot(node_1,node_2,ans)
                elif op == sp.acoth:
                    return self.indent(level)+self.acoth(node_1,node_2,ans)
                
        else:
            raise ValueError("Length of the bnode is wrong!")

    
    def generate_tsm_core(self,level,ccodefilepath):
        '''
        Functions that generate taylor series method's core C code.

        Parameters
        ----------

        level : code indentation level of the generated C code.
        
        ccodefilepath: filepath of the C code.
        '''
        with open(ccodefilepath,'a',encoding='utf-8')as file:
            for i in range(0,len(self.__tree.var)):
                file.write(f"{self.indent(level)}mpfr_set({self.__tree.var[i]}[0],{self.__tree.var[i]}_p,GMP_RNDN);\n".replace('\\','').replace('(t)',''))
            file.write(f"{self.indent(level)}/*------------{self.methodname}------------*/\n")
            file.write(f"{self.indent(level)}for (i = 0; i < order; i++){{\n")
            for i in range(0,len(self.__tree.bnode)):
                file.write(self.__str_tsm_general_code_c(i,level+1).replace('\\',''))
            for i in range(0,len(self.__tree.var)):
                file.write(f"{self.indent(level+1)}mpfr_div_ui({self.__tree.var[i]}[i+1],saver_{self.__tree.var_saver[i]}[i],(i+1),GMP_RNDN);\n".replace('\\','').replace('(t)',''))
            file.write(f"{self.indent(level)}}}\n")
    
    def generate_tsm_header(self,mpi_flag,adaptive_flag,headerfilepath):
        '''
        Generate a header that contains the list of the functions, the
        header is a .h file and full of marco define.
        
        Parameters
        ----------

        mpi_flag : flag of parallel calculation.
        
        headerfilepath: filepath of the head file.
        '''
        with open(headerfilepath,'a',encoding='utf-8')as file:
            if mpi_flag:
                file.write(f"#include \"gmp.h\"\n"+
                           f"#include \"mpfr.h\"\n"+
                           f"#include \"mpi.h\"\n"+
                           f"#include \"mpf2mpfr.h\"\n"+
                           f"#include \"mpi_mpf.h\"\n\n")
                file.write(f"#define tsm_init();\\\n"+
                           f"{self.indent(1)}double st, et;\\\n"+
                           f"{self.indent(1)}int _i, _j,_k;\\\n"+
                           f"{self.indent(1)}int _id, _size;\\\n"+
                           f"{self.indent(1)}MPI_Init(&argc, &argv);\\\n"+
                           f"{self.indent(1)}MPI_Comm_size(MPI_COMM_WORLD, &_size);\\\n"+
                           f"{self.indent(1)}MPI_Comm_rank(MPI_COMM_WORLD, &_id);\\\n"+
                           f"{self.indent(1)}mpfr_set_default_prec(_prec);\\\n"+
                           f"{self.indent(1)}commit_mpf(&(MPI_MPF),_prec,MPI_COMM_WORLD);\\\n"+
                           f"{self.indent(1)}create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);\\\n"+
                           f"{self.indent(1)}mpfr_t _son, _father,_zero;\\\n"+
                           f"{self.indent(1)}mpfr_inits2(_prec, _son, _father,_zero, (mpfr_ptr) 0);\\\n"+
                           f"{self.indent(1)}void *packed_1, *packed_rec1,*packed_2, *packed_rec2; \\\n"+
                           f"{self.indent(1)}packed_1 = allocbuf_mpf(_prec, 1);\\\n"+
                           f"{self.indent(1)}packed_rec1 = allocbuf_mpf(_prec, 1);\\\n"+
                           f"{self.indent(1)}packed_2 = allocbuf_mpf(_prec, 1);\\\n"+
                           f"{self.indent(1)}packed_rec2 = allocbuf_mpf(_prec, 1);\n\n")
                file.write(f"#define tsm_change_prec(PREC);\\\n"+
                           f"{self.indent(1)}mpfr_set_default_prec(PREC);\\\n"+
                           f"{self.indent(1)}commit_mpf(&(MPI_MPF), PREC, MPI_COMM_WORLD);\\\n"+
                           f"{self.indent(1)}create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_son,PREC,GMP_RNDN);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_father,PREC,GMP_RNDN);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_zero,PREC,GMP_RNDN);\\\n"+
                           f"{self.indent(1)}packed_1 = allocbuf_mpf(PREC, 1);\\\n"+
                           f"{self.indent(1)}packed_rec1 = allocbuf_mpf(PREC, 1);\\\n"+
                           f"{self.indent(1)}packed_2 = allocbuf_mpf(PREC, 1);\\\n"+
                           f"{self.indent(1)}packed_rec2 = allocbuf_mpf(PREC, 1);\n\n")
                file.write(f"#define tsm_fina();\\\n"+
                           f"{self.indent(1)}mpfr_clears(_son, _father,_zero, NULL);\\\n"+
                           f"{self.indent(1)}free_mpf_op(&(MPI_MPF_SUM));\\\n"+
                           f"{self.indent(1)}free_mpf(&(MPI_MPF));\\\n"+
                           f"{self.indent(1)}free(packed_1);\\\n"+
                           f"{self.indent(1)}free(packed_rec1);\\\n"+
                           f"{self.indent(1)}free(packed_2);\\\n"+
                           f"{self.indent(1)}free(packed_rec2);\\\n"+
                           f"{self.indent(1)}MPI_Finalize();\n\n")
                file.write(f"#define tsm_clock_start();\\\n"+
                           f"{self.indent(1)}st=MPI_Wtime();\n\n")
                file.write(f"#define tsm_clock_end();\\\n"+
                           f"{self.indent(1)}if(_id==0)\\\n"+
                           f"{self.indent(1)}{{\\\n"+
                           f"{self.indent(2)}et=MPI_Wtime();\\\n"+
                           f"{self.indent(2)}printf(\"CPU's time is %fs\", et-st);\\\n"+
                           f"{self.indent(1)}}}\n\n")
                file.write(f"#define tsm_block_start(x);\\\n"+
                           f"{self.indent(1)}if(_id==x){{\n")
                file.write(f"#define tsm_block_end();\\\n"+
                           f"{self.indent(1)}}}\n\n")
                for fun in list(set([self.__tree.bnode[i][0] for i in range(0,len(self.__tree.bnode)) ])):
                    file.write(self.__str_tsm_general_code_h(fun,mpi_flag))
            else:
                file.write(f"#include \"gmp.h\"\n"+
                           f"#include \"mpfr.h\"\n"+
                           f"#include <sys/time.h>\n\n")
                file.write(f"double timer()\n"+
                           f"{self.indent(1)}{{\n"+
                           f"{self.indent(2)}double t;\n"+
                           f"{self.indent(2)}struct timeval tv;\n"+
                           f"{self.indent(2)}gettimeofday(&tv, NULL);\n"+
                           f"{self.indent(2)}t = (double)tv.tv_sec*1.0+(double)tv.tv_usec*1e-6;\n"+
                           f"{self.indent(2)}return t;\n"+
                           f"{self.indent(1)}}}\n\n")
                file.write(f"#define tsm_init();\\\n"+
                           f"{self.indent(1)}double st, et;\\\n"+
                           f"{self.indent(1)}int _i, _j,_k;\\\n"+
                           f"{self.indent(1)}mpfr_t _son, _father,_zero;\\\n"+
                           f"{self.indent(1)}mpfr_inits2(_prec, _son, _father,_zero, (mpfr_ptr) 0);\n\n")
                file.write(f"#define tsm_change_prec(PREC);\\\n"+
                           f"{self.indent(1)}mpfr_set_default_prec(PREC);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_son,PREC,GMP_RNDN);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_father,PREC,GMP_RNDN);\\\n"+
                           f"{self.indent(1)}mpfr_prec_round(_zero,PREC,GMP_RNDN);\n\n")
                file.write(f"#define tsm_fina();\\\n"+
                           f"{self.indent(1)}mpfr_clears(_son, _father,_zero, NULL);\n\n")
                file.write(f"#define tsm_clock_start();\\\n"+
                           f"{self.indent(1)}st=timer();\n\n")
                file.write(f"#define tsm_clock_end();\\\n"+
                           f" {self.indent(1)}et=timer();\\\n"+
                           f" {self.indent(1)}printf(\"CPU's time is %fs\", et-st);\n\n")
                for fun in list(set([self.__tree.bnode[i][0] for i in range(0,len(self.__tree.bnode)) ])):
                    file.write(self.__str_tsm_general_code_h(fun,mpi_flag))

            if mpi_flag:
                if adaptive_flag:
                    file.write("\\\n".join(["#define tsm_sum2(A,b,C,k);","\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\tfor (_i = _id; _i < k+1; _i+=_size)","\t{","\t\tmpfr_set_str(_son,\"1.0\",10,GMP_RNDN);","\t\tfor(_j=1;_j<_i+1;_j++)","\t\t\tmpfr_mul(_son,_son,b,GMP_RNDN);","\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\tmpfr_add(_father, _father, _son, GMP_RNDN);","\t}","\tpack_mpf(_father, 1, packed_1);","\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\tunpack_mpf(packed_rec1, _father, 1);\\\n"])+"\tmpfr_set(C,_father,GMP_RNDN);\n")
                else:
                    file.write("\\\n".join(["#define tsm_sum(A,B,C);","\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\tfor (_i = _id; _i < _order+1; _i+=_size)","\t{","\t\tmpfr_mul(_son, A[_i], B[_i], GMP_RNDN);","\t\tmpfr_add(_father, _father, _son, GMP_RNDN);","\t}","\tpack_mpf(_father, 1, packed_1);","\tMPI_Reduce(packed_1, packed_rec1, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);","\tMPI_Bcast(packed_rec1, 1, MPI_MPF, 0, MPI_COMM_WORLD);","\tunpack_mpf(packed_rec1, _father, 1);\\\n"])+"\tmpfr_set(C,_father,GMP_RNDN);\n")
            else:
                if adaptive_flag:
                    file.write("\\\n".join(["#define tsm_sum2(A,b,C,k);","\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\tfor (_i = 0; _i < k+1; _i++)","\t{","\t\tmpfr_set_str(_son,\"1.0\",10,GMP_RNDN);","\t\tfor(_j=1;_j<_i+1;_j++)","\t\t\tmpfr_mul(_son,_son,b,GMP_RNDN);","\t\tmpfr_mul(_son, A[_i], _son, GMP_RNDN);","\t\tmpfr_add(_father, _father, _son, GMP_RNDN);","\t}\\\n"])+"\tmpfr_set(C,_father,GMP_RNDN);\n")
                else:
                    file.write("\\\n".join(["#define tsm_sum(A,B,C);","\tmpfr_set_str(_father, \"0.0\", 10, GMP_RNDN);","\tfor (_i = 0; _i < _order+1; _i++)","\t{","\t\tmpfr_mul(_son, A[_i], B[_i], GMP_RNDN);","\t\tmpfr_add(_father, _father, _son, GMP_RNDN);","\t}\\\n"])+"\tmpfr_set(C,_father,GMP_RNDN);\n")                    