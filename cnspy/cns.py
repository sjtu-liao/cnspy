import os
import sys
import shutil
import mpmath
import sympy as sp
import numpy as np
from math import ceil

from .funcs import trees
from .method import method
from .method import tsm
from .printer import Foldermaker,Cprinter_helper
from .adaptive import adaptive,adaptive_parameters
from .lyap import lyapunov
from .method import RK4
from scipy.integrate import solve_ivp

class cns(Foldermaker,Cprinter_helper):
    """
    Main class of the strategies of Clean Numerical Simulation.
    
    Parameters
    ==========
    
    *name : name of the system.

    *funcs : A list of all odinary differential equations, the left hands of
    the equations must has a one order differential with time.
    
    *paras : A list of parameters and its values.
    
    *inits : A list of varables and its initial values.
    
    *tspan : A tuple of integration time range.
    
    method : A class that contain the order, number of significant digits~(10e-nosd) and wordsize of this method.
    
    adaptive : True or False. If adaptive is true, the integration will take 
    automatic steps. Default is True.
    
    *dt : Time steps of the integration. Default is 0.01.
    
    saveat : Time steps of the points that saving to answers. Default is 0.1.
    
    *mpi : True or False. If mpi is true, the integration will be parallel computing.
    Only TS() has the mpi version. Default is False.
    
    *makefile : True or False. If makefile is True, it will genrate a makefile file to
    compile the C file. Default is True.

    run_os : The operating system of the running environment. Default is the system's os.
    Or you can set it to 'nt' (Windows) or 'posix' (Unix).

    callback : A function event that can change or terminate the integration process.
    Default is none.

    Yields
    ======
    A .cns file, a .c file , a .h head file, and a makefile file.

    Examples
    ========
    
    >>> import cns
    
    >>> import sympy as sp
        
    >>> a,b,c,d,e,f,g = sp.symbols(r'a,b,c,d,e,f,g')
        
    >>> funs = [a*(b+c)+sp.sin(d), -sp.ln(f)/g, h*sp.cos(i)*(j+k)/(l+m*n)]

    >>> cns(funcs = funs,inits=[[a,1],[b,1],[c,1]],paras=[[d,1],[e,1],[f,1],[g,1]],tspan=(0,10000))
    
    """
    __funcs: list
    __inits: list
    __paras: list
    __tspan: tuple
    __method: method
    __adaptive: adaptive
    __dt: float
    __saveat: float
    __mpi: bool
    __makefile: bool
    __lyap: float
    #__callback: function
    

    
    def __init__(self,**kwargs):
        self.windowspathdelimiter = '\\'
        self.unixpathdelimiter = '/'
        self.name = kwargs.get('name')
        Foldermaker.__init__(self,os.path.join(os.getcwd(),self.name))
        Cprinter_helper.__init__(self)
        self.__inits = kwargs.get('inits')
        self.__paras = kwargs.get('paras')
        self.__tspan = kwargs.get('tspan')
        self.__method = kwargs.get('method')
        self.__prefuncs = kwargs.get('funcs')
        if kwargs.get('run_os') == None:
            self.__os = os.name
        else:
            self.__os = kwargs.get('run_os')
        self.__dim = len(self.__inits)
        self.__funcs = []
        #for save the constant
        mpmath.mp.dps = self.__method.wordsize*np.log10(2)
        #To substitute the parameters into its value
        for i in range(0,len(self.__prefuncs)):
            self.__funcs.append(self.__prefuncs[i].subs([(self.__paras[j][0],mpmath.mp.mpf(self.__paras[j][1])) for j in range(0,len(self.__paras))]).evalf(n=ceil(self.__method.wordsize*np.log10(2)),subs={sp.pi: sp.pi.evalf(n=ceil(self.__method.wordsize*np.log10(2)))}))

        if self.__method.methodname == 'TSM' or self.__method.methodname == 'TS' or self.__method.methodname == 'Tayor Series Method':
            self.__tree=trees(self.__funcs,[self.__inits[i][0] for i in range(0,len(self.__inits))],[self.__paras[i][0] for i in range(0,len(self.__paras))])
            self.__methodprinter = tsm(self.__tree,self.__method.order,self.__method.wordsize)
        else:
            pass
        if kwargs.get('numofprocessor') == None:
            self.__numofprocessor = 1
        else:
            self.__numofprocessor = kwargs.get('numofprocessor')
        if kwargs.get('lyap') == None:
            
            self.__lyap = lyapunov(self.__prefuncs,self.__tspan,0.1,self.__paras,self.__inits,int(self.__tspan[1]*0.1)).max_lyapunov_exponent()
            print(f"Lyapunov exponent is {self.__lyap}.")
            if self.__lyap < 16*np.log(10.0)/np.log2(10)/(self.__tspan[1]-self.__tspan[0])/1.1:
                print(f"Lyapunov exponent is too small, the self-adaptive method may not working.")
        else:
            self.__lyap = kwargs.get('lyap')
        if kwargs.get('adaptive') == None:
            self.__adaptive = adaptive()
        else:
            self.__adaptive = kwargs.get('adaptive')
        if kwargs.get('dt') == None:
            self.__dt = 0.01
        else:
            self.__dt = kwargs.get('dt')
        if kwargs.get('saveat') == None:
            self.__saveat = self.__dt
        else:
            self.__saveat = kwargs.get('saveat')

        if kwargs.get('mpi') == None and (self.__method.methodname == 'Taylor Series Method' or self.__method.methodname == 'TS' or self.__method.methodname == 'TSM'):
            self.__mpi = False
        elif  kwargs.get('mpi') == None and (self.__method.methodname != 'Taylor Series Method' and self.__method.methodname != 'TS' and self.__method.methodname != 'TSM'):
            raise Exception(f"{self.__method.methodname} doesn't support mpi module.")
        elif  kwargs.get('mpi') != None and (self.__method.methodname != 'Taylor Series Method' and self.__method.methodname != 'TS' and self.__method.methodname != 'TSM'):
            raise Exception(f"{self.__method.methodname} doesn't support mpi module.")
        else:
            self.__mpi = kwargs.get('mpi')
        if kwargs.get('makefile') == None:
            self.__makefile = True
        else:
            self.__makefile = kwargs.get('makefile')
        if kwargs.get('callback') == None:
            self.__callback = None
        else:
            self.__callback = kwargs.get('callback')

        self.adptv_para = adaptive_parameters(self.__mpi,self.__tspan[1],self.__adaptive,self.__lyap,self.__tspan,self.__numofprocessor)
        if kwargs.get('magic_number') == None:
            pass
        else:
            self.adptv_para.reset_adaptive_order_parameters(kwargs.get('magic_number'))
        self.__delta_T = self.adptv_para.delta_T
        self.__remaining_T = self.adptv_para.remaining_T
        self.__magic_number = self.adptv_para.magic_number
        self.__gamma = self.adptv_para.gamma
        self.__adaptive_t_flag = self.adptv_para.adptv_t_flag
        self.__adaptive_order_flag = self.adptv_para.adptv_order_flag
        self.__adaptive_wordsize_flag = self.adptv_para.adptv_wordsize_flag
        pre_nosd = int(self.__gamma*self.__lyap*(self.__tspan[1]-self.__tspan[0])/np.log(10.0))
        #pre_order = int(self.__magic_number*pre_nosd) + 10 
        pre_order = int(self.__magic_number[0]*pre_nosd*pre_nosd+self.__magic_number[1]*pre_nosd+self.__magic_number[2])
        
        if self.__method.wordsize == 64:
            self.__method.wordsize = int(pre_nosd*np.log2(10))+1
            self.__method.nosd = pre_nosd
        if self.__method.order == 4:
            self.__method.order = pre_order
        print(f"Predefine wordsize is {self.__method.wordsize}.")
        print(f"Predefine order is {self.__method.order}.")

        if kwargs.get('printprecision') == None:
            self.__printprecision = 15
        elif kwargs.get('printprecision') == 'full':
            self.__printprecision = self.__method.nosd
        else:
            self.__printprecision = kwargs.get('printprecision')
  
        if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
            pre_nosd = int(self.__gamma*self.__lyap*(self.__tspan[1]-self.__tspan[0])/np.log(10.0)) + 20
            #pre_order = int(self.__magic_number*pre_nosd) + 10
            pre_order = int(self.__magic_number[0]*pre_nosd*pre_nosd+self.__magic_number[1]*pre_nosd+self.__magic_number[2])+5
            if self.__adaptive_wordsize_flag:
            #if self.__method.wordsize < pre_nosd:
                self.__method.wordsize = int(np.log2(10)*pre_nosd)+1
            if self.__adaptive_order_flag:
            #if self.__method.order < pre_order:
                self.__method.order = pre_order
        self.__funcs = []
        #for save the constant
        mpmath.mp.dps = self.__method.wordsize*np.log10(2)
        #To substitute the parameters into its value
        for i in range(0,len(self.__prefuncs)):
            self.__funcs.append(self.__prefuncs[i].subs([(self.__paras[j][0],mpmath.mp.mpf(self.__paras[j][1])) for j in range(0,len(self.__paras))]).evalf(n=ceil(self.__method.wordsize*np.log10(2)),subs={sp.pi: sp.pi.evalf(n=ceil(self.__method.wordsize*np.log10(2)))}))

        if self.__method.methodname == 'TSM' or self.__method.methodname == 'TS' or self.__method.methodname == 'Tayor Series Method':
            self.__tree=trees(self.__funcs,[self.__inits[i][0] for i in range(0,len(self.__inits))],[self.__paras[i][0] for i in range(0,len(self.__paras))])
            self.__methodprinter = tsm(self.__tree,self.__method.order,self.__method.wordsize)
        else:
            pass

        
    def __str__(self):
        return f'File has been generated.\nFunctions are {self.__method.funcs}.\nMethod is {self.__method.methodname} with order of {self.__method.order};\n'
    
    def __repr__(self):
        return f'<State: Success! Class: {self.__method.methodname}; Order: {self.__method.order}>'

    def generate_c_code(self):
        #Copy necessary header files
        shutil.copy2(os.path.join(os.path.dirname(os.path.abspath(__file__)),'c','mpi_gmp.h'), self.ccodepath)
        shutil.copy2(os.path.join(os.path.dirname(os.path.abspath(__file__)),'c','mpi_mpf.h'), self.ccodepath)
        #print .h file
        hfilepath=os.path.join(self.ccodepath,f"{self.name}.h")
        with open(hfilepath,'w',encoding='utf-8')as file:
            file.write(self.brief_info())
            file.write(f"#include <stdio.h>\n"+
                       f"#include <string.h>\n"+
                       f"#include <stdlib.h>\n"+
                       f"#include <stdbool.h>\n"+
                       f"#include <math.h>\n\n")
        self.__methodprinter.generate_tsm_header(self.__mpi,self.__adaptive_t_flag,hfilepath)

        #print .c file
        cfilepath=os.path.join(self.ccodepath,f"{self.name}.c")
        if self.__os == 'nt':
            ansfilepath=os.path.join(self.anspath,f"result_{self.name}.txt")
        else:
            ansfilepath=os.path.join(self.anspath,f"result_{self.name}.txt").replace(self.windowspathdelimiter,self.unixpathdelimiter)
        with open(cfilepath,'w',encoding='utf-8')as file:
            file.write(self.brief_info())
            file.write(f"{self.indent(0)}#include \"{self.name}.h\"\n")
            file.write(f"{self.indent(0)}#define _prec {self.__method.wordsize}\n")
            file.write(f"{self.indent(0)}#define _order {self.__method.order}\n")
            file.write(f"{self.indent(0)}int main(int argc, char *argv[]){{\n")
            file.write(f"{self.indent(1)}tsm_init();\n")
            file.write(f"{self.indent(1)}tsm_clock_start();\n")
            file.write(f"{self.indent(1)}FILE *fp=fopen(\"{ansfilepath}\",\"w\");\n")
            file.write(f"{self.indent(1)}int i, j, order;\n")
            file.write(f"{self.indent(1)}order = _order;\n")
            # print declare variables          
            file.write(f"{self.indent(1)}mpfr_t ")
            for i in range(0,len(self.__inits)-1):
                file.write(f"{self.__inits[i][0]}[order+1], ".replace('\\','').replace('(t)',''))
            file.write(f"{self.__inits[len(self.__inits)-1][0]}[order+1];\n".replace('\\','').replace('(t)',''))
            file.write(f"{self.indent(1)}mpfr_t ")
            for i in range(0,self.__tree.maxsavernum-1):
                if i < 0:
                    break
                else:                
                    file.write(f"saver_{i}[order+1], ")
            file.write(f"saver_{self.__tree.maxsavernum-1}[order+1];\n")
            if self.__tree.constantnum != 0:
                file.write(f"{self.indent(1)}mpfr_t ")
                for i in range(0,self.__tree.constantnum-1):
                    if i < 0:
                        break
                    else:
                        file.write(f"constant_{i}, ")
                if self.__tree.constantnum == 1:
                    file.write(f"constant_{self.__tree.constantnum-1};\n")
                else:
                    file.write(f"constant_{self.__tree.constantnum-1};\n")
            if self.__adaptive_t_flag:
                file.write(f"{self.indent(1)}mpfr_t t[order+1];\n")
            else:
                file.write(f"{self.indent(1)}mpfr_t t[order+1], timeseries[order+1];\n")
            file.write(f"{self.indent(1)}mpfr_t inc_time, final_time;\n")
            file.write(f"{self.indent(1)}mpfr_t ")
            for i in range(0,len(self.__inits)-1):
                file.write(f"{self.__inits[i][0]}_p, ".replace('\\','').replace('(t)',''))
            file.write(f"{self.__inits[len(self.__inits)-1][0]}_p;\n".replace('\\','').replace('(t)',''))   
            #print initial
            file.write(f"{self.indent(1)}mpfr_inits2(_prec, inc_time, final_time, (mpfr_ptr) 0);\n")
            file.write(f"{self.indent(1)}mpfr_inits2(_prec, ")
            for i in range(0,len(self.__inits)):
                file.write(f"{self.__inits[i][0]}_p, ".replace('\\','').replace('(t)',''))
            for i in range(0,self.__tree.constantnum):
                file.write(f"constant_{i}, ")
            file.write(f"(mpfr_ptr) 0);\n")
            file.write(f"{self.indent(1)}for(i=0; i<order+1; i++){{\n")
            if self.__adaptive_t_flag:
                file.write(f"{self.indent(2)}mpfr_inits2(_prec, t[i], ")                
            else:
                file.write(f"{self.indent(2)}mpfr_inits2(_prec, t[i], timeseries[i], ")
            for i in range(0,len(self.__inits)):
                file.write(f"{self.__inits[i][0]}[i], ".replace('\\','').replace('(t)',''))
            for i in range(0,self.__tree.maxsavernum):
                file.write(f"saver_{i}[i], ")
            file.write(f"(mpfr_ptr) 0);\n{self.indent(1)}}}\n")
            #print set value
            file.write(f"{self.indent(1)}mpfr_set_str(final_time, \"{self.__tspan[1]}\", 10, GMP_RNDN);\n")
            file.write(f"{self.indent(1)}mpfr_set_str(inc_time, \"{self.__dt}\", 10, GMP_RNDN);\n")
            file.write(f"{self.indent(1)}for(i=0; i<order+1; i++){{\n")
            file.write(f"{self.indent(2)}mpfr_set_str(t[i], \"0.0\", 10, GMP_RNDN);\n{self.indent(1)}}}\n")
            
            if self.__adaptive_t_flag:
                pass           
            else:
                file.write(f"{self.indent(1)}mpfr_set_str(timeseries[0], \"1.0\", 10, GMP_RNDN);\n")
                file.write(f"{self.indent(1)}for(i=1; i<order+1; i++){{\n")
                file.write(f"{self.indent(2)}mpfr_pow_ui(timeseries[i], inc_time, i, GMP_RNDN);\n{self.indent(1)}}}\n")
            if ceil(self.__method.wordsize*np.log10(2))+10 > 640:
                sys.set_int_max_str_digits(ceil(self.__method.wordsize*np.log10(2))+10)
            
            for i in range(0,self.__tree.constantnum):
                file.write(f"{self.indent(1)}mpfr_set_str(constant_{i}, \"{sp.sympify(mpmath.mp.mpf(self.__tree.constants[i])).evalf(n=ceil(self.__method.wordsize*np.log10(2))) }\", 10, GMP_RNDN);\n")
            
            file.write(f"{self.indent(1)}mpfr_set_str(t[0], \"{self.__tspan[0]}\", 10, GMP_RNDN);\n")
            file.write(f"{self.indent(1)}mpfr_set_str(t[1], \"1.0\", 10, GMP_RNDN);\n")
            for i in range(0,len(self.__inits)):
                file.write(f"{self.indent(1)}mpfr_set_str({self.__inits[i][0]}_p, \"{self.__inits[i][1]}\", 10, GMP_RNDN);\n".replace('\\','').replace('(t)',''))
            
            #print adaptive with precsion and order initial variables
            if self.__adaptive_t_flag or self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                file.write(f"{self.indent(1)}int fprint_flag=1;\n")
                file.write(f"{self.indent(1)}mpfr_t print_time, print_time_inc, fake_time;\n")
                file.write(f"{self.indent(1)}mpfr_inits2(_prec, print_time, print_time_inc, fake_time, (mpfr_ptr) 0);\n")                
                file.write(f"{self.indent(1)}mpfr_set_str(print_time, \"{self.__saveat}\", 10, GMP_RNDN);\n")
                file.write(f"{self.indent(1)}mpfr_set_str(print_time_inc, \"{self.__saveat}\", 10, GMP_RNDN);\n")
                file.write(f"{self.indent(1)}mpfr_set_str(fake_time, \"0.0\", 10, GMP_RNDN);\n")
            else:
                if self.__saveat >= self.__dt:
                    file.write(f"{self.indent(1)}int fprint_flag=1;\n")
                    file.write(f"{self.indent(1)}mpfr_t print_time, print_time_inc, fake_time;\n")
                    file.write(f"{self.indent(1)}mpfr_inits2(_prec, print_time, print_time_inc, fake_time, (mpfr_ptr) 0);\n")                
                    file.write(f"{self.indent(1)}mpfr_set_str(print_time, \"{self.__saveat}\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(1)}mpfr_set_str(print_time_inc, \"{self.__saveat}\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(1)}mpfr_set_str(fake_time, \"0.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(1)}double print_time_d, fake_time_d;\n")
                    file.write(f"{self.indent(1)}print_time_d = mpfr_get_d(print_time,GMP_RNDN);\n")
                    file.write(f"{self.indent(1)}fake_time_d = mpfr_get_d(fake_time,GMP_RNDN);\n")
                elif self.__saveat < self.__dt:
                    raise Exception(f"Save time {self.__saveat} should be larger than or euqal to time step {self.__dt}.")
                else:
                    pass
                
            if self.__adaptive_t_flag:
                if self.__magic_number[0] == 0 and self.__magic_number[1] != 1.15:
                    file.write(f"{self.indent(1)}mpfr_t tol, max_Mm1, max_M, abs_tmp_1, abs_tmp_2, CM, VSk1, VSk2, VSY1, VSY2;\n")
                    file.write(f"{self.indent(1)}mpfr_inits2(_prec, tol, max_Mm1, max_M, abs_tmp_1, abs_tmp_2, CM, VSk1, VSk2, VSY1, VSY2, (mpfr_ptr) 0);\n")
                else:
                    file.write(f"{self.indent(1)}mpfr_t tol, max_Mm1, max_M, max_M0, one, dotseven, abs_tmp_1, abs_tmp_2, rhom_1 ,rhom_2, rhom;\n")
                    file.write(f"{self.indent(1)}mpfr_inits2(_prec, tol, max_Mm1, max_M, max_M0, one, dotseven, abs_tmp_1, abs_tmp_2, rhom_1 ,rhom_2, rhom, (mpfr_ptr) 0);\n")
                    file.write(f"{self.indent(1)}mpfr_set_str(one, \"1.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(1)}mpfr_set_str(dotseven, \"0.7\", 10, GMP_RNDN);\n")
                
            if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                #file.write(f"{self.indent(1)}double magic_number={self.__magic_number};\n")
                file.write(f"{self.indent(1)}int PRE=(int)ceil(_prec*log10(2));\n")
                file.write(f"{self.indent(1)}double delta_T={self.__delta_T}, remaining_T={self.__remaining_T}, gamma={self.__gamma}, lyap={self.__lyap};\n")
                file.write(f"{self.indent(1)}int m=0;\n")
                file.write(f"{self.indent(1)}double print_time_d, T_s, T_c;\n")
                file.write(f"{self.indent(1)}print_time_d = mpfr_get_d(print_time,GMP_RNDN);\n")
                file.write(f"{self.indent(1)}T_s = m*print_time_d;\n")
                file.write(f"{self.indent(1)}T_c = mpfr_get_d(final_time,GMP_RNDN);\n")
                
            if self.__adaptive_t_flag and not self.__adaptive_order_flag and not self.__adaptive_wordsize_flag:
                file.write(f"{self.indent(1)}mpfr_set_str(tol,\"10.0\",10,GMP_RNDN);\n")
                file.write(f"{self.indent(1)}mpfr_pow_si(tol, tol, -_prec, GMP_RNDN);\n")

            file.write(f"{self.indent(1)}while (mpfr_less_p(t[0], final_time)){{\n")
            
            #print adaptive with precsion and order
            if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                file.write(f"{self.indent(2)}if (fprint_flag && m%(int)(delta_T/print_time_d)==0 && m<=(int)((T_c-remaining_T)/print_time_d)){{\n")
                if self.__adaptive_wordsize_flag:
                    file.write(f"{self.indent(3)}PRE=(int)ceil(gamma*lyap*(T_c-T_s)/log(10.0));\n")
                if self.__adaptive_order_flag:
                    #file.write(f"{self.indent(3)}order=(int)ceil(magic_number*PRE);\n")
                    #if self.__magic_number[0] != 0:
                    #    file.write(f"{self.indent(3)}if(PRE<_size*log2(_size)/log(10)*(1/2+lyap)){{\n")
                    #    file.write(f"{self.indent(4)}order=(int)ceil(1.15*PRE+1);\n")
                    #    file.write(f"{self.indent(3)}}}\n{self.indent(3)}else{{\n")
                    file.write(f"{self.indent(4)}order=(int)ceil({self.__magic_number[0]}*PRE*PRE+{self.__magic_number[1]}*PRE+{self.__magic_number[2]});\n")
                    #if self.__magic_number[0] != 0:
                    #    file.write(f"{self.indent(3)}}}\n")
                if self.__adaptive_t_flag:
                    file.write(f"{self.indent(3)}mpfr_set_str(tol,\"10.0\",10,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_pow_si(tol, tol, -PRE, GMP_RNDN);\n")
                if self.__adaptive_wordsize_flag:
                    file.write(f"{self.indent(3)}PRE = (int)ceil(log2(10)*PRE)+1;\n")
                    file.write(f"{self.indent(3)}tsm_change_prec(PRE);\n")
                    file.write(f"{self.indent(3)}for(i=0; i<order+1; i++){{\n")
                    file.write(f"{self.indent(4)}mpfr_prec_round(t[i], PRE, GMP_RNDN);\n")
                    for i in range(0,len(self.__inits)):
                        file.write(f"{self.indent(4)}mpfr_prec_round({self.__inits[i][0]}[i], PRE, GMP_RNDN);\n".replace('\\','').replace('(t)',''))
                    for i in range(0,self.__tree.maxsavernum):
                        file.write(f"{self.indent(4)}mpfr_prec_round(saver_{i}[i], PRE, GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}}}\n")
                    if self.__adaptive_t_flag:
                        if self.__magic_number[0] == 0 and self.__magic_number[1] != 1.15:
                            for mpfr_variable in ["tol","inc_time","max_Mm1","max_M","abs_tmp_1","abs_tmp_2","CM","VSk1","VSk2","VSY1","VSY2","print_time","fake_time","final_time"]:
                                file.write(f"{self.indent(3)}mpfr_prec_round({mpfr_variable},PRE,GMP_RNDN);\n")
                        else:
                            for mpfr_variable in ["tol","inc_time","max_Mm1","max_M","max_M0","one","abs_tmp_1","abs_tmp_2","rhom_1","rhom_2","rhom","print_time","fake_time","final_time"]:
                                file.write(f"{self.indent(3)}mpfr_prec_round({mpfr_variable},PRE,GMP_RNDN);\n")
                    for i in range(0,len(self.__inits)):
                        file.write(f"{self.indent(3)}mpfr_prec_round({self.__inits[i][0]}_p,PRE,GMP_RNDN);\n".replace('\\','').replace('(t)',''))
                    for i in range(0,self.__tree.constantnum):
                        file.write(f"{self.indent(3)}mpfr_prec_round(constant_{i},PRE,GMP_RNDN);\n")
                file.write(f"{self.indent(2)}}}\n")
            
        self.__methodprinter.generate_tsm_core(2,cfilepath)
        
        with open(cfilepath,'a',encoding='utf-8')as file:
            #print adaptive with time steps
            if  self.__adaptive_t_flag:
                if self.__magic_number[0] == 0 and self.__magic_number[1] != 1.15:
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_0[order-1],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_2,x_1[order-1],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_max(max_Mm1, abs_tmp_1, abs_tmp_2, GMP_RNDN);\n")
                    for i_remain in range(2,len(self.__inits)):
                        file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_{i_remain}[order-1],GMP_RNDN);\n")
                        file.write(f"{self.indent(2)}mpfr_max(max_Mm1, max_Mm1, abs_tmp_1, GMP_RNDN);\n")

                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_0[order],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_2,x_1[order],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_max(max_M, abs_tmp_1, abs_tmp_2, GMP_RNDN);\n")
                    for i_remain in range(2,len(self.__inits)):
                        file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_{i_remain}[order],GMP_RNDN);\n")
                        file.write(f"{self.indent(2)}mpfr_max(max_M, max_M, abs_tmp_1, GMP_RNDN);\n")

                    file.write(f"{self.indent(2)}mpfr_set_str(CM, \"1.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_div_si(CM, CM, order, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_pow(VSk1, tol, CM, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_set_str(CM, \"-1.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_div_si(VSY1, CM, order-1, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_set_str(CM, \"1.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_div_si(CM, CM, order+1, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_pow(VSk2, tol, CM, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_set_str(CM, \"-1.0\", 10, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_div_si(VSY2, CM, order, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_pow(max_Mm1, max_Mm1, VSY1, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_mul(max_Mm1, VSk1, max_Mm1, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_pow(max_M, max_M, VSY2, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_mul(max_M, VSk2, max_M, GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_min(inc_time, max_Mm1, max_M, GMP_RNDN);\n")
                else:
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_0[0],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_2,x_1[0],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_max(max_M0, abs_tmp_1, abs_tmp_2, GMP_RNDN);\n")
                    for i_remain in range(2,len(self.__inits)):
                        file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_{i_remain}[0],GMP_RNDN);\n")
                        file.write(f"{self.indent(2)}mpfr_max(max_M0, max_M0, abs_tmp_1, GMP_RNDN);\n")
                        
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_0[order-1],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_2,x_1[order-1],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_max(max_Mm1, abs_tmp_1, abs_tmp_2, GMP_RNDN);\n")
                    for i_remain in range(2,len(self.__inits)):
                        file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_{i_remain}[order-1],GMP_RNDN);\n")
                        file.write(f"{self.indent(2)}mpfr_max(max_Mm1, max_Mm1, abs_tmp_1, GMP_RNDN);\n")

                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_0[order],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_2,x_1[order],GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_max(max_M, abs_tmp_1, abs_tmp_2, GMP_RNDN);\n")
                    for i_remain in range(2,len(self.__inits)):
                        file.write(f"{self.indent(2)}mpfr_abs(abs_tmp_1,x_{i_remain}[order],GMP_RNDN);\n")
                        file.write(f"{self.indent(2)}mpfr_max(max_M, max_M, abs_tmp_1, GMP_RNDN);\n")
                    
                    file.write(f"{self.indent(2)}if (mpfr_less_p(max_M0, one)){{\n")
                    file.write(f"{self.indent(3)}mpfr_log(rhom_1,max_Mm1,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_div_si(rhom_1,rhom_1,1-order,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_log(rhom_2,max_M,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_div_si(rhom_2,rhom_2,-order,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_min(rhom,rhom_1,rhom_2,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}}}\n")
                    file.write(f"{self.indent(2)}else{{\n")
                    file.write(f"{self.indent(3)}mpfr_log(abs_tmp_1,max_M0,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_log(abs_tmp_2,max_Mm1,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_sub(rhom_1,abs_tmp_1,abs_tmp_2,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_div_si(rhom_1,rhom_1,order-1,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_log(abs_tmp_1,max_M0,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_log(abs_tmp_2,max_M,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_sub(rhom_2,abs_tmp_1,abs_tmp_2,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_div_si(rhom_2,rhom_1,order,GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_min(rhom,rhom_1,rhom_2,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}}}\n")
                    # file.write(f"{self.indent(2)}mpfr_log(inc_time,tol,GMP_RNDN);\n")
                    # file.write(f"{self.indent(2)}mpfr_div_ui(inc_time,inc_time,(order+1),GMP_RNDN);\n")
                    # file.write(f"{self.indent(2)}mpfr_add(inc_time,inc_time,rhom,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_div_ui(inc_time,dotseven,1-order,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_sub_ui(inc_time,inc_time,2,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_add(inc_time,rhom,inc_time,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}mpfr_exp(inc_time,inc_time,GMP_RNDN);\n")

                #if print time is not the same as the current time
                #then the time step will be force to be print time minus the previous time
                file.write(f"{self.indent(2)}fprint_flag = 0;\n")
                file.write(f"{self.indent(2)}mpfr_add(fake_time,t[0],inc_time,GMP_RNDN);\n")
                file.write(f"{self.indent(2)}if (mpfr_greater_p(fake_time, print_time)){{\n")
                file.write(f"{self.indent(3)}mpfr_sub(inc_time, print_time, t[0], GMP_RNDN);\n")
                file.write(f"{self.indent(3)}mpfr_add(print_time, print_time, print_time_inc, GMP_RNDN);\n")
                file.write(f"{self.indent(3)}fprint_flag = 1;\n")
                if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                    file.write(f"{self.indent(3)}m = m + 1;\n")
                    file.write(f"{self.indent(3)}T_s = m*print_time_d;\n")
                file.write(f"{self.indent(2)}}}\n")
            else:
                if self.__dt == self.__saveat:
                    if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                        file.write(f"{self.indent(2)}m = m + 1;\n")
                        file.write(f"{self.indent(2)}T_s = m*print_time_d;\n") 
                elif self.__saveat/self.__dt - int(self.__saveat/self.__dt) == 0.0:
                    file.write(f"{self.indent(2)}fprint_flag = 0;\n")
                    file.write(f"{self.indent(2)}mpfr_add(fake_time,t[0],inc_time,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}print_time_d = mpfr_get_d(print_time,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}fake_time_d = mpfr_get_d(fake_time,GMP_RNDN);\n")
                    #file.write(f"{self.indent(2)}if (mpfr_greater_p(fake_time, print_time)){{\n")
                    #using mpfr to compare the time will introduce the round-off error, so we use the double to compare the time
                    file.write(f"{self.indent(2)}if (fake_time_d >= print_time_d){{\n")
                    file.write(f"{self.indent(3)}mpfr_add(print_time, print_time, print_time_inc, GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}fprint_flag = 1;\n")                                     
                    if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                        file.write(f"{self.indent(3)}m = m + 1;\n")
                        file.write(f"{self.indent(3)}T_s = m*print_time_d;\n")
                    file.write(f"{self.indent(2)}}}\n")     
                else:
                    file.write(f"{self.indent(2)}fprint_flag = 0;\n")
                    file.write(f"{self.indent(2)}mpfr_add(fake_time,t[0],inc_time,GMP_RNDN);\n")
                    file.write(f"{self.indent(2)}if (mpfr_greater_p(fake_time, print_time)){{\n")
                    file.write(f"{self.indent(3)}mpfr_sub(inc_time, print_time, t[0], GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}mpfr_add(print_time, print_time, print_time_inc, GMP_RNDN);\n")
                    file.write(f"{self.indent(3)}fprint_flag = 1;\n")                                     
                    if self.__adaptive_order_flag or self.__adaptive_wordsize_flag:
                        file.write(f"{self.indent(3)}m = m + 1;\n")
                        file.write(f"{self.indent(3)}T_s = m*print_time_d;\n")
                    file.write(f"{self.indent(2)}}}\n")                
            #summate the taylor series coefficients
            for i in range(0,len(self.__inits)):
                if self.__adaptive_t_flag:
                    file.write(f"{self.indent(2)}tsm_sum2({self.__inits[i][0]}, inc_time, {self.__inits[i][0]}_p, order);\n".replace('\\','').replace('(t)',''))
                else:
                    file.write(f"{self.indent(2)}tsm_sum({self.__inits[i][0]}, timeseries, {self.__inits[i][0]}_p);\n".replace('\\','').replace('(t)',''))
            file.write(f"{self.indent(2)}mpfr_add(t[0],t[0],inc_time,GMP_RNDN);\n")
            if self.__adaptive_t_flag or self.__adaptive_order_flag or self.__adaptive_wordsize_flag or self.__saveat > self.__dt:
                file.write(f"{self.indent(2)}if (fprint_flag){{\n")
            if self.__mpi:
                file.write(f"{self.indent(2)}tsm_block_start(0);\n")
                file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%.2Rf\\t\",t[0]);\n")
                for i in range(0,len(self.__inits)-1):
                    file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%1.{self.__printprecision}Re\\t\",{self.__inits[i][0]}_p);\n".replace('(t)',''))
                file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%1.{self.__printprecision}Re\\n\",{self.__inits[len(self.__inits)-1][0]}_p);\n".replace('(t)',''))
                #file.write(f"{self.indent(1)}tsm_clock_end();\n")
                file.write(f"{self.indent(2)}tsm_block_end();\n")
            else:
                file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%.2Rf\\t\",t[0]);\n")
                for i in range(0,len(self.__inits)-1):
                    file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%1.{self.__printprecision}Re\\t\",{self.__inits[i][0]}_p);\n".replace('(t)',''))
                file.write(f"{self.indent(2)}mpfr_fprintf(fp,\"%1.{self.__printprecision}Re\\n\",{self.__inits[len(self.__inits)-1][0]}_p);\n".replace('(t)',''))
            if self.__adaptive_t_flag or self.__adaptive_order_flag or self.__adaptive_wordsize_flag or self.__saveat > self.__dt:
                file.write(f"{self.indent(2)}}}\n")
            file.write(f"{self.indent(1)}}}\n")
            file.write(f"{self.indent(1)}tsm_clock_end();\n")
            file.write(f"{self.indent(1)}tsm_fina();\n")
            file.write(f"{self.indent(1)}return 0;\n")
            file.write(f"{self.indent(0)}}}\n")
        
        #print makefile
        if self.__makefile:
            makefilepath=os.path.join(self.ccodepath,"makefile")
            with open(makefilepath,'w',encoding='utf-8')as file:
                file.write(f"{self.name}:{self.name}.c\n")
                if self.__mpi:    
                    file.write(f"\tmpicc {self.name}.c -o {self.name} -lmpfr -lgmp -lm")
                else:
                    file.write(f"\tgcc {self.name}.c -o {self.name}  -lmpfr -lgmp -lm")

    def quick_run(self,init=None,tspan=None,num_steps=None,dt=None,using_scipy=False, scipy_method='RK45', scipy_tol=1e-6):
        '''
        Run the generated function using double precision and RK4 method.
        '''
        qr_tspan = self.__tspan if tspan is None else tspan
        qr_dt = self.__dt if dt is None else dt
        qr_num_steps = int((qr_tspan[1]-qr_tspan[0])/qr_dt) if num_steps is None else num_steps
        t = sp.symbols(r't')
        u = [self.__inits[i][0] for i in range(self.__dim)]
        X = sp.Matrix(self.__prefuncs).subs([(self.__paras[i][0], self.__paras[i][1]) for i in range(len(self.__paras))])
        du = []
        for i in range(self.__dim):
            du.append(sp.lambdify((t,u),X[i]))
        def numerical_funcs(t, state, p):
            return [du[i](t,state) for i in range(self.__dim)]
        if self.__saveat > self.__dt:
            using_scipy = True
            qr_dt = self.__saveat
        if using_scipy:
            self.qrsol = solve_ivp(numerical_funcs,[qr_tspan[0], qr_tspan[1]],[float(self.__inits[i][1]) for i in range(self.__dim)] if init is None else init,args=(0,),method=scipy_method, t_eval=np.arange(qr_tspan[0], qr_tspan[1]+qr_dt, qr_dt),rtol=scipy_tol,atol=scipy_tol)
        else:
            self.qrsol = RK4(numerical_funcs, [float(self.__inits[i][1]) for i in range(self.__dim)] if init is None else init, [qr_tspan[0],qr_tspan[1],qr_num_steps], (0,))
    
    def run(self,ssh=None,SL=None):
        if ssh is None:
            os.system(f"cd {self.ccodepath} && make")
            if np == 1:
                os.system(f"cd {self.ccodepath} && ./{self.name}")
            else:
                os.system(f"cd {self.ccodepath} && mpirun -n {np} ./{self.name}")
        else:
            try:
                ssh.connect()
                SL.clear_script(self.ccodepath)
                SL.write_script(self.ccodepath,mlx4_0=False)
                target_parrentpath = os.path.join(ssh.working_directory(),self.name).replace(self.windowspathdelimiter,self.unixpathdelimiter)
                target_ccodepath = os.path.join(target_parrentpath,'ccode').replace(self.windowspathdelimiter,self.unixpathdelimiter)
                ssh.execute(f'mkdir -p {target_parrentpath.replace(self.windowspathdelimiter,self.unixpathdelimiter)}')
                ssh.upload_directory(self.parrentpath.replace(self.windowspathdelimiter,self.unixpathdelimiter),target_parrentpath.replace(self.windowspathdelimiter,self.unixpathdelimiter))
                ssh.execute(f'cd {target_ccodepath.replace(self.windowspathdelimiter,self.unixpathdelimiter)}; make; {SL.submit_command()}')
            except:
                raise Exception("SSH connection failed, please check the ssh object.")
    class CNSsol():
        def __init__(self, name, inits, dim):
            self.name = name
            self.inits = inits
            self.dim = dim
            self.y = []
            self.t = []  
            self.t.append(0)
            self.y.append([mpmath.mp.mpf(self.inits[i][1]) for i in range(self.dim)])
            self.loadtxt_mp(os.path.join(os.path.join(self.name,"ans"),f"result_{self.name}.txt")) 
        def loadtxt_mp(self, filename, delimiter="\t"):
            """Reads a file and splits each line into a list based on the delimiter."""

            with open(filename, 'r') as file:

                for line in file:
                    line = line.strip()  # Remove leading/trailing whitespace
                    row = line.split(delimiter)
                    row = list(map(mpmath.mp.mpf,row))
                    self.t.append(row[0])
                    self.y.append(row[1:])
            self.y = np.array(self.y).transpose()
            self.t = np.array(self.t).astype(float)

    def read_cnssol(self):
        self.cnssol = self.CNSsol(self.name, self.__inits, self.__dim)
                

                    
    def analysis(self,filepath):
        pass

    @property
    def tspan(self):
        return self.__tspan
    @property
    def tree(self):
        return self.__tree
    @property
    def mpi(self):
        return self.__mpi
    @property
    def lyap(self):
        return self.__lyap
    @property
    def funcs(self):
        return self.__funcs
    @property
    def adptv(self):
        return self.adptv_para
    @property
    def numofprocessor(self):
        return self.__numofprocessor
    @property
    def printprecision(self):
        return self.__printprecision    
    
   
    

        
