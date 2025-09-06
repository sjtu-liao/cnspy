import sympy as sp
import numpy as np
from .printer import Cprinter_helper

class adaptive():
    """
    Father method class of the adaptive methods.
    """
    def __init__(self, *args, delta_T=None, remaining_T=None, magic_number=None, gamma=1.1):
        allowed_args = {'t', 'order', 'wordsize'}
        invalid_args = [arg for arg in args if arg not in allowed_args]
        if invalid_args:
            raise ValueError(f"Invalid adaptive arguments: {invalid_args}")
        # Initialize all adaptive variables with False by default
        self.__adaptive_variables = {arg: False for arg in allowed_args}
        # Set the provided adaptive variables to True
        for arg in args:
            self.__adaptive_variables[arg] = True
        self.delta_T = delta_T
        self.remaining_T = remaining_T
        self.magic_number = magic_number
        self.gamma = gamma
        

    def get_adaptive_variables(self):
        return self.__adaptive_variables
    
class adaptive_parameters():    
    def __init__(self, mpiflag:bool, ttotal, adptv:adaptive, lyap,tspan, numofprocessor):
        self.__adptv = adptv
        self.__mpiflag = mpiflag
        self.__ttotal = ttotal
        self.__delta_T = self.__adptv.delta_T
        self.__remaining_T = self.__adptv.remaining_T
        self.__magic_number = self.__adptv.magic_number
        self.__gamma = self.__adptv.gamma
        self.adptv_t_flag = self.__adptv.get_adaptive_variables()['t']
        self.adptv_order_flag = self.__adptv.get_adaptive_variables()['order']
        self.adptv_wordsize_flag = self.__adptv.get_adaptive_variables()['wordsize']
        self.__lyap = lyap
        self.max_wordsize = int(np.log2(10)*self.__gamma*self.__lyap*(tspan[1]-tspan[0])/np.log(10.0)) + 10
        self.numofprocessor = numofprocessor
        self.set_common_parameters()
        self.set_adaptive_order_parameters()
        self.set_adaptive_wordsize_parameters()
    
    def set_common_parameters(self):
        """
        Set the common parameters of the adaptive methods.
        """
        self.delta_T = (0.005*self.__ttotal) if self.__delta_T is None else self.__delta_T
        self.remaining_T = (0.05*self.__ttotal) if self.__remaining_T is None else self.__remaining_T
            
    def set_adaptive_order_parameters(self):
        """
        Set the order adaptive parameters of the adaptive methods.
        """
        # Linear formula 
        # Set magic number
        # Default values of the magic number for series method is ln(10)/2 
        # Default values of the magic number for parallel method is 1.5
        self.magic_number = [0, np.log(10)/2, 1] if not self.__mpiflag else [0, 1.5, 0] if self.__magic_number is None else self.__magic_number
        #self.magic_number = np.log(10)/2 if not self.__mpiflag else 1.5 if self.__magic_number is None else self.__magic_number
        # quadratic formula
        if self.adptv_order_flag:
            c,d,p,rho,eps,M,mp,Ns = sp.symbols(r'c,d,p,\rho_m,\epsilon,M_m,mp,N_s')
            phi = sp.symbols(r'\phi', cls=sp.Function)
            #phi = (c*((p+1)**2/mp*1/16**2*c*Ns**2+(p+1)*sp.log(mp,2)*1/16*Ns)+d*(p+1)*1/16**2*c*Ns**2)/rho/(eps/M)**(1/(p+1))
            phi = (((p+1)**2/mp+c*(p+1)*sp.log(mp,2)))/rho/(eps/M)**(1/(p+1))
            phi_d = sp.diff(phi, p)
            #phi_d = phi_d.subs(M,1).subs(eps,10**(-Ns)).subs(rho,1).subs(c,np.log(10)/2*Ns/(mp-1)/mp/sp.log(mp,2)).simplify()
            phi_d = phi_d.subs(M,1).subs(eps,10**(-Ns)).subs(rho,1).subs(c,1).simplify()
            ps = []
            Nss = []
            for i in [i for i in range(16,self.max_wordsize,self.delta_T)]:
                ii = (self.__ttotal-i)*self.__gamma*self.__lyap/np.log(10)
                try:
                    ps.append(sp.nsolve(phi_d.subs(mp,self.numofprocessor*(0.99)**self.numofprocessor).subs(Ns,ii), p, ii))
                    Nss.append(ii)
                except:
                    pass
            self.magic_number[0], self.magic_number[1], self.magic_number[2] = np.polyfit(np.array(Nss,float),np.array(ps,float),2)
            #self.magic_number[1], self.magic_number[2] = np.polyfit(np.array(Nss,float),np.array(ps,float),1)
            #if self.__lyap < 0.4:
            #    self.magic_number[0] = 0.0 
            #    self.magfic_number[1] = np.log(10)/2
            #    self.magic_number[2] = 1.0 
            #else:
            #    self.magic_number[0], self.magic_number[1], self.magic_number[2] = np.polyfit(np.array(Nss,float),np.array(ps,float),2)
            if not self.__mpiflag:
                self.magic_number[0] = 0
                self.magic_number[1] = np.log(10)/2
                self.magic_number[2] = 1.0
        
    def set_adaptive_wordsize_parameters(self):
        """
        Set the wordsize adaptive parameters of the adaptive methods.
        """
        # Set safe coefficient gamma
        # Default values of the safe coefficient gamma is 1.1
        self.gamma = 1.1 if self.__gamma is None else self.__gamma
    
    def reset_adaptive_order_parameters(self,reset_magic_number):
        self.magic_number = reset_magic_number