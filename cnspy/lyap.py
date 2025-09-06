import sympy as sp
import numpy as np
import re
from scipy.integrate import solve_ivp
from cnspy.printer import progress_bar


class lyapunov():
    """
    Using the method from Using the method from G. Benettin et al., Phys. Rev. A 14, pp 2338 (1976).
    
    Parameters
    ==========

    *funcs : A list of all odinary differential equations, the left hands of
    the equations must has a one order differential with time.

    *tspan : A tuple of integration time range.
    
    *dt : The time increment of lyapunov exponent series.
    
    *paras : A list of parameters and its values.
    
    *inits : A list of varables and its initial values.
    
    *transient_steps : Extra "transient" steps to evolve the system before calculating the lyapunov exponent.

    """
    def __init__(self, funcs, tspan, dt, paras, inits, transient_steps, safety_factor=1.0):
        self.__dim = len(inits)
        self.__iters = int((tspan[1]-tspan[0])/dt)
        self.__dt = dt
        self.__transient_steps = transient_steps
        self.__safety_factor = safety_factor
        t = sp.symbols(r't')
        u = [inits[i][0] for i in range(self.__dim)]
        self.u0 = [float(inits[i][1]) for i in range(self.__dim)]
        X = sp.Matrix(funcs).subs([(paras[i][0], paras[i][1]) for i in range(len(paras))])
        Y = sp.Matrix(u)
        self.__du = []
        for i in range(self.__dim):
            self.__du.append(sp.lambdify((t,u),X[i]))
        self.J = X.jacobian(Y)
        self.J_f = sp.lambdify((t,u), self.J)
        self.L = []
        self.calculate_lyapunov_exponents()
        
    def pack(self, u, U):
        return np.hstack((u, U.flatten()))

    def unpack(self, packed_uU):
        return packed_uU[:self.__dim], packed_uU[self.__dim:].reshape(self.__dim, self.__dim) 
    
    def funcs_for_lyapunov(self, t, state, p):
        u, U = self.unpack(state)
        du = [self.__du[i](t,u) for i in range(self.__dim)]
        dU = self.J_f(t,u) @ U
        return self.pack(du,dU)

    def calculate_lyapunov_exponents(self):
        u = self.u0.copy()
        U = np.eye(self.__dim)

        for _ in range(self.__iters):
            progress_bar(_, self.__iters-1, prefix='Calculating lyapunov exponents:', length=50)
            sol = solve_ivp(
                self.funcs_for_lyapunov,
                [0, self.__dt],
                self.pack(u,U),
                t_eval=[self.__dt],
                method='DOP853',
                rtol=1e-10,
                atol=1e-10,
                args=(0,),
            )
            u, U = self.unpack(sol.y.flatten())
            U, R = np.linalg.qr(U)
            self.L.append(np.log(abs(R.diagonal())) / self.__dt)

        self.L = np.average(self.L[self.__transient_steps:], axis=0)

    def roundToNearestZero(self, number, digits=2):
        n = abs(number)
        if n < 1:
            # Find the first non-zero digit.
            # We want 3 digits, starting at that location.
            s = f'{n:.99f}'
            index = re.search('[1-9]', s).start()
            new_float = float(s[:index + digits])
            if float(s[index + digits ]) > 4:
                new_float += 10 ** -(index - 1 + digits - 1)
            return new_float
        else:
            # We want 2 digits after decimal point.
            return round(n, digits-1)

    def max_lyapunov_exponent(self,significance_digit = 2):
        return self.roundToNearestZero(max(self.L)*self.__safety_factor, digits=significance_digit)