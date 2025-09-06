
Theory
======

Basic ideas of the CNS
----------------------

Clean Numerical Simulation (CNS) :cite:`liao2008reliability` is a set of strategies with a purpose to control "numerical noises" arbitrarily. Its main point is decreasing the truncation error and the round-off error to a required level.

Truncation errors come from the discretization of continuous systems. Numerical methods have the following general form:

.. math::

   f(t+h) = f(t) + h \cdot RHS(t)

where :math:`RHS(t)` is the right-hand side, it varies according to the numerical methods. For the N-th order Taylor series method, the right-hand side can written as:

.. math::

    RHS(t) = \sum\limits_{i=1}^{N}\dfrac{f^{(i)}(t)}{i!}h^{i-1} + \mathcal{O}(h^{N})

where :math:`\mathcal{O}(h^{N})` is the order of the global truncation errors.

The CNS suggests that reducing the timestep :math:`h` or increasing :math:`N` with high order methods can decrease truncation errors to a small level which will not damage the long-term prediction. With the help of extended precision formats to store extra digits in computers, round-off errors can be largely avoided. Consequently, the numerical noises of any simulation can be controlled arbitrarily small. 

Based on the ability to control numerical noises arbitrarily small, the CNS can be used to estimate the critical predictable time :math:`T_c` of a chaotic system. The critical predictable time :math:`T_c` is defined as the time when the averaged level of numerical noise reaches a certain level :math:`\epsilon_c`.

The averaged level of numerical noise should increase exponentially within a temporal interval [0, :math:`T_c`], i.e.

.. math::

    \varepsilon(t) = \varepsilon_0 \exp(\kappa t), \quad 0 \leq t \leq T_c

So the critical predictable time :math:`T_c` can be estimated by:

.. math::

    T_c = \dfrac{1}{k}\ln(\dfrac{\varepsilon_c}{\varepsilon_0})

For more information, please refer to Liao's review article :cite:`liao2017clean`.

The old way
------------

In the past studies, using CNS to solve the ordinary differential equations (ODEs) is to differentiate both sides of the equations to get the recursion formulas of the variables. For example, for the Lorenz equations:

.. math::

   \dot{x}(t) & = \sigma (y(t) - x(t))

   \dot{y}(t) & = R x(t) - y(t) -x(t)z(t)

   \dot{z}(t) & = x(t)y(t) - b z(t)

where :math:`\sigma, R, b` are physical parameters,  the dot denotes the differentiation with respect to the time :math:`t`, respectively.

Then, differentiating both sides of the Lorenz eqs  :math:`(N-1)` times with respect to :math:`t` and then dividing them by :math:`n!`, we obtain the high-order temporal derivatives:

.. math::

    x^{[n+1]} & = \dfrac{\sigma}{n+1}(y^{[n]}-x^{[n]})

    y^{[n+1]} & = \dfrac{1}{n+1}(Rx^{[n]}-y^{[n]}-\sum\limits_{j=0}^{n}x^{[j]}z^{[n-j]})

    z^{[n+1]} & = \dfrac{1}{n+1}(bz^{[n]}+\sum\limits_{j=0}^{n}x^{[j]}y^{[n-j]})

where :math:`x^{[n]}, y^{[n]}, z^{[n]}` are the coefficients of the Taylor series:

.. math::

    x^{[n]} & = \dfrac{1}{n!}\dfrac{\partial^n x(t)}{\partial t^n}

    y^{[n]} & = \dfrac{1}{n!}\dfrac{\partial^n y(t)}{\partial t^n} 

    z^{[n]} & = \dfrac{1}{n!}\dfrac{\partial^n z(t)}{\partial t^n} 

To achieve iterative solutions for simple equations, manual methods may suffice. However, for more complex equations, manually solving the iterative process can be tedious. In such cases, it becomes necessary to adopt automated methods for derivative finding.

A new way
----------

Automatic Differentiation (AD) is a set of techniques to evaluate the derivate of a function specified by a computer program. In the method of Taylor series method (TSM), it defines basic arithmetic of Taylor series, and use these basic arithmetic of Taylor series, one can automatical get the arbitrary derivative of the function.

The basic arithmetic of Taylor series are derivated by Leibniz rule:

for the :math:`a(t) = b(t)c(t)`, it gives the following formula:


.. math::

    a^{[n]}(t) = \sum\limits_{j=0}^{n}b^{[j]}(t)c^{[n-j]}(t) 


where  :math:`a^{[n]}, b^{[n]}, c^{[n]}` are the coefficients of the Taylor series defined in the last section.

Besides, the derivative of the :math:`a(t), b(t), c(t)` can be write as:

.. math::

    a^{'}(t) = \sum\limits_{j=0}^{n-1}(j+1)a^{[j+1]}\Delta t^{j}

    b^{'}(t) = \sum\limits_{j=0}^{n-1}(j+1)b^{[j+1]}\Delta t^{j}

    c^{'}(t) = \sum\limits_{j=0}^{n-1}(j+1)c^{[j+1]}\Delta t^{j}
 

Based on those fomulas, we can obtain the different basic functions of the taylor series.

Addition and subtraction
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a(t) = b(t) \pm c(t)`, we have the iterative fomulas:

.. math::

    a^{[n]}(t) = b^{[n]}(t) \pm c^{[n]}(t)

Multiplication
>>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a(t) = b(t)c(t)`, we have the iterative fomulas:

.. math::

    a^{[n]}(t) = \sum\limits_{j=0}^{n}b^{[j]}(t)c^{[n-j]}(t)

Division
>>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a(t) = b(t)/c(t)`, which constantly equals :math:`a(t)c(t) = b(t)`, we have the iterative fomulas:

.. math::

    a^{[n]}(t)=\frac{1}{c^{[0]}(t)}\left[b^{[n]}(t)-\sum_{j=1}^n a^{[j]}(t) c^{[n-j]}(t)\right]

Logarithm
>>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a(t) = \ln{b(t)}`, we have the iterative fomulas:

.. math::

    a^{[n]}(t)=\frac{1}{n b^{[0]}(t)}\left[n b^{[n]}(t)-\sum_{j=1}^{n-1} j  a^{[j]}(t) b^{[n-j]}(t) \right]

For the :math:`a(t) = \log_{c} b(t)`, we can use the logarithm change of base rule :math:`\log_{c} b(t) = \dfrac{\ln b(t)}{\ln c}`, which gives:

.. math::

    a^{[n]}(t)=\frac{1}{nb^{[0]}(t)\ln c }\left[n b^{[n]}(t)-\sum_{j=1}^{n-1} j   a^{[j]}(t) b^{[n-j]}(t) \right]

Exponential
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a(t) = e^{b(t)}`, we have the iterative fomulas:

.. math::

    a^{[n]}(t)=\frac{1}{n} \sum_{j=0}^{n-1} (n-j) a^{[j]}(t) b^{[n-j]}(t)

For the :math:`a(t) = c^{b(t)}`, we can use :math:`c^{b(t)} = e^{b(t)\ln{c}}`, which gives:

.. math::

    a^{[n]}(t)=\frac{\ln{c}}{n} \sum_{j=0}^{n-1} (n-j) a^{[j]}(t) b^{[n-j]}(t)

Sine and cosine
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a_s(t) = \sin(b(t))` and the :math:`a_c(t) = \cos(b(t))` we have the iterative fomulas:

.. math::

    a_s^{[n]}(t) & = \dfrac{1}{n}\sum_{j=0}^{n-1}(n-j)a_c^{[j]}(t)b^{[n-j]}(t)

    a_c^{[n]}(t) & = -\dfrac{1}{n}\sum_{j=0}^{n-1}(n-j)a_s^{[j]}(t)b^{[n-j]}(t)

Hyperbolic sine and hyperbolic cosine
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

For the :math:`a_{sh}(t) = \sinh(b(t))` and the :math:`a_{ch}(t) = \cosh(b(t))` we have the iterative fomulas:

.. math::

    a_{sh}^{[n]}(t) & = \dfrac{1}{n}\sum_{j=0}^{n-1}(n-j)a_{ch}^{[j]}(t)b^{[n-j]}(t)

    a_{ch}^{[n]}(t) & = \dfrac{1}{n}\sum_{j=0}^{n-1}(n-j)a_{sh}^{[j]}(t)b^{[n-j]}(t)

Tangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \tan (b(t))`, let:

.. math::

    h(b(t)) = \dfrac{1}{\cos^{2}(b(t))}

we have:

.. math::

    a^{'}(t) & = \tan^{'} (b(t)) b^{'}(t)

    & = h(b(t)) b^{'}(t)

    h^{'}(b(t)) & = \left( \dfrac{1}{\cos^{2}(b(t))} \right)^{'}

    & = \dfrac{2\sin(b(t))b^{'}(t)}{\cos^{3}(b(t))}

    & = 2 a(t) h(b(t)) b^{'}(t)

    & = 2 a(t) a^{'}(t)

.. note::
    The `'` of a quantity denote the derivative of the this quantity about the quantity in its parentheses. For example, :math:`a^{'}(t)` means :math:`\dfrac{\partial a(t)}{\partial t}` and :math:`\tan^{'}(b(t))` means :math:`\dfrac{\partial tan(b(t))}{\partial b(t)}`.

with the multiplication of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)h^{[j]}(t)b^{[n-j]}(t)

    h^{[n]}(t) & = \dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[j]}(t)a^{[n-j]}(t)

Hyperbolic tangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \tanh (b(t))`, let:

.. math::

    h(b(t)) = \dfrac{1}{\cosh^{2}(b(t))}

we have:

.. math::

    a^{'}(t) & = \tanh^{'} (b(t)) b^{'}(t)

    & = h(b(t)) b^{'}(t)

    h^{'}(b(t)) & = \left( \dfrac{1}{\cosh^{2}(b(t))} \right)^{'}

    & = -\dfrac{2\sinh(b(t))b^{'}(t)}{\cosh^{3}(b(t))}

    & = -2 a(t) h(b(t)) b^{'}(t)

    & = -2 a(t) a^{'}(t)

with the multiplication of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)h^{[j]}(t)b^{[n-j]}(t)

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[j]}(t)a^{[n-j]}(t)


Cotangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \cot (b(t))`, let:

.. math::

    h(b(t)) = -\dfrac{1}{\sin^{2}(b(t))}

we have:

.. math::

    a^{'}(t) & = \cot^{'} (b(t)) b^{'}(t)

    & = h(b(t)) b^{'}(t)

    h^{'}(b(t)) & = \left( -\dfrac{1}{\sin^{2}(b(t))} \right)^{'}

    & = \dfrac{2\cos(b(t))b^{'}(t)}{\sin^{3}(b(t))}

    & = -2 a(t) h(b(t)) b^{'}(t)

    & = -2 a(t) a^{'}(t)

with the multiplication of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)h^{[j]}(t)b^{[n-j]}(t)

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[j]}(t)a^{[n-j]}(t)

Hyperbolic cotangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \coth (b(t))`, let:

.. math::

    h(b(t)) = -\dfrac{1}{\sinh^{2}(b(t))}

we have:

.. math::

    a^{'}(t) & = \cot^{'} (b(t)) b^{'}(t)

    & = h(b(t)) b^{'}(t)

    h^{'}(b(t)) & = \left( \dfrac{1}{\sinh^{2}(b(t))} \right)^{'}

    & = \dfrac{2\cosh(b(t))b^{'}(t)}{\sinh^{3}(b(t))}

    & = -2 a(t) h(b(t)) b^{'}(t)

    & = -2 a(t) a^{'}(t)

with the multiplication of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)h^{[j]}(t)b^{[n-j]}(t)

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[j]}(t)a^{[n-j]}(t)

Inverse sine
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \arcsin (b(t))`, let:

.. math::

    h(b(t)) = \sqrt{1-b^2(t)}

we have:

.. math::

    a^{'}(t) & = \arcsin^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))}

    h^{'}(b(t)) & = \left( \sqrt{1-b^2(t)} \right)^{'}

    & = -\dfrac{b^{'}(t)b(t)}{\sqrt{1-b^2(t)}}

    & = -a^{'}(t)b(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{nh^{[0]}(t)}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = -\dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[n-j]}(t)b^{[j]}(t)

Inverse hyperbolic sine
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \operatorname{arcsinh} (b(t))`, let:

.. math::

    h(b(t)) = \sqrt{1+b^2(t)}

we have:

.. math::

    a^{'}(t) & = \operatorname{arcsinh}^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))}

    h^{'}(b(t)) & = \left( \sqrt{1+b^2(t)} \right)^{'}

    & = \dfrac{b^{'}(t)b(t)}{\sqrt{1+b^2(t)}}

    & = a^{'}(t)b(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{nh^{[0]}(t)}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[n-j]}(t)b^{[j]}(t)

Inverse cosine
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \arccos (b(t))`, let:

.. math::

    h(b(t)) = -\sqrt{1-b^2(t)}

we have:

.. math::

    a^{'}(t) & = \arccos^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))}

    h^{'}(b(t)) & = \left( -\sqrt{1-b^2(t)} \right)^{'}

    & = \dfrac{b^{'}(t)b(t)}{\sqrt{1-b^2(t)}}

    & = -a^{'}(t)b(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{nh^{[0]}(t)}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = -\dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[n-j]}(t)b^{[j]}(t)

Inverse hyperbolic cosine
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \operatorname{arccosh} (b(t))`, let:

.. math::

    h(b(t)) = \sqrt{b^2(t)-1}

we have:

.. math::

    a^{'}(t) & = \operatorname{arccosh}^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))}

    h^{'}(b(t)) & = \left( \sqrt{b^2(t)-1} \right)^{'}

    & = \dfrac{b^{'}(t)b(t)}{\sqrt{1-b^2(t)}}

    & = a^{'}(t)b(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{nh^{[0]}(t)}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = \dfrac{1}{n} \sum\limits_{j=0}^{n-1}(n-j)a^{[n-j]}(t)b^{[j]}(t)


Inverse tangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \arctan (b(t))`, let:

.. math::

    h(b(t)) = b^2(t)

we have:

.. math::

    a^{'}(t) & = \arctan^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))+1}

    h^{'}(b(t)) & = \left( b^2(t) \right)^{'}

    & = 2b(t)b^{'}(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{n(1+h^{[0]}(t))}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = \dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)b^{[j]}(t)b^{[n-j]}(t)

Inverse hyperbolic tangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \operatorname{arctanh} (b(t))`, let:

.. math::

    h(b(t)) = -b^2(t)

we have:

.. math::

    a^{'}(t) & = \operatorname{arctanh}^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))+1}

    h^{'}(b(t)) & = \left( -b^2(t) \right)^{'}

    & = -2b(t)b^{'}(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{n(1+h^{[0]}(t))}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)b^{[j]}(t)b^{[n-j]}(t)

Inverse cotangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \operatorname{arccot} (b(t))`, let:

.. math::

    h(b(t)) = -b^2(t)

we have:

.. math::

    a^{'}(t) & = \operatorname{arccot}^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))-1}

    h^{'}(b(t)) & = \left( -b^2(t) \right)^{'}

    & = -2b(t)b^{'}(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{n(-1+h^{[0]}(t))}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)b^{[j]}(t)b^{[n-j]}(t)

Inverse hyperbolic cotangent
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

For :math:`a(t) = \operatorname{arccoth} (b(t))`, let:

.. math::

    h(b(t)) = -b^2(t)

we have:

.. math::

    a^{'}(t) & = \operatorname{arccoth}^{'} (b(t)) b^{'}(t)

    & = \dfrac{b^{'}(t)}{h(b(t))+1}

    h^{'}(b(t)) & = \left( -b^2(t) \right)^{'}

    & = -2b(t)b^{'}(t)

with the multiplication and division of taylor series, we have the final iterative fomulas:

.. math::

    a^{[n]}(t) & = \frac{1}{n(1+h^{[0]}(t))}\left[nb^{[n]}(t)-\sum_{j=1}^{n-1}j a^{[j]}(t) h^{[n-j]}(t)\right]

    h^{[n]}(t) & = -\dfrac{2}{n} \sum\limits_{j=0}^{n-1}(n-j)b^{[j]}(t)b^{[n-j]}(t)


Self-Adaptive Algorithm
-----------------------

Using the above basic arithmetic of Taylor series, we can get the arbitrary derivative of the function. When the total calculation time of the dynamical system is decided, the order of the Taylor series method and word size of the program is settled. But when the calculation progresses over time, the average noise will growth which make the order and word size are meaningless to maintain with the same value. So we need to use the self-adaptive algorithm to speed up the calculation.

Following Qin and Liao, the word size of the program is defined as:

.. math::

    N_s=\left\lceil\frac{\gamma \kappa\left(T_c-t^{\prime}\right)}{\ln 10}-\log _{10} \varepsilon_c\right\rceil

where the :math:`\gamma` is safety factor greater than 1, :math:`\kappa` is the noise-growing exponent, :math:`T_c` is the critical predictable time, :math:`t^{\prime}` is the current time, :math:`\varepsilon_c` is the critical predictable error, and :math:`\lceil \cdot \rceil` is ceiling function.

the order of the Taylor series method is defined as:

.. math::

    M = \left\lceil -1.5\log_{10}(tol) \right\rceil = \left\lceil 1.5N_s \right\rceil

It should be noted that the number 1.5 in order of the Taylor series method is a empirical value, which can be changed according to the specific situation. According to :cite:`jorba2005software`, for non-multi process the order of the Taylor series method is defined as:

.. math::

    M = \left\lceil 1.15N_s + 1 \right\rceil

Following Barrio et al. :cite:`barrio2005vsvo`, the optimal time stepsize is given by:

.. math::

    \Delta t=\min \left(\frac{\operatorname{tol}^{\frac{1}{M}}}{\left\|x_i^{[M-1]}(t)\right\|_{\infty}^{\frac{1}{M-1}}}, \frac{\operatorname{tol}^{\frac{1}{M+1}}}{\left\|x_i^{[M]}(t)\right\|_{\infty}^{\frac{1}{M}}}\right)

References
----------
.. bibliography:: theory.bib
   :style: unsrt
   :labelprefix: T
