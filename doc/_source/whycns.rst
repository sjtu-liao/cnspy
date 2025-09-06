

Why CNS?
==========

A chaotic dynamical system has a popular property called sensitive dependence on initial conditions (SDIC), which means that the system is very sensitive to the initial conditions. It not only has the sensitivity dependence on initial conditions (SDIC) but also possesses the sensitivity dependence on numerical algorithms (SDNA), as reported by Lorenz :cite:`lorenz1989computational`. Thus, for a chaotic dynamical system, calculated trajectories of computer simulations are very sensitive to the initial conditions and numerical algorithms. This is the reason why we need to develop a new numerical algorithm for chaotic dynamical systems.

The results of chaotic dynamical systems given by the CNS can be regarded as a “clean” benchmark solution, which is the main purpose of proposing
this CNS strategy. By contrast, solely adopting the Taylor series method :cite:`corliss1982solving` to solve a chaotic system for high precision, one usually does not focus on the “critical predictable time” Tc and thus obtain a mixture of the “true” physical solution and the “false” numerical noise, which are mostly at the same order of magnitude, since the background numerical noise of a simulation of chaos should increase exponentially (and quickly) until to the same level of “true” physical solution, which is not considered by traditional numerical strategies.

Works using CNS
------------------

The CNS has already been used in varies chaotic dynamcial systems. 

Classical Lorenz 69 system
>>>>>>>>>>>>>>>>>>>>>>>>>>>

A reproducible and convergent numerical solution of the Lorenz system is obtained within :math:`t \in [0,10000]` by Liao and Wang :cite:`liao2014mathematically`. 

Three-body problem
>>>>>>>>>>>>>>>>>>>>>>>>>>>

Significantly, the CNS strategy was applied to investigate the periodic orbits of the famous three-body problem, and more than 2000 brand-new families of periodic orbits were discovered successfully by Li et al :cite:`li2017more,li2018over,li2019collisionless`. (For more gifs and pictures please refer to `Periodic orbits for Newtonian planar three-body problem <https://numericaltank.sjtu.edu.cn/three-body/three-body.htm>`_)

Spatiotemporal chaos
>>>>>>>>>>>>>>>>>>>>>>>>>>>

As for spatiotemporal chaos,  Hu & Liao :cite:`hu2020risks` and Qin & Liao :cite:`qin2020influence` proposed an efficient CNS strategy utilized in physical space to numerically solve the 1D complex Ginzburg-Landau equation (CGLE) and the damped driven sine-Gordon equation (SGE), respectively.

Turbulence
>>>>>>>>>>>>>>>>>>>>>>>>>>>

Qin & Liao :cite:`qin2022large` provide rigorous evidence that numerical noises as a kind of tiny artificial stochastic disturbances have quantitatively and qualitatively large-scale influences on a sustained turbulence.

.. raw:: html

   <video width="640" height="360" controls>
      <source src="https://raw.githubusercontent.com/sjtu-liao/RBC/main/RBC-mv.mp4" type="video/mp4">
   </video>

References
----------
.. bibliography:: whycns.bib
   :style: unsrt
   :labelprefix: W