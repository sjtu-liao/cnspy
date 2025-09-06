"""
Clean Numerical Simulation code generater using sympy as symbolic operation tools.

Bo Zhang 2023 <zilpher@sjtu.edu.cn>
"""
from .func import tree
from .funcs import trees
from .method import method, tsm
from .adaptive import adaptive, adaptive_parameters
from .printer import Cprinter_helper, Foldermaker
from .lyap import lyapunov
from .cns import cns


__all__=['tree','trees','method','tsm','adaptive','adaptive_parameters','Cprinter_helper','Foldermaker','lyapunov','cns']
__version__ = '1.0.0'