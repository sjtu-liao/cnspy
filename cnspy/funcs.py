import sympy as sp
import numpy as np
import re
import copy
from .func import tree

class trees():
    """
    This class will transfer a set of functions into a list of functions that only
    contains binary operation and two nodes.
    
    After initialization of the class, one can get the list of new functions.
    
    For example: 
        
        >>> import funcs
        
        >>> import sympy as sp
        
        >>> a,b,c,d,e,f,g,h,i,j,k,l,m,n = sp.symbols(r'a,b,c,d,e,f,g,h,i,j,k,l,m,n')
        
        >>> funs = [a*(b+c)+sp.sin(d), -sp.ln(f)/g, h*sp.cos(i)*(j+k)/(l+m*n)]
        
        >>> trees = funcs.trees(funs)
        
        >>> trees.functree[0].bnode
        
    """
    __funcnum: int
    __varnum: int
    __parnum: int
    __maxsavernum: int
    __functree: list
    __bnode: list
    __slist: list
    __var: list
    __par: list
    
    def __init__(self,funcs:list,vars:list,paras:list):
        self.__funcnum=len(funcs)
        self.__varnum=len(vars)
        self.__parnum=len(paras)
        self.__var = []
        self.__par = []
        for i in range(0,self.__varnum):
            self.__var.append(vars[i])
        for i in range(0,self.__parnum):        
            self.__par.append(paras[i])
        Firsttime = True
        startnum = 0
        incnum = 0
        self.__functree = []
        self.__bnode = []
        self.__var_saver = []
        for fun in funcs:
            obj = tree(fun)
            obj.evaluate()
            if Firsttime:
                incnum = obj.maxsavernum
                self.__functree.append(obj)
                self.__var_saver.append(obj.bnode[incnum-1][-1])
                Firsttime = False
            else:
                startnum += incnum
                incnum = obj.maxsavernum
                obj.replace_saver(startnum)
                self.__functree.append(obj)
                self.__var_saver.append(obj.bnode[incnum-1][-1])
        self.__maxsavernum = startnum + incnum
        for i in range(0,self.__funcnum):
            for j in range(0,len(self.__functree[i].bnode)):
                self.__bnode.append(self.__functree[i].bnode[j])
        self.__slist=[]
        for i in range(0,len(self.__bnode)):
            self.__slist.append(copy.deepcopy(self.__bnode[i]))
        for i in range(0,len(self.__slist)):
            for j in range(1,len(self.__slist[i])):
                self.__slist[i][j] = False
        self.reorde_bnode()
        self.series_list()
        self.constant()
        self.__maxsavernum = int(re.findall(r'saver_\d+', str(self.__bnode[-1][-1]))[0].split("_")[1])+1
        #self.var_saver()

    def __reorder_savers(self,lst):
        savers = []
        for expr in lst:
            for node in expr:
                if re.findall(r'saver_\d+', str(node)) != []:
                    savers.append(re.findall(r'saver_\d+', str(node))[0])
        int_values = [int(s.split("_")[1]) for s in savers]
        #to find the variables
        self.__int_oldlist = copy.deepcopy(int_values)
        int_newlist = copy.deepcopy(int_values)
        #int_newlist = int_values
        tmp_value = 0
        free_list = []
        global_i = 0
        for i in range(0,len(int_values)):
            #value is not equal to the golbal index (natural number order) and not in the free list
            if int_values[i] != global_i and i not in free_list:
                tmp_value = copy.deepcopy(int_values[i])
                int_newlist[i] = copy.deepcopy(global_i)         
                for j in range(i,len(int_values)):
                    if int_values[j] == tmp_value:
                        int_newlist[j] = copy.deepcopy(global_i)
                        free_list.append(j)
                global_i += 1 
            #value is equal to the global index and not in the free list
            elif int_values[i] == global_i and i not in free_list:
                global_i += 1
                for j in range(i,len(int_values)):
                    if int_values[j] == int_values[i]:
                        free_list.append(j)
        return int_newlist
                
    def reorde_bnode(self):
        """
        Reording bnode, since the series of the `sin(x)` and `cos(x)`, the `sinh(x)` and 
        `cosh(x)`, we need to check if there are the same `x`, and make them be together
        in the first of the bnode list. If find a function that need asistant function, it will
        change the bnode list to add the asistant series.
        
        [ ....                                   [ ....
         [sympy.core.add.Add,a,b,c]               [sympy.core.add.Add,a,b,c]
         [sin,c,d]                                [sin,c,d]
         [sympy.core.add.Mul,e,f,g]               [cos,c,i]
         [cosh,g,h]                    ---->      [sinh,g,j]
         [cos,c,i]                                [cosh,g,h]
         [sympy.core.add.Add,i,j,k]               [sympy.core.add.Mul,e,f,g]
          ....                                    [sympy.core.add.Add,i,j,k]
        ]                                          ....]     

        Examples
        ========
            
        >>> self.reorder_bnode()
        
        """
        #firstly, to reording the function pair.
        for func in [[sp.sin,sp.cos],[sp.sinh,sp.cosh]]:
            i = 0
            while i<len(self.__bnode):
                not_found = True
                #if the first found func is sin,sinh.
                if self.__bnode[i][0] == func[0]:
                    #then search for the cos,cosh pair,and locate their locations.
                    for j in range(i+1,len(self.__bnode)):
                        # only if they share the same input.
                        if self.__bnode[j][0] == func[1] and self.__bnode[j][1] == self.__bnode[i][1]:
                            #remove the latter pair from the previous location and add it after the first pair. 
                            insertposition = i+1
                            self.__bnode.insert(insertposition,self.__bnode.pop(j))                                                                   
                            i += 1
                            not_found = False
                            break
                    if not_found:
                        #no pairs, add a dumy saver after the first pair.
                        insertposition = i+1
                        dumysaver = sp.symbols(f"saver_{self.__maxsavernum}")
                        self.__bnode.insert(insertposition,[func[1],self.__bnode[i][1],dumysaver])
                        self.__maxsavernum += 1
                        i += 1
                # if the first found func is cos,cosh
                elif self.__bnode[i][0] == func[1]:
                    #then search for the sin,sinh pair, and locate their locations.
                    for j in range(i+1,len(self.__bnode)):
                        # only if their share the same input.
                        if self.__bnode[j][0] == func[0] and self.__bnode[j][1] == self.__bnode[i][1]:
                            #if the first pair is at the begining of the bnode list, the insert position is 0, else the insert position is i-1
                            insertposition = i
                            #remove the latter pair from the previous location and add it before the first pair.                                 
                            self.__bnode.insert(insertposition,self.__bnode.pop(j))
                            i += 1
                            not_found = False
                            break
                    if not_found:
                        #if the first pair is at the begining of the bnode list, the insert position is 0, else the insert position is i-1
                        insertposition = i
                        #no pairs, add a dumy saver before the first pair.
                        dumysaver = sp.symbols(f"saver_{self.__maxsavernum}")
                        self.__bnode.insert(insertposition,[func[0],self.__bnode[i][1],dumysaver])
                        self.__maxsavernum += 1
                        i += 1
                #no iteration function, pass to the next iteration
                else:
                    pass
                i += 1
        
        #Then, we rearrange the function that need asistant function: 
        #tan, tanh, cot, coth, asin, asinh, acos, acosh, atan, atanh, acot, acoth
        for func in [sp.tan,sp.tanh,sp.cot,sp.coth,sp.asin,sp.asinh,sp.acos,sp.acosh,sp.atan,sp.atanh,sp.acot,sp.acoth]:
            i = 0
            while i<len(self.__bnode):
                if self.__bnode[i][0] == func:
                    dumysaver = sp.symbols(f"saver_{self.__maxsavernum}")
                    #The functions here are all only contains one input, so the insert position is fixed at 2.
                    self.__bnode[i].insert(2,dumysaver)
                    self.__maxsavernum += 1
                i += 1
                

        #delete the same bnodes
        index_of_duplicate = []
        for i in range(0,len(self.__bnode)):
            for j in range(i+1,len(self.__bnode)):
                if len(self.__bnode[i]) == 4 and len(self.__bnode[j]) == 4:
                    if self.bnode[i][0] == self.bnode[j][0] and self.bnode[i][1] == self.bnode[j][1] and self.bnode[i][2] == self.bnode[j][2]:
                        for _i in range(i+1,len(self.__bnode)):
                            for _j in range(1,len(self.__bnode[_i])-1): 
                              if self.__bnode[_i][_j] == self.__bnode[j][3]:
                                  self.__bnode[_i][_j] = self.__bnode[i][3]
                        index_of_duplicate.append(j)
                elif len(self.__bnode[i]) == 3 and len(self.__bnode[j]) == 3:
                    if self.bnode[i][0] == self.bnode[j][0] and self.bnode[i][1] == self.bnode[j][1]:
                        for _i in range(i+1,len(self.__bnode)):
                            for _j in range(1,len(self.__bnode[_i])-1): 
                              if self.__bnode[_i][_j] == self.__bnode[j][2]:
                                  self.__bnode[_i][_j] = self.__bnode[i][2]
                        index_of_duplicate.append(j)
        #offset = 0
        #for i in index_of_duplicate:
        #    self.__bnode.pop(i-offset)
        #    offset += 1
        #    self.__maxsavernum -= 1
        # just find the duplicate saver and delete it without changing the index of the bnode list.
        popsaver = []
        for i in index_of_duplicate:
            popsaver.append(self.__bnode[i][-1])
        for _ in range(len(index_of_duplicate)):
            for i in range(len(self.__bnode)):
                if self.__bnode[i][-1] in popsaver:
                    self.__bnode.pop(i)
                    self.__maxsavernum -= 1
                    break

        #delete the duplicate binary output node
        #for i in range(0,len(self.__bnode)-2):
        #    if self.__bnode[i][-1] == self.__bnode[i+1][-1]:
        #        newsaver = sp.symbols(f"saver_{self.__maxsavernum}")
        #        self.__bnode[i+1][-1] = newsaver
        #        self.__maxsavernum += 1
        #        for j in range(i+2,len(self.__bnode)):
        #            for k in range(1,len(self.__bnode[j])):
        #                if self.__bnode[i][-1] == self.__bnode[j][k]:
        #                    self.__bnode[j][k] = newsaver
        
        #finally, we reindexing the saver matrix.
        int_new_list = self.__reorder_savers(self.__bnode)
        
        #this part can find the saver that represents the variable
        var_saver_num = [int(str(s).split("_")[1]) for s in self.__var_saver]
        self.var_saver = []
        for j in range(0,len(var_saver_num)):
            for i in range(0,len(self.__int_oldlist)):
                if var_saver_num[j] == self.__int_oldlist[i]:
                    self.var_saver.append(int_new_list[i])

        #this part just do the reordering the saver
        int_new_list_index = 0
        for i in range(0,len(self.__bnode)):
            for j in range(1,len(self.__bnode[i])):
                if re.findall(r'saver_\d+', str(self.__bnode[i][j])) != []:
                    self.__bnode[i][j] = sp.symbols(f"saver_{int_new_list[int_new_list_index]}")
                    int_new_list_index += 1

                    
        self.__slist=[]
        for i in range(0,len(self.__bnode)):
            self.__slist.append(copy.deepcopy(self.__bnode[i]))
        for i in range(0,len(self.__slist)):
            for j in range(1,len(self.__slist[i])):
                self.__slist[i][j] = False

    #give a list of each variable that are series(connected to time and variables)
    def series_list(self):
        #for each elements in bnode list
        for i in range(0,len(self.__bnode)):
            for j in range(1,len(self.__bnode[i])):
                #for each variable
                for k in range(0,self.__varnum):
                    #if is variable, then series list is True
                    if self.__bnode[i][j] == self.__var[k]:
                        self.__slist[i][j] = True
                for k in range(0,self.__parnum):
                    #if is parameter, then series list is False
                    if self.__bnode[i][j] == self.__par[k]:
                        self.__slist[i][j] = False
                if self.__bnode[i][j] == sp.Symbol('t'):
                    self.__slist[i][j] = True
        #if the operation contains help functions, the help functions will be series, and change it in the following bnode list. 
        for i in range(0,len(self.__bnode)):
            if self.__bnode[i][0] in [sp.tan,sp.tanh,sp.cot,sp.coth,sp.asin,sp.asinh,sp.acos,sp.acosh,sp.atan,sp.atanh,sp.acot,sp.acoth]:
                self.__slist[i][2] = True
        for i in range(0,len(self.__bnode)):
            #if self.__bnode[i][0] == sp.core.add.Add or self.__bnode[i][0] == sp.core.mul.Mul:
            for j in range(1,len(self.__bnode[i])):
                #if the operation contains series, the out put will be series, and change it in the following bnode list.
                if self.__slist[i][j]:
                    self.__slist[i][-1] = True
                    for _i in range(i+1,len(self.__bnode)):
                        for _j in range(1,len(self.__bnode[_i])):
                            if self.__bnode[_i][_j] == self.__bnode[i][-1]:
                                self.__slist[_i][_j] = True
                else:
                    pass

    def constant(self):
        constants_list = []
        constants_index = []
        saver_list = sp.symarray('saver', self.__maxsavernum)
        for i in range(0,len(self.__bnode)):
            for j in range(1,len(self.__bnode[i])):
                if self.__slist[i][j] == False and self.__bnode[i][j] not in saver_list and self.__bnode[i][j] not in self.__par :
                    constants_list.append(str(self.__bnode[i][j]))
                    constants_index.append([i,j])
        self.__constants = list(set(constants_list))
        self.__constantnum = len(self.__constants)
        constant_sp_list = sp.symarray('constant', self.__constantnum)
        for i in range(0,self.__constantnum):
            for j in range(0,len(constants_index)):
                if constants_list[j] == self.__constants[i]:
                    self.__bnode[constants_index[j][0]][constants_index[j][1]] = constant_sp_list[i]
    
    #def var_saver(self):
    #    self.var_saver_list = []
    #    flag = True
    #    for i in range(0,len(self.__bnode)):
    #        for j in range(i+1,len(self.__bnode)):
    #            for k in range(1,len(self.__bnode[j])):
    #                if self.__bnode[i][-1] == self.__bnode[j][k]:
    #                    flag = False
    #        if flag:
    #            self.var_saver_list.append(i)
    #        flag = True
    #    for i in range(0,self.var_saver_list-1):
    #        if self.var_saver_list[i]+1 == self.var_saver_list[i+1] and (self.__bnode[i][0] == sp.sin or self.__bnode[i][0] == sp.sinh) and self.__bnode[i][0] not in self.__var_op:
    #            self.var_saver_list.remove(i)
    #            
    #    for i in self.var_saver_list:
    #        if self.__bnode[i][0] not in self.__var_op:
    #           self.var_saver_list.remove(i)
            

    @property
    def functree(self):
        return self.__functree
    @property
    def funcnum(self):
        return self.__funcnum
    @property
    def maxsavernum(self):
        return self.__maxsavernum
    @property
    def bnode(self):
        return self.__bnode  
    @property
    def slist(self):
        return self.__slist    
    @property
    def var(self):
        return self.__var
    @property
    def par(self):
        return self.__par
    @property
    def constantnum(self):
        return self.__constantnum
    @property
    def constants(self):
        return self.__constants
        
        

