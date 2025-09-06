import sympy as sp
import numpy as np

class tree():
    """
    This class will transfer one function into a list of functions that only contains
    binary operation and two nodes.
    
    After initialization of the class, one can evaluate it to get the list of new 
    functions.
    
    Examples
    ======== 
        
        >>> import func
        
        >>> a,b,c,d,e,f,g,h,i,j,k,l,m,n = sp.symbols(r'a,b,c,d,e,f,g,h,i,j,k,l,m,n')
        
        >>> fun = a*(b+c)+sp.sin(d)-sp.ln(f)/g+h*sp.cos(i)*(j+k)/(l+m*n)
        
        >>> tree = func.tree(fun)
        
        >>> tree.evaluate()
        
        >>> tree.bnode
        
    """
    __cl: list
    __node: list
    __bnode: list
    __level: list
    __maxnum: int
    __maxsavernum: int
    __maxlevel: int
    __tmplevel: int
    
    def __init__(self,func):
        cl = []
        node = []
        level = []
        for arg in sp.postorder_traversal(func):
            cl.append(arg.func)
            node.append(arg)
            level.append(0)
        self.__cl=cl
        self.__node=node
        self.__level=level
        self.__maxnum = len(node)
        self.__tmplevel = 1
        self.__bnode = []

    def __which_node(self,subnode,index):
        """
        Find subnode's global index.
        
        Examples
        ======== 
            
        >>> self.__which_node(node.args[0],3)
            
        """
        if index == self.__maxnum:
            return self.__maxnum
        for i in range(0,index):
            if subnode == self.__node[i]:
                return i
        
    def __find_level(self,node,index):
        """
        Recursive private function to give the node's level.
        
        Examples
        ========
            
        >>> self.__find_level(self.__node[i].args,index)
            
        """
        for i in range(0,len(node)):
            if node[i]!=() and self.__tmplevel <= self.__level[self.__which_node(node[i],index)]:
                self.__tmplevel = self.__level[self.__which_node(node[i],index)]+1
    
    def __break_nonbinary(self):
        """
        Since the node Add and Mul in sympy can have more than two subnodes, 
        so this function will break the Add and Mul note in to binary subnodes.
        The return will fill the __bnode list which mean binary_node.
        
        Examples
        ========
            
        >>> self.break_nonbinary()
            
        """
        num = 0
        for le in range(1,self.__maxlevel+1):
            for i in range(0,self.__maxnum):
                if self.__level[i] == le:
                    if len(self.__node[i].args) >= 2:
                        num += len(self.__node[i].args)-1
                    elif len(self.__node[i].args) == 1:
                        num += 1
        saver=sp.symarray('saver', num)
        self.__maxsavernum = num
        savers = []
        saverindex = 0     
        for le in range(1,self.__maxlevel+1):
            for i in range(0,self.__maxnum):
                if self.__level[i] == le:
                    if len(self.__node[i].args) > 2:
                        if self.__cl[i] == sp.core.add.Add:
                            tmp = self.__node[i].args[0]+self.__node[i].args[1]
                            savers.append(tmp)
                            self.__bnode.append([self.__cl[i],self.__node[i].args[0],self.__node[i].args[1],saver[saverindex]])
                            for j in range(1,len(self.__node[i].args)-1):
                                tmp += self.__node[i].args[j+1]
                                savers.append(tmp)
                                self.__bnode.append([self.__cl[i],saver[saverindex],self.__node[i].args[j+1],saver[saverindex+1]])
                                saverindex += 1
                            saverindex += 1
                        elif self.__cl[i] == sp.core.mul.Mul:
                            tmp = self.__node[i].args[0]*self.__node[i].args[1]
                            savers.append(tmp)
                            self.__bnode.append([self.__cl[i],self.__node[i].args[0],self.__node[i].args[1],saver[saverindex]])
                            for j in range(1,len(self.__node[i].args)-1):
                                tmp *= self.__node[i].args[j+1]
                                savers.append(tmp)
                                self.__bnode.append([self.__cl[i],saver[saverindex],self.__node[i].args[j+1],saver[saverindex+1]])
                                saverindex += 1
                            saverindex += 1
                        else:
                            raise Exception("An operation has multi-subnodes that not considered. In the present, we only consider the sympy.core.add.Add and sympy.core.mul.Mul.")
                    elif len(self.__node[i].args) == 2:
                        savers.append(self.__node[i])
                        self.__bnode.append([self.__cl[i],self.__node[i].args[0],self.__node[i].args[1],saver[saverindex]])
                        saverindex += 1
                    elif len(self.__node[i].args) == 1:
                        savers.append(self.__node[i])
                        self.__bnode.append([self.__cl[i],self.__node[i].args[0],saver[saverindex]])
                        saverindex += 1
                    else:
                        raise Exception("The level {le} has zero arg")
        for i in range(0,num):
            for j in range(1,len(self.__bnode[i])-1):
                if len(self.__bnode[i][j].args) != 0:
                    for k in range(0,num):
                        if self.__bnode[i][j] == savers[k]:
                            self.__bnode[i][j]=self.__bnode[k][-1] 
            
    @property
    def cl(self):
        return self.__cl
    @property
    def node(self):
        return self.__node
    @property
    def bnode(self):
        return self.__bnode
    @property
    def level(self):
        return self.__level
    @property
    def maxnum(self):
        return self.__maxnum
    @property
    def maxsavernum(self):
        return self.__maxsavernum
    @property
    def maxlevel(self):
        return self.__maxlevel
    
    def replace_saver(self,startnum):
        """
        Replace the saver variable number in different functions from a start index number..
        
        Examples
        ========
            
        >>> self.replace_saver(10)
        
        """
        saver = sp.symarray('saver', self.__maxsavernum+startnum)
        newsaver = sp.symarray('newsaver', self.__maxsavernum+startnum)
        for i in range(0,self.__maxsavernum):
            for j in range(1,len(self.__bnode[i])):
                for k in range(0,i+1):
                    self.__bnode[i][j] = self.__bnode[i][j].subs(saver[k],newsaver[k+startnum])
        
        for i in range(0,self.__maxsavernum):
            for j in range(1,len(self.__bnode[i])):
                for k in range(startnum,i+1+startnum):
                    self.__bnode[i][j] = self.__bnode[i][j].subs(newsaver[k],saver[k])
        
    def evaluate(self):
        """
        Evaluating the function and fullfill the unknown properties.
        
        """
        for i in range(0,self.__maxnum):
            if self.__node[i].args!=():
                self.__level[i] = 1
        for i in range(0,self.__maxnum):
            self.__tmplevel = self.__level[i]
            if self.__level[i] == 1:
                self.__find_level(self.__node[i].args,i) 
                self.__level[i] = self.__tmplevel
        self.__maxlevel = max(self.__level)
        self.__break_nonbinary()

        if self.__maxsavernum == 0:
            num = 1
            saver=sp.symarray('saver', num)
            self.__maxsavernum = num
            self.__bnode.append([sp.core.relational.Eq,self.__node[0],saver[0]])
            
        