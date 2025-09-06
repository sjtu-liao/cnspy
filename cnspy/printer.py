from importlib.metadata import version
import sympy as sp
import os

class Cprinter_helper():
    """
    Help function of C code generater.
    """
    def __init__(self):
        pass
    def indent(self,level):
        '''
        Return a string with the appropriate amount of whitespace
        for the specified indentation level (in number of tabs).
        '''
        return '\t' * level
    def brief_info(self):
        '''
        Return a string that contains brief information of the file.
        '''
        return  f"/* This file is part of the cnspy({version('cnspy')}).\n" +\
                f" * Copyright (C) 2023 Shanghai Jiao Tong University and Authors.\n" +\
                f" * Contact: Bo Zhang <zilpher@sjtu.edu.cn>\n" +\
                f" * Licensed under the <MPL-2.0>;\n" +\
                f" * you may not use this file except in compliance with the License.\n" +\
                f" * You may obtain a copy of the License at <https://www.mozilla.org/en-US/MPL/2.0/>.\n" +\
                f" */\n\n"
    
class Foldermaker():
    """
    Creating the Folder.
    """
    def __init__(self,dirpath):
        self.__parrentpath = dirpath
        self.__ccodepath = os.path.join(self.__parrentpath,"ccode")
        self.__anspath = os.path.join(self.__parrentpath,"ans")
        self.__figpath = os.path.join(self.__parrentpath,"fig")
        if os.path.exists(self.__parrentpath):
            if os.path.exists(self.__ccodepath):
                pass
            else:
                os.mkdir(self.__ccodepath)
            if os.path.exists(self.__anspath):
                pass
            else:
                os.mkdir(self.__anspath)
            if os.path.exists(self.__figpath):
                pass
            else:
                os.mkdir(self.__figpath)
        else:
            os.mkdir(self.__parrentpath)
            os.mkdir(self.__ccodepath)
            os.mkdir(self.__anspath)
            os.mkdir(self.__figpath)
        self.__anspath = os.path.join("..","ans")
    @property
    def parrentpath(self):
        return self.__parrentpath
    @property
    def ccodepath(self):
        return self.__ccodepath
    @property
    def anspath(self):
        return self.__anspath
    @property
    def figpath(self):
        return self.__figpath

def progress_bar(iterator, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    A progress bar.
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iterator / float(total)))
    filledLength = int(length * iterator // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    if iterator == total:
        print()
