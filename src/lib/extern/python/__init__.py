

import libmpqc_py
from libmpqc_py import *

def __colonize__(keyword):
    r = keyword.replace('__',':')
    return r

def __keyval_assign_list__(prefix,seq,akv):
    for i in range(len(seq)):
        value = seq[i]
        prefix2 = prefix+':'+str(i)
        if type(value) == list:
            __keyval_assign_list__(prefix2,value,akv)
        elif type(value).__class__ == Molecule.__class__:
            akv.assigndescribedclass(__colonize__(kw),value)
        else:
            akv.assignstring(__colonize__(prefix2), str(value))

def __init_keyval__(self,keywords):
    keyval = AssignedKeyVal();
    keys = keywords.keys()
    for kw in keys:
        value = keywords[kw]
        if type(value) == list:
            __keyval_assign_list__(kw,value,keyval)
        elif type(value).__class__ == Molecule.__class__:
            ref_describedclass_value = value.ref()
            keyval.assigndescribedclass(__colonize__(kw),
                                                ref_describedclass_value)
        else:
            keyval.assignstring(__colonize__(kw),str(value))
    self.__init_boost__(keyval)

def __init_general__(self,*arguments,**keywords):
    if len(arguments) != 0:
        # If there are arguments, then only try the boost ctor.
        self.__init_boost__(*arguments,**keywords)
    else:
        # If there are no arguments, then first try to find a boost
        # constructor, then try the keyval constructor.
        try:
            self.__init_boost__(*arguments,**keywords)
        except:
            self.__init_keyval__(keywords)

for i in dir(libmpqc_py):
    # If this is a boost class, then make constructors specifying
    # keyword arguments use AssignedKeyVal to initialize.
    if libmpqc_py.__dict__[i].__class__ == libmpqc_py.Molecule.__class__:
       libmpqc_py.__dict__[i].__init_boost__ = libmpqc_py.__dict__[i].__init__
       libmpqc_py.__dict__[i].__init_keyval__ = __init_keyval__
       libmpqc_py.__dict__[i].__init__ = __init_general__

for i in [SCMatrix, SymmSCMatrix, DiagSCMatrix, SCVector]:
    i.__add__ = i.operator_add
    i.__mul__ = i.operator_mul
