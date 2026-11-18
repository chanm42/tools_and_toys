import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.random as npr
import gc
import math
from   datetime import datetime

#USE THIS FOR FUNCTION SEPARATORS
##=====================================================================================##

#USE FOR SECTION SEPARATORS
#----

##=====================================================================================##
def debug_exit(mesg):
##=====================================================================================##
    print("*** DEBUG EXIT ***")
    print(mesg)
    exit()

    return None

##=====================================================================================##
def vec_to_str(v):
##=====================================================================================##
    str_v = str(v[0])

    for i in range(2,len(v)):
        str_v = str_v + " " + str(v[i])

    return str_v

##=====================================================================================##
def freemem(v):
##=====================================================================================##
    v = None
    gc.collect()
    return v


##=====================================================================================##
def infile_reader(infile):
##=====================================================================================##
    xvec = []
    yvec = []
    with open(infile,"r") as f:
        for line in f:
            cols = line.split(line)
            xvec.append(float(cols[0]))
            yvec.append(float(cols[1]))
    return xvec,yvec

##=====================================================================================##
def pretty_print(M,precision=4):
##=====================================================================================##
    nrow,ncol = M.shape
    for row in range(nrow):
        row_str = []
        for col in range(ncol):
            x            = M[row,col]
            row_str.append(f"{x:12.{precision}f}")
            #assert isinstance(row_str,str), "Argument of wrong type"
        row_out      = " ".join(row_str)
        print(row_out)
    return None

##=====================================================================================##
def get_2x2_evecs_evals(M):
##=====================================================================================##
    nrow,ncol = M.shape
    assert nrow == 2
    assert nrow == ncol
    
    i = M[1,1]
    j = M[2,2]
    k = M[2,1]
    l = M[1,2]

    a = 1
    b = i + j
    c = i*j - l*k

    eval1 = -b + (np.sqrt(1 - 4*a*c))/(2*a)
    eval2 = -b - (np.sqrt(1 - 4*a*c))/(2*a)
    evals = [eval1,eval2]

    vec_root1 = i - j + np.sqrt((i-j)^2 + 4*k*l) 
    vec_root2 = i - j - np.sqrt((i-j)^2 + 4*k*l) 
   
    evec1 = [vec_root1,2*l]
    evec2 = [vec_root,2*l]

    evecs = [evec1,evec2]

    return evecs, evals 

##=====================================================================================##
def make_1d_tb_hamil(nsites,diag,off_diag):
##=====================================================================================##
    H = np.zeros((nsites,nsites))
    for i in range(0,nsites):
        for j in range(0,nsites):
            if i == j:
                H[i,j] = diag
            elif abs(i-j) == 1:
                H[i,j] = off_diag
    return H

##=====================================================================================##
def test1A():
##=====================================================================================##
    nsites   = 8
    diag     = 2
    off_diag = -1
    H        = make_1d_tb_hamil(nsites,diag,off_diag)
    print(">> TB Hamiltonian initialized:")
    pretty_print(H)
    ngroup   = 2 # set the size of the block for RG method
    nblock   = math.floor(nsites / ngroup)
    for iblock in range(0,nblock):
        idx1 = 2*iblock
        idx2 = 2*iblock+1
        sub_mat     = H[idx1:idx2,idx1:idx2]
        evecs,evals = get_2x2_evecs_evals(sub_mat)

    return None



##=====================================================================================##
## MAIN CODE ##
##=====================================================================================##
if __name__ == "__main__":
#----
    print()
    print(f"Starting main from file: {__file__}")
    print("---------------------------------------------")
    print(f"Start time: {datetime.now()}")
    print("---------------------------------------------")
#----
    test1A()    
#----
    print()
    print("---------------------------------------------")
    print(f"End time: {datetime.now()}")
    print("---------------------------------------------")




