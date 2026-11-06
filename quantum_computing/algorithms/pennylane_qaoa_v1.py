import scipy as sp
import matplotlib.pyplot as plt
import numpy.random as npr
import gc
import pennylane as qp
from pennylane import qaoa
from pennylane import numpy as np
from networkx import Graph
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
def pretty_print(M):
##=====================================================================================##
    nrow,ncol = M.shape
    for row in range(nrow):
        row_str = np.zeros(ncol)
        for col in range(ncol):
            x            = M[row,col]
            row_str[col] = str(x)
            row_out      = "\t".join(row_str)
        print(row_out)
    return None

##=====================================================================================##
def circuit(wires,nlayer,params,cost_h,mixer_h):
##=====================================================================================##
    for w in wires:
        qp.Hadamard(wires=w)
    for layer in range(nlayer):
        gamma = params[layer][0]
        alp   = params[layer][1]

        qaoa_layer(gamma,alp,cost_h,mixer_h)

    return None

##=====================================================================================##
def qaoa_layer(gamma,alpha,cost_h,mixer_h):
##=====================================================================================##
    qaoa.cost_layer(gamma,cost_h)
    qaoa.mixer_layer(alpha,mixer_h)
    return None

##=====================================================================================##
def cost_func(cost_h):
##=====================================================================================##
    eng = qp.expval(cost_h)
    return eng

##=====================================================================================##
def test1A():
##=====================================================================================##
    # Test function solving the max cut problem using qaoa in pennylane 
    # From https://docs.pennylane.ai/en/stable/code/qp_qaoa.html#solving-the-maxcut-problem-using-qaoa

    nlayer = 5
    params = np.full((nlayer, 2), 0.1,requires_grad=True)

    wires  = range(3)
    graph  = Graph([(0,1),(1,2),(2,0)])

    cost_h, mixer_h = qaoa.maxcut(graph)
    dev = qp.device('default.qubit',wires=len(wires))
    
    

    @qp.qnode(dev) 
    def cost_func(params):
        circuit(wires,nlayer,params,cost_h,mixer_h)
        return qp.expval(cost_h)

    opt   = qp.GradientDescentOptimizer(stepsize=0.01)
    nstep = 1000
    eng   = 0.0e0
    energies = []

    for step in range(nstep):
        params,eng = opt.step_and_cost(cost_func,params)
        #print(f"step= {step} eng= {eng} params= {params}")
        print(f"step= {step} eng= {eng}")
        energies.append(eng)

    print(f"final_params= {params} final_eng= {eng}")
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




