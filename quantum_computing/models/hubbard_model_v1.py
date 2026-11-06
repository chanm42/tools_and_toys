from unittest import result

import numpy as np
import scipy as sp
import numpy.random as npr
import gc
from qiskit_nature.second_q.hamiltonians.lattices import (
        LineLattice,
        BoundaryCondition,
    )

from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from datetime import datetime
#from . import quantum_gates as qg


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
        row_str = []
        for col in range(ncol):
            row_str.append(str(M[row,col]))
            row_out = "\t".join(row_str)
        print(row_out)
    return None

##=====================================================================================##
def print_state(state,norb):
##=====================================================================================##
    binary_state = format(state,f'0{norb}b')
    print("#----")
    print(f"Integer state: {state}")
    print(f"Binary state: {binary_state}")
    print("#----")

    orbital = 0
    while orbital < norb:
        occupation = occ(state,orbital)
        print(f"Orbital: {orbital}  Occupation: {occupation}")
        orbital += 1

    print("#----")
    return None

##=====================================================================================##
def print_wavefunction(psi,norb,tol=1.0e-8):
##=====================================================================================##
    print("#----")
    print("Wavefunction components:")
    print("#----")
    state = 0
    while state < len(psi):
        amp  = psi[state]
        prob = np.abs(amp)**2
        if abs(amp) > tol:
            binary_state = format(state,f'0{norb}b')
            print(f"Basis state: | {binary_state}")
            print(f"Amplitude: {amp:.4f}")
            print(f"Probability: {prob:.4f}")
            print("#----")
            orbital = 0
            while orbital < norb:
                occupation = occ(state,orbital)
                print(f"Orbital: {orbital}  Occupation: {occupation}")
                orbital += 1
        state += 1

    print("#----")
    return None

##=====================================================================================##
def check_normalization(state):
##=====================================================================================##
    print("#----")
    norm = np.sum(np.abs(state)**2)
    print(f"Normalization: {norm:.6f}")
    print("#----")

    return None

##=====================================================================================##
def x_gate():
##=====================================================================================##
    return np.array([[0, 1], [1, 0]],dtype=complex)

##=====================================================================================##
def y_gate():
##=====================================================================================##
    return np.array([[0, -1j], [1j, 0]],dtype=complex)

##=====================================================================================##
def z_gate():
##=====================================================================================##
    return np.array([[1, 0], [0, -1]],dtype=complex)

##=====================================================================================##
def cnot_gate():
##=====================================================================================##
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 0, 1],
                     [0, 0, 1, 0]],dtype=complex)

##=====================================================================================##
def rotation_y(theta):
##=====================================================================================##
    return np.array(
        [[np.cos(theta/2), -np.sin(theta/2)],
         [np.sin(theta/2), np.cos(theta/2)]],dtype=complex
    )

##=====================================================================================##
def apply_single_qubit_gate(nqubit, qubit, gate):
##=====================================================================================##
    assert qubit < nqubit, "Qubit index out of range"
    op = 1
    I  = np.eye(2,dtype=complex)
    for i in range(nqubit):
        if i == qubit:
            op = np.kron(op,gate)
        else:
            op = np.kron(op,I)
    return op

##=====================================================================================##
def apply_two_qubit_gate(nqubit, qubit1, qubit2, gate):
##=====================================================================================##
    assert qubit1 < nqubit, "Qubit1 index out of range"
    assert qubit2 < nqubit, "Qubit2 index out of range"
    assert qubit1 != qubit2, "Qubit1 and Qubit2 must be different"
    
    # Apply two qubit gates on adjacent sites for Hubbard model hopping
    # Remove this condition for general Hamiltonians
    if qubit2 not in [qubit1-1, qubit1+1]:
        raise ValueError("Qubit1 and Qubit2 must be adjacent")
    
    op = 1
    q  = 0
    I  = np.eye(2,dtype=complex)
    while q < nqubit:
        if q == qubit1:
            op = np.kron(op,gate)
            q  += 2
        else:
            op = np.kron(op,I)
            q += 1
    return op

##=====================================================================================##
def occ(state,orbital):
##=====================================================================================##
    # Shift the bitstring so that the indicated state is on the far left
    # Then read out 0 or 1
    shift_state = state >> orbital
    orb_occ     = shift_state & 1
    return orb_occ

##=====================================================================================##
def count_left_occ(state,orb):
##=====================================================================================##
    # Count the number of occupied states to the left of the indicated state in the bitstring
    count    = 0
    for left_orb in range(orb):
        occup = occ(state,left_orb)
        
        if occup == 1:
            count += 1
        
    return count

##=====================================================================================##
def parity_sign(state,orbital):
##=====================================================================================##
    nleft  = count_left_occ(state,orbital)
    parity = 0

    if nleft % 2 == 0:
        parity = 1

    else:
        parity = -1

    return parity

##=====================================================================================##
def create(state,orbital):
##=====================================================================================##
    new_state = None
    sign      = 0
    # Can't create if already occupied
    if occ(state,orbital) == 1:
        new_state = None
        sign      = 0
        return new_state, sign

    sign      = parity_sign(state,orbital)
    new_state = state | (1 << orbital)

    return new_state, sign

##=====================================================================================##
def majorana_plus(j,nqubit):
##=====================================================================================##
    c  = creation(j,nqubit)
    a  = annihilation(j,nqubit)
    mp = c + a
    return mp

##=====================================================================================##
def majorana_minus(j,nqubit):
##=====================================================================================##
    c  = creation(j,nqubit)
    a  = annihilation(j,nqubit)
    mm = -1j * (c - a)
    return mm

##=====================================================================================##
def majorana(k,nqubit):
##=====================================================================================##
    orb = k // 2
    maj_op = None
    if k%2 == 0:
        maj_op = majorana_plus(orb,nqubit)
    else:
        maj_op = majorana_minus(orb,nqubit)

    return maj_op

##=====================================================================================##
def number_majorana(orb,nqubit):
##=====================================================================================##
    I = np.eye(2**nqubit,dtype=complex)
    num = 0.5 * (I - 1j * majorana(2*orb,nqubit) @ majorana((2*orb)+1,nqubit))
    return num

##=====================================================================================##
def annihilate(state,orbital):
##=====================================================================================##
    new_state = None
    sign      = 0
    # Can't annihilate if already empty
    if occ(state,orbital) == 0:
        new_state = None
        sign      = 0
        return new_state, sign

    sign      = parity_sign(state,orbital)
    mask      = ~(1 << orbital)
    new_state = state & mask

    return new_state, sign

##=====================================================================================##
def creation(j,nqubit):
##=====================================================================================##
    a_dag = 1
    I     = np.eye(2,dtype=complex)
    for i in range(nqubit):
        # Note matrix a is an argument of np.kron - creating strings of Pauli Z's + a single combination of x and y gates
        if i < j:
            a_dag = np.kron(a_dag,z_gate())
        elif i == j:
            a_dag = np.kron(a_dag,((0.5)*(x_gate() - (1j * y_gate()))))
        else:
            a_dag = np.kron(a_dag,I)

    return a_dag

##=====================================================================================##
def annihilation(j,nqubit):
##=====================================================================================##
    a = 1
    I = np.eye(2,dtype=complex)
    for i in range(nqubit):
        # Note matrix a is an argument of np.kron - creating strings of Pauli Z's + a single combination of x and y gates
        if i < j:
            a = np.kron(a,z_gate())
        elif i == j:
            a = np.kron(a,((0.5)*(x_gate() + (1j * y_gate()))))
        else:
            a = np.kron(a,I)

    return a

##=====================================================================================##
def number(j,nqubit):
##=====================================================================================##
    num = creation(j,nqubit) @ annihilation(j,nqubit)
    return num


##=====================================================================================##
def single_hop(H,state,orb_i,orb_j,t):
##=====================================================================================##
    int_state, sign1 = annihilate(state,orb_j)
    
    if int_state is None:
        return H
    
    fin_state, sign2 = create(int_state,orb_i)

    if fin_state is None:
        return H
    
    tot_sign = sign1 * sign2
    mat_ele  = -t * tot_sign
    # <final|H|initial>
    H[fin_state,state] += mat_ele

    return H

##=====================================================================================##
def build_interaction_matrix(H,num_sites,U):
##=====================================================================================##
    assert H.shape[0] == H.shape[1], "H must be square"
    
    dimH = H.shape[0]

    state = 0
    for state in range(dimH):
        int_eng = 0.0
        site    = 0
        for site in range(num_sites):
            up_orb   = orb(site,0,num_sites)
            down_orb = orb(site,1,num_sites)

            nup     = occ(state,up_orb)
            ndown   = occ(state,down_orb)
            dub_occ = nup * ndown

            int_eng += (U * dub_occ)
            # end for site
        
        H[state,state] = int_eng
        # end for state
    
    return H

##=====================================================================================##
def build_hopping_matrix(H,num_sites,t):
##=====================================================================================##
    assert H.shape[0] == H.shape[1], "H must be square"
    
    dimH = H.shape[0]

    for state in range(dimH):
        
        for site in range(num_sites-1):
            up_i = orb(site,0,num_sites)
            up_j = orb(site+1,0,num_sites)

            H = single_hop(H,state,up_i,up_j,t)
            H = single_hop(H,state,up_j,up_i,t)

            down_i = orb(site,1,num_sites)
            down_j = orb(site+1,1,num_sites)

            H = single_hop(H,state,down_i,down_j,t)
            H = single_hop(H,state,down_j,down_i,t)
            
            # end for site
        
        # end for state
    
    return H

##=====================================================================================##
def majorana_hopping_matrix(i,j,nqubit,t):
##=====================================================================================##
    g2i   = majorana(2*i,nqubit)
    g2ip1 = majorana((2*i)+1,nqubit)
    g2j   = majorana(2*j,nqubit)
    g2jp1 = majorana((2*j)+1,nqubit)

    mat_element = 0.5j*t*(g2i@g2jp1 - g2ip1@g2j)
    return mat_element

def majorana_interaction(orb,nqubit):
    return None

##=====================================================================================##
def orb(site,spin,num_sites):
##=====================================================================================##
    orb = None
    if spin == 0:
        orb = site
    else:
        orb = site + num_sites
    assert orb is not None, "Invalid spin value"
    return orb

##=====================================================================================##
def neighbors(site,num_sites,boundary="open"):
##=====================================================================================##
    if boundary == "open":
        if site > 0:
            yield site-1
        if site < num_sites-1:
            yield site+1

##=====================================================================================##
def build_qubit_matrix(num_orb,num_sites,t,U):
##=====================================================================================##
    dimH = 2**num_orb
    H    = np.zeros((dimH,dimH),dtype=complex)

    #---------------------------
    # Interaction term
    #---------------------------
    for site in range(num_sites):

        up = orb(site,0,num_sites)
        dn = orb(site,1,num_sites)

        H += U * (
            number(up,num_orb)
            @
            number(dn,num_orb)
        )

    #---------------------------
    # Hopping term
    #---------------------------
    for site in range(num_sites-1):

        for spin in range(2):

            orb_i = orb(site,spin,num_sites)
            orb_j = orb(site+1,spin,num_sites)

            a_dag_i = creation(orb_i,num_orb)
            a_j     = annihilation(orb_j,num_orb)

            a_dag_j = creation(orb_j,num_orb)
            a_i     = annihilation(orb_i,num_orb)

            H += -t * (
                (a_dag_i @ a_j)
                +
                (a_dag_j @ a_i)
            )

    return H

##=====================================================================================##
def build_qubit_matrix_majorana(num_orb,num_sites,t,U):
##=====================================================================================##
    dimH = 2**num_orb
    H    = np.zeros((dimH,dimH),dtype=complex)

    #---------------------------
    # Interaction term
    #---------------------------
    for site in range(num_sites):

        up = orb(site,0,num_sites)
        dn = orb(site,1,num_sites)

        H += U * (
            number_majorana(up,num_orb)
            @
            number_majorana(dn,num_orb)
        )

    #---------------------------
    # Hopping term
    #---------------------------
    for site in range(num_sites-1):

        for spin in range(2):
            orb_i = orb(site,spin,num_sites)
            orb_j = orb(site+1,spin,num_sites)
            Hop2  = majorana_hopping_matrix(orb_i,orb_j,num_orb,t)
            H    += Hop2
    return H

##=====================================================================================##
def build_ansatz(params,nqubits):
##=====================================================================================##
    dim = 2**nqubits
    psi = np.zeros(dim,dtype=complex)
    psi[0] = 1.0

    for qubit in range(nqubits):
        gate = rotation_y(params[qubit])
        ent  = cnot_gate()
        U    = apply_single_qubit_gate(nqubits,qubit,gate)
        psi  = U @ psi

    for qubit in range(nqubits-1):
        ent = cnot_gate()
        #print(f"Qubit1: {qubit}  Qubit2: {qubit+1}")
        U   = apply_two_qubit_gate(nqubits,qubit,qubit+1,ent)
        psi = U @ psi
    return psi

##=====================================================================================##
def eng_expectation(params,H,nqubits):
##=====================================================================================##
    psi = build_ansatz(params,nqubits)
    eng = np.vdot(psi,H @ psi)
    eng = np.real(eng)
    return eng

##=====================================================================================##
def test1A():
##=====================================================================================##

    norb = 4

    state = 15

    print()
    print("INITIAL STATE")
    print("NOTE: Orbital 0 is the rightmost bit")
    print_state(state,norb)

    orbital = 1

    print()
    print(f"Applying annihilation operator at orbital {orbital}")

    new_state,sign = annihilate(state,orbital)

    print(f"Fermionic sign = {sign}")

    if new_state is not None:

        print()
        print("FINAL STATE")
        print_state(new_state,norb)

    return None

##=====================================================================================##
def test2A(t,U):
##=====================================================================================##
    num_sites = 2  
    nspin     = 2

    norb = num_sites * nspin
    dimH = 2**norb

    H = np.zeros((dimH,dimH))
    H = build_interaction_matrix(H,num_sites,U)
    H = build_hopping_matrix(H,num_sites,t)

    Hdiff = H - H.T
    maxdiff = np.max(np.abs(Hdiff))
    print()
    print(">>Beginning exact calculation<<")
    print(f"Max difference between H and H^T: {maxdiff:.2e}")

    evals,evecs = np.linalg.eigh(H)
    print()
    print("Eigenvalues:")
    for i, eval in enumerate(evals):
        print(f"  {i}: {eval:.6f}")
    print()
    print("Ground state eigenvector:")
    ground_state = evecs[:,0]
    print(ground_state)
    print_wavefunction(ground_state,norb)
    pretty_print(H)
    #print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in H]))
    check_normalization(ground_state)

    for i in range(dimH):
        for j in range(dimH):
            diff = H[i,j] - H[j,i]
            if abs(diff) > 1e-12:
                print("#----")
                print(f"ASYMMETRY DETECTED")
                print(f"H[{i},{j}] = {H[i,j]}")
                print(f"H[{j},{i}] = {H[j,i]}")
                print("#----")

    return None

##=====================================================================================##
def test3A(t,U):
##=====================================================================================##
    num_runs  = 10
    num_sites = 2
    nspin     = 2
    nqubits   = num_sites * nspin

    norb = num_sites * nspin
    dimH = 2**norb

    H    = build_qubit_matrix(norb,num_sites,t,U)
    Hmaj = build_qubit_matrix_majorana(norb,num_sites,t,U)

    Hdiff      = H - H.T
    maxdiff    = np.max(np.abs(Hdiff))
    Hmajdiff   = Hmaj - Hmaj.T
    maxdiffmaj = np.max(np.abs(Hmajdiff))
    diff       = np.max(np.abs(H - Hmaj))
    print()
    print(">>Beginning qubit matrix calculation<<")
    print(f"Max difference between H and H^T: {maxdiff:.2e}")
    print(f"Max difference between H_maj and H_maj^T: {maxdiffmaj:.2e}")
    print(f"Max difference between H_maj and H: {diff:.2e}")
    evals,evecs = np.linalg.eigh(H)
    majevals,majevecs = np.linalg.eigh(Hmaj)
    print()
    print("Standard qubit representation:")
    print("Eigenvalues:")
    for i, eval in enumerate(evals):
        print(f"  {i}: {eval:.6f}")
    print()
    print("Ground state eigenvector:")
    ground_state = evecs[:,0]
    print_wavefunction(ground_state,norb)
    #pretty_print(H)
    print()
    print("Majorana qubit representation:")
    print("Eigenvalues:")
    for i, eval in enumerate(majevals):
        print(f"  {i}: {eval:.6f}")
    print()
    print("Ground state eigenvector:")
    ground_state = majevecs[:,0]
    print_wavefunction(ground_state,norb)
    #pretty_print(Hmaj)
    #print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in H]))
    check_normalization(ground_state)

    for i in range(dimH):
        for j in range(dimH):
            diff = H[i,j] - H[j,i]
            if abs(diff) > 1e-12:
                print("#----")
                print(f"ASYMMETRY DETECTED")
                print(f"H[{i},{j}] = {H[i,j]}")
                print(f"H[{j},{i}] = {H[j,i]}")
                print("#----")
    print()
    print("#----")
    print(">>Beginning VQE optimization<<")
    print(f"Using COBYLA optimizer and a rotation + CNOT entangler ansatz for {num_runs} runs with random initial parameters")

    res_vec   = []
    param_vec = []
    print("For normal Hamiltonian:")
    for i in range(num_runs):
        params0  = np.random.rand(nqubits)
        result   = sp.optimize.minimize(
        eng_expectation,
        params0,
        args=(H,nqubits),
        method="COBYLA"
        )
        res_vec.append(result.fun)
        param_vec.append(result.x)
        print(f"Run {i+1}: Optimal energy = {result.fun:.6f}  Optimal parameters = {result.x}")
    print("Average Optimal parameters:")
    print(np.mean(param_vec,axis=0))
    print("Average Optimal energy:")
    print(np.mean(res_vec))
    print("Standard deviation of optimal energy:")
    print(np.std(res_vec))

    print("Exact ground state energy:")
    print(evals[0])
    print("#----")


    res_vec   = []
    param_vec = []
    print("For Majorana Hamiltonian:")
    for i in range(num_runs):
        params0  = np.random.rand(nqubits)
        result   = sp.optimize.minimize(
        eng_expectation,
        params0,
        args=(H,nqubits),
        method="COBYLA"
        )
        res_vec.append(result.fun)
        param_vec.append(result.x)
        print(f"Run {i+1}: Optimal energy = {result.fun:.6f}  Optimal parameters = {result.x}")
    print("Average Optimal parameters:")
    print(np.mean(param_vec,axis=0))
    print("Average Optimal energy:")
    print(np.mean(res_vec))
    print("Standard deviation of optimal energy:")
    print(np.std(res_vec))

    print("Exact ground state energy:")
    print(evals[0])
    print("#----")

    return None 

##=====================================================================================##
def test4A(t,U):
##=====================================================================================##
    num_sites = 2

    #----
    # Build lattice
    #----

    lattice = LineLattice(
        num_nodes=num_sites,
        boundary_condition=BoundaryCondition.OPEN,
    )

    #----
    # Assign hopping parameters
    #----

    lattice = lattice.uniform_parameters(
        uniform_interaction=-t,
        uniform_onsite_potential=0.0,
    )

    #----
    # Build Hubbard Hamiltonian
    #----

    model = FermiHubbardModel(
        lattice,
        onsite_interaction=U,
    )

    fermionic_op = model.second_q_op()

    print("\nFermionic Hamiltonian:")
    print(fermionic_op)


    mapper = JordanWignerMapper()

    qubit_op = mapper.map(fermionic_op)

    print("\nQubit Hamiltonian:")
    print(qubit_op)

    H = qubit_op.to_matrix()

    print("\nHamiltonian matrix:")
    pretty_print(H)

    evals, evecs = np.linalg.eigh(H)

    print("\nEigenvalues:")
    for i, eval in enumerate(evals):
        print(f"  {i}: {eval:.6f}")
    print()
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
    t = 1.0
    U = 2.0
    #test1A()  
    #test2A(t,U)  
    test3A(t,U)
    #test4A(t,U)
#----
    print()
    print("---------------------------------------------")
    print(f"End time: {datetime.now()}")
    print("---------------------------------------------")




