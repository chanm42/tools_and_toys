using Random
using LinearAlgebra
using SpecialFunctions
using Dates
using Base: Float64, Integer
using Distributions
using Printf

#USE THIS FOR SECTION SEPARATOR: 
  #---

#USE THIS FOR FUNCTION MARKERS: 
##================================================================##


##================================================================##
function debug_exit(mesg)
##================================================================##
    println("****DEBUG_EXIT****")
    println(mesg)
    exit()
end
##================================================================##
function vec_to_str(v)
##================================================================##
  #---
    str_v = string(v[1])
  #---
    for i=2:length(v)
        str_v = str_v * " " * string(v[i])
    end
    return str_v
end
##================================================================##
function freemem(v)
##================================================================##
    v = nothing
    GC.gc()
    return v
end

##================================================================##
function infile_reader(infile)
##================================================================##
    x_vec = Float64[]
    y_vec = Float64[]

    open(infile, "r") do file
        for line in eachline(file)
            cols = split(line)
            x_val = push!(x_vec, parse(Float64, cols[1]))
            y_val = push!(y_vec, parse(Float64, cols[2]))
        end
    end
    return x_vec, y_vec
end

##================================================================##
function build_sym_mat(n)
##================================================================##
    println(">> Building symmetric A matrix")
    A = rand(Int8,n,n) 
    A = (A * A') / 2
    println(">> A constructed and symmetrized")
    if true
        print(A)
    end
    return A
end

##================================================================##
function svd_v1(A) # compute SVD using the most basic implementation of M'*M and M*M'
##================================================================##
    println(">> Computing the SVD of A using v1")
    # Compute AA' and A'A
    A1 = A*A'
    A2 = A'*A
    # Get eigenvalues and eigenvectors of A1 and A2
    eval1,evec1 = eigen(A1)
    eval2,evec2 = eigen(A2)
    nleft       = length(eval2)
    nright      = length(eval1)

    # Construct U and V' -> Columns ov V are evec1 and columns of U are evec2
    U = stack(evec2)
    V = stack(evec1)

    svals = eigvals(A1)
    svals = sqrt.(svals)
    S     = Matrix{Float64}(I,size(A))
    S     = S .* svals

    return U,S,V
end

##================================================================##
function test1A()
##================================================================##
    n::Int = 10

    A     = build_sym_mat(n)
    U,S,V = svd_v1(A)

    # Do some verification
    eval,evec = eigen(A)
    A_new     = U*S*V'
    eval_new,evec_new = eigen(A_new)
    n = length(eval)
    for i=1:n
        e1 = eval[i]
        e1 = sqrt(e1^2)
        e2 = eval_new[i]
        println(">>> Original_eval_new_eval: $e1 | $e2")
    end

    
    return nothing
end

##================================================================##
## MAIN CODE ##
##================================================================##
if abspath(PROGRAM_FILE) == @__FILE__
  #---
    println()
    println("Starting main from file: ", PROGRAM_FILE)
    println("--------------------------------------------")
    println("start_time= ",now())
    println("--------------------------------------------")
    println()
  #---
        test1A()
  #---
    println()
    println("--------------------------------------------")
    println("end_time= ",now())
    println("--------------------------------------------")
  #---
end #abspath









