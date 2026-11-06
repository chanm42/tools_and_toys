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
function pretty_print(M,V)
##================================================================##
    atmvec = Vector{String}(undef, length(V))

    for i in eachindex(V)
        if V[i] == 1
            atmvec[i] = "Li"
        elseif V[i] == 2
            atmvec[i] = "Co"
        end
    end

    println("   ", join([@sprintf("%10s", x) for x in atmvec], " "))
    println("     ____________________________________________________________________________________________________________________________________________________________________________________")
    for (cnt, row) in enumerate(eachrow(M))

        typ = ""
        if V[cnt] == 1
            typ = "Li"
        elseif V[cnt] == 2
            typ = "Co"
        end

        rowstr = join([@sprintf("%10.3e", x) for x in row], " ")

        println(@sprintf("%2s | %s", typ, rowstr))
    end

    return nothing
end

##================================================================##
function infile_reader(infile)
##================================================================##
    atmvec = Vector{Int}(undef,16)
    atmmat = Matrix{Float64}(undef,16,3)
    open(infile, "r") do file
        for line in eachline(file)
            if contains(line,"_atom_site_occupancy")
            # hard code this part just for testing
                for iatm=1:16
                    atm  = 0
                    line = readline(file)
                    cols = split(line)
                    typ  = cols[1]
                    if typ == "Li"
                        atm = 1
                    elseif typ == "Co"
                        atm = 2
                    end
                    idx  = cols[2]
                    occ  = parse(Int,cols[3])
                    x    = parse(Float64,cols[4])
                    y    = parse(Float64,cols[5])
                    z    = parse(Float64,cols[6])
                    atmvec[iatm]   = atm
                    atmmat[iatm,1] = x
                    atmmat[iatm,2] = y
                    atmmat[iatm,3] = z
                end
            end # contains
        end
    end
    pretty_print(atmmat,atmvec)
    return atmvec,atmmat
end

##================================================================##
function calc_rpos(posvec)
##================================================================##
    @assert length(posvec) == 3 "posvec must be 3d! $posvec"
    sq = 0.0e0
    for i=1:3
        sq += posvec[i]^2
    end
    rpos = sqrt(sq)
    return rpos
end

##================================================================##
function calc_dist_sq(pos1,pos2)
##================================================================##
    ndim = length(pos1)
    distsq = 0.0e0
    for idim=1:ndim
        distsq = distsq + ((pos1[idim]-pos2[idim])^2)
    end
    return distsq
end

##================================================================##
function check_symmetry(M,V)
##================================================================##
    nrow,ncol = size(M)
    for i=1:nrow
        for j=1:nrow
            val12 = M[i,j]
            val21 = M[j,i]
            if i != j
                if val12 == val21
                    println("Symmetry in index [$i,$j] -> [$j,$i]")
                    println("M[$i,$j]($(V[i]))_M[$j,$i]($(V[j])) | $(M[i,j]) $(M[j,i])")
                end
            end
        end
    end
    return nothing
end

##================================================================##
function test1A()
##================================================================##
    coef     = 8.988e9 # coulomb coefficient in N*m^2 / C^2
    e_charge = 1.607e-19 # in C
    coef     = coef*(e_charge)^2
    infile   = "LiCo_toy_model2x2x2.cif"
    atmvec,atmmat = infile_reader(infile)
    natm = length(atmvec)
    cmat = zeros((natm,natm))
    li_q = 1
    co_q = 3
    for iatm=1:natm
        q1   = 0.0e0
        ityp = atmvec[iatm]
        ityp == 1 ? q1 = li_q : q1 = co_q
        ipos = atmmat[iatm,:]
        #println("pos_of_atm $iatm: $ipos")
        ipos = ipos .* 10.0e-10 # convert from angstrom to meters
        r1   = calc_rpos(ipos)
        for jatm=1:natm
            q2   = 0.0e0
            jtyp = atmvec[jatm]
            jtyp == 1 ? q2 = li_q : q2 = co_q
            jpos = atmmat[jatm,:]
            #println("pos_of_atm $jatm: $jpos")
            jpos = ipos .* 10.0e-10 # convert from angstrom to meters
            r2   = calc_rpos(jpos)
            r12  = calc_dist_sq(ipos,jpos)
            if iatm != jatm
                if isapprox(ipos,jpos)
                    println("*** WARNING ***")
                    println("Atm1 and Atm2 have same position but are not the same atom!")
                    println("Atm1_#Atm1|pos: $ityp $iatm | $ipos")
                    println("Atm2_#Atm2|pos: $jtyp $jatm | $jpos")
                    println()
                end

                c_eng = coef * ((q1 * q2) / (r12 + (0.5^2))) # 0.5 is a random estimate for a soft coulomb parameter
                c_eng = round(c_eng,sigdigits=3)
                cmat[iatm,jatm] = c_eng
            end # if
        end # jatm
    end # iatm
    pretty_print(cmat,atmvec)
    #check_symmetry(cmat,atmvec)
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









