using Random
using LinearAlgebra
using SpecialFunctions
using Dates
using Base: Float64, Integer
using Distributions
using Printf
using Parameters

#USE THIS FOR SECTION SEPARATOR: 
  #---

#USE THIS FOR FUNCTION MARKERS: 
##================================================================##

@with_kw struct Shell
    ang_mom_m ::Int
    alpvec_m  ::Vector{Float64}
    coefvec_m ::Vector{Float64}
    repr_m    ::Int
end

@with_kw struct Atom
    type_m  ::String
    center_m::Vector{Float64}
    shells_m::Vector{Shell}
end

##================================================================##
function api_build_atom(type,center,shells)
##================================================================##
    this_atom = Atom(type_m   = type,
                     center_m = center,
                     shells_m = shells
                    )
    return this_atom
end

##================================================================##
function api_build_shell(l,alpvec,coefvec,repr)
##================================================================##
    @assert length(alpvec) == length(coefvec) "Alpha and coef vectors aren't same length! $alpvec | $coefvec"
    this_shell = Shell(ang_mom_m = l,
                       alpvec_m  = alpvec,
                       coefvec_m = coefvec,
                       repr_m    = repr
                      )
    return this_shell
end

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
function get_ang_mom(type)
##================================================================##
    l::Int    = 0
    # use 0 for cartesian representation and 1 for spherical
    repr::Int = 1
    if type == "S"
        l = 0
    elseif type == "P"
        l = 1
    elseif type == "D"
        l = 2
    elseif type == "F"
        l = 3
    elseif type == "G"
        l = 4
    elseif type == "H"
        l = 5
    elseif type == "I"
        l = 6
    end
    return l,repr
end

##================================================================##
function update_nbas(nbas,l,repr)
##================================================================##
    """
    Takes in current number of basis function nbas, angular momentum (orbital type) l, and representation (cartesian vs spherical) to update nbas
    """
    if repr == 1
    # spherical harmonic basis functions
        nbas += (2*l + 1)
    elseif repr == 0
    # Cartesian Gaussian basis functions
        nbas += ((l+1)(l+2)) / 2
    end
    @assert typeof(nbas) == Int "Number of basis functions must be an integer!"
    return nbas
end

##================================================================##
function read_ebasis(infile)
##================================================================##
    atmlist      = []
    centcnt::Int = 0
    
    open(infile, "r") do file
        header  = readline(file)
        ncent   = readline(file)
        ncent   = parse(Int,ncent)
        # create the atom list
        atmlist = Vector{Atom}(undef,ncent)
        # pre-allocate variables for the atom struct
        atmtyp     = ""
        cent       = zeros(3)
        shell_list = []
        for line in eachline(file)
            line = strip(line)
            if occursin("---------------",line)
                cols    = split(line)
                atmtyp  = cols[2]
            end # atmtyp
            if line == "CENTER"
                centcnt += 1
                line = readline(file)
                cols = split(line)
                @assert length(cols) == 3 "Center vector is not in 3d! $line | $cols"
                x,y,z = parse.(Float64,cols)
                cent[1] = x
                cent[2] = y
                cent[3] = z
            end # center
            if line == "SYMBOLS"
                orbl   = readline(file)
                orbl   = strip(orbl)
                nshell = parse(Int,orbl)
                # pre-allocate the shell list variable for the atom struct
                shell_list = Vector{Shell}(undef,nshell)
                for ishell=1:nshell
                    line         = readline(file)
                    shelltyp,npg = split(line)
                    npg          = parse(Int,npg)
                    # pre-allocate the alpvec and coefvec variable for the shell struct
                    alpvec     = zeros(npg)
                    coefvec    = zeros(npg)
                    for pg=1:npg
                        line = readline(file)
                        idx,alp,coef = split(line)
                        idx   = parse(Int,idx)
                        alp   = parse(Float64,alp)
                        coef  = parse(Float64,coef)
                        alpvec[pg]  = alp
                        coefvec[pg] = coef
                    end # primitives
                    # get angular momentum value and representation (Cartesian/Spherical) for shell struct
                    l,repr = get_ang_mom(shelltyp)
                    shell  = api_build_shell(l,alpvec,coefvec,repr)
                    shell_list[ishell] = shell
                end # shell loop
                # create the atom struct
                atmcent   = deepcopy(cent)
                atmshells = deepcopy(shell_list)
                atom      = api_build_atom(atmtyp,atmcent,atmshells)
                atmlist[centcnt] = atom
            end # symbols
        end
    end
    return atmlist
end

##================================================================##
function read_mo_coef(infile)
##================================================================##
    mo_coef_mat = nothing
    open(infile,"r") do file
        header  = readline(file)
        line    = readline(file)
        nao,nmo = parse.(Int,split(line))
        header  = readline(file)
        mo_coef_mat = zeros(nao,nmo)
        for line in eachline(file)
            cols = split(line)
            iao  = parse(Int,cols[1])
            imo  = parse(Int,cols[2])
            @assert 1 <= iao <= nao "iao outside of nao range!"
            @assert 1 <= imo <= nmo "imo outside of nmo range!"
            coef = parse(Float64,cols[3])
            mo_coef_mat[iao,imo] = coef
        end
    end
    return mo_coef_mat
end

##================================================================##
function read_grid_info(infile)
##================================================================##
    ndim     = 0
    grid_inf = nothing
    open(infile,"r") do file
        header   = readline(file)
        ndim     = parse(Int,readline(file))
        header   = readline(file)
        grid_inf = zeros(ndim,3) # ipt | imin | imax
        for idim=1:ndim
            line = readline(file)
            cols = split(line)
            npt  = parse(Int,cols[1])
            imin = parse(Float64,cols[2])
            imax = parse(Float64,cols[3])
            grid_inf[idim,1] = ndim
            grid_inf[idim,2] = imin
            grid_inf[idim,3] = imax
        end
    end
    return grid_inf
end

##================================================================##
function build_grid(grid_inf)
##================================================================##
    ndim,ninf = size(grid_inf)
    npt_l     = zeros(ndim)
    rangemat  = zeros(ndim,2)
    for idim=1:ndim
        npt  = grid_inf[idim,1]
        imin = grid_inf[idim,2]
        imax = grid_inf[idim,3]
        npt_l[idim] = npt
        rangemat[idim,1] = imin
        rangemat[idim,2] = imax
    end
    grid = zeros(Float64,npt,npt,npt)
    xr   = range(rangemat[1,1],rangemat[1,2],length=npt_l[1])
    yr   = range(rangemat[2,1],rangemat[2,2],length=npt_l[2])
    zr   = range(rangemat[3,1],rangemat[3,2],length=npt_l[3])
    cnt_x = 0
    cnt_y = 0
    cnt_z = 0
    for x in xr
        cnt_x += 1
        for y in yr
            cnt_y += 1
            for z in zr
                cnt_z += 1
                ans = build_mo_at_gridpoint()
                grid[cnt_x,cnt_y,cnt_z] = build_mo 
            end # zr
        end # yr
    end # xr
    return grid
end

##================================================================##
function construct_basis_funcs()
##================================================================##
    infile_ebas  = "./ebasis.inp"
    infile_mo    = "./mo_to_ao.inp"
    infile_grid  = "./grid_maker_info.inp"
    atmlist  = read_ebasis(infile_ebas)
    mo_mat   = read_mo_coef(infile_mo)
    grid_inf = read_grid_info(infile_grid)
    #grid     = build_grid(grid_inf)
    grid     = zeros(Float64,npt,npt,npt)
    natm     = length(atmlist)
    ndim,inf = size(grid_inf)
    for iatm=1:natm
        atom   = atmlist[iatm]
        typ    = atom.type_m
        cent   = atom.cent_m
        shells = atom.shells_m
        nshell = length(shells)
        for ishell=1:nshell
            shell   = shells[ishell]
            ang_mom = shell.ang_mom_m
            repr    = shell.repr_m
            alpvec  = shell.alpvec_m
            coefvec = shell.coefvec_m
            nprim   = length(alpvec)
            for iprim=1:nprim
                alp  = alpvec[iprim]
                coef = coefvec[iprim]
                for idim=1:ndim
                    npt = grid_inf[idim,1]
                    min = grid_inf[idim,2]
                    max = grid_inf[idim,3]
                    loc = cent[idim]
                    cnt_x = 0
                    cnt_y = 0
                    cnt_z = 0
                    for x in range(min,max,length=npt)
                        cnt_x += 1
                        for y in range(min,max,length=npt)
                            cnt_y += 1
                            for z in range(min,max,length=npt)
                                cnt_z += 1
                                gaus = coef * exp(-alp * (x - loc)^2)
                            end # z
                        end # y
                    end # x
                end # ndim
            end # prim
        end # shell
    end # atom
    return nothing
end

##================================================================##
function test1A()
##================================================================##
    ebasis_infile = "./ebasis.inp"
    atmlist = read_ebasis(ebasis_infile)
    # test to check ebasis reader
    natm = length(atmlist)
    nbas = 0
    ngas = 0
    println(">Total atoms in the list: $natm")
    println(">===============================================<")
    for iatm=1:natm
        atom   = atmlist[iatm]
        type   = atom.type_m
        cent   = atom.center_m
        shells = atom.shells_m
        nshell = length(shells)
        println(">>Type for atom $iatm: $type")
        println(">>Atom $iatm centered at $cent")
        println(">>Atom $iatm has $nshell shells")
            for ishell=1:nshell
            shell   = shells[ishell]
            ang_mom = shell.ang_mom_m
            alpvec  = shell.alpvec_m
            coefvec = shell.coefvec_m
            repr    = shell.repr_m
            nbas    = update_nbas(nbas,ang_mom,repr)
            nprim   = length(alpvec)
            @assert length(coefvec) == nprim "Different numbers of alpha values and coefficients! $alpvec | $coefvec"
            println(">>>Angular momentum for shell $ishell : $ang_mom")
            println(">>>Shell $ishell has $nprim primitives:")
            for iprim=1:nprim
                ngas += 1
                alp   = alpvec[iprim]
                coef  = coefvec[iprim]
                println(">>>>Primitive $iprim : alp=$alp coef=$coef")
            end
        end
    end
    println(">===============================================<")
    println(">Total number of basis functions: $nbas")
    println(">Total number of Gaussians to construct: $ngas")
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









