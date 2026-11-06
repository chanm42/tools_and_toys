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
function infile_checker(infile1,infile2)
##================================================================##
    lines1 = readlines(infile1)
    lines2 = readlines(infile2)

    cl1    = length(lines1)
    cl2    = length(lines2)
    if cl1 != cl2
        println("$infile1 and $infile2 have different line counts!")
        println("$infile1 linecount= $cl1")
        println("$infile2 linecount= $cl2")
        return nothing
    end
    println("$infile1 linecount= $cl1")
    println("$infile2 linecount= $cl2")

    for lc=1:cl1
        line1 = lines1[lc]
        line2 = lines2[lc]

        line1 = split(strip(line1))
        line2 = split(strip(line2))

        len1  = length(line1)
        len2  = length(line2)
        if len1 != len2
            println("line $lc in $infile1 and $infile2 are not the same length!")
            println("$infile1 = $line1")
            println("$infile2 = $line2")
            return nothing
        end

        for i=1:len1
            li1 = line1[i]
            li2 = line2[i]

            if li1 != li2
                println("Difference found at line $lc element $i in $infile1 and $infile2")
                println("$infile1 - $li1 -- $line1")
                println("$infile2 - $li2 -- $line2")
                return nothing
            end
        end
    end
    println("Files $infile1 and $infile2 are the same by this checking")
    return nothing
end

##================================================================##
function test1A()
##================================================================##
    filelist1 = ["ebasis_jl.inp","coords_jl.inp","gto_jl.inp","movec_jl.inp","roots_jl.inp","mo_to_ao_jl.inp","moeng_jl.inp","occupancy_jl.inp","spin_jl.inp"]
    filelist2 = ["ebasis_py.inp","coords_py.inp","gto_py.inp","movec_py.inp","roots_py.inp","mo_to_ao_py.inp","moeng_py.inp","occupancy_py.inp","spin_py.inp"]
    for ifile=1:length(filelist1)
        infile1 = filelist1[ifile]
        infile2 = filelist2[ifile]
        infile_checker(infile1,infile2)
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









