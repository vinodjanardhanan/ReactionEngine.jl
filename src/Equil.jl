module Equil
using LightXML, Printf
using LinearAlgebra, Roots, RowEchelon
include("IdealGas.jl")
include("Utils.jl")

function equil(input_file::AbstractString, lib_dir::AbstractString)
    #lib directory    
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)
    

    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    thermo_file = lib_dir*"therm.dat"
    thermo_obj = IdealGas.create_thermo(gasphase, thermo_file)
    all_species_thermo = thermo_obj.thermo_all
    #create elements from the species composition
    elements = Array{String,1}()
    for sp in all_species_thermo        
        append!(elements,collect(keys(sp.composition)))
    end
    unique!(elements)    
    mole_fracs  = get_molefraction_from_xml(xmlroot,thermo_obj.molwt,gasphase)
    local T = get_value_from_xml(xmlroot,"T")
    local p = get_value_from_xml(xmlroot,"p")
    
    #Total number of moles under the given conditions
    n_total = p/R/T
    #Number of moles of each component    
    moles = n_total*mole_fracs
    #Order the species according to inlet moles
    order_species(gasphase,moles)
    #Save the ordered list for later use
    gasphase_in = copy(gasphase)
    moles_in = copy(moles)
    
    #reconstruct the thermo obj based on the ordered list
    thermo_obj = IdealGas.create_thermo(gasphase, thermo_file)
    all_species_thermo = thermo_obj.thermo_all

    #Create the stoichiometric coefficients for the indepdendent reactions
    stc = stoichiometric_coeffient_matrix(gasphase,moles,elements,all_species_thermo)   
    stc = round.(stc,digits=2)
    #The above method removes the inerts from the list. Therefore store the inlet moles of inerts for final calculation 
    inert_species = setdiff(gasphase_in,gasphase)
    inert_moles = Array{Float64,1}()    
    for i in inert_species
        append!(inert_moles,moles_in[get_index(i,gasphase_in)])
    end
    #Create the thermo_obj for the new list of species (i.e without inerts)
    thermo_obj = IdealGas.create_thermo(gasphase, thermo_file)
    all_species_thermo = thermo_obj.thermo_all

    #Uncomment this for printing the independent reactions
    # print_reactions(stc,gasphase)

    #get the Gibb's free energy of all species to calculate Kp
    H_all = IdealGas.H_all(thermo_obj,T)
    S_all = IdealGas.S_all(thermo_obj,T)
    G_all = H_all - T*S_all        
    
    # #Calculate the Del G for each reaction 
    DelG = Array{Float64,1}()
    for i in 1:size(stc,1)        
        push!(DelG,sum(stc[i,:] .* G_all))
    end    
    molefracs = moles_in/sum(moles_in)
    println("\nInititial condition:\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase_in)
        @printf("%10s \t %.4e \t %.4e \n", gasphase_in[k], moles_in[k], molefracs[k])
    end
    Kp = exp.(-DelG/R/T)
    # Minimize the Gibb's energy of the mixture
    minimize_G!(moles, stc, thermo_obj, Kp, T, p)
    #Account for the inert moles 
    moles_all = copy(moles)
    append!(moles_all,inert_moles)
    molefracs = moles_all/sum(moles_all)
    gasphase_in = copy(gasphase)
    append!(gasphase_in,inert_species)
    println("\nEquilibrium composition @ T= $T K and p=$p Pa\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase_in)
        @printf("%10s \t %.4e \t %.4e\n", gasphase_in[k], moles_all[k],molefracs[k])        
    end

    return Symbol("Success")
end

"""
Bubble sort the species with the highest to lower order of moles
"""
function order_species(gasphase::Array{String}, moles::Array{Float64})
    for i in 1:length(gasphase)-1
        for j in i+1:length(gasphase)
            if moles[i] < moles[j]
                moles[i], moles[j] = moles[j], moles[i]
                gasphase[i], gasphase[j] = gasphase[j], gasphase[i]
            end
        end
    end
end

"""
This function creates the stoichiometric coefficients for 
    linearly indepdendent reactions based on null space
"""
function stoichiometric_coeffient_matrix(gasphase, moles, elements, all_species_thermo)
    #constuct the B matrix
    B = zeros(length(gasphase),length(elements))
    for i in eachindex(all_species_thermo)
        for (k,v) in all_species_thermo[i].composition
            j = get_index(String(k),elements)
            B[i,j] = v
        end
    end
    B = B'
    # Eliminate the inert element from B
    nr = Array{Int64,1}()
    for i in 1:size(B,1)
        if count(e->e>0, B[i,:]) < 2
            append!(nr,i)
        end
    end    
    ncount = 0    
    for i in nr
        B = B[1:end .!=i-ncount,:]
        ncount += 1
    end    
    # Eliminate the element constituting the inert from the elements list 
    deleteat!(elements,nr)    
    # Eliminate the empty columns (delete the inert)
    nc = Array{Int64,1}()        
    for i in 1:size(B,2)        
        if count(e->e>0,B[:,i]) < 1
            append!(nc,i)
        end
    end 
    ncount = 0
    for i in nc
        B = B[:,1:end .!= i-ncount]
        ncount += 1
    end    
    #The inert needs to be removed from the gasphase list and the moles as well
    deleteat!(gasphase,nc)
    deleteat!(moles,nc)
    
    if size(B,1) != size(B,2)
        n_null_space = length(gasphase)-rank(B)    
        A = zeros(n_null_space, length(gasphase)-n_null_space)
        AI = Matrix(I,n_null_space, n_null_space)
        B = vcat(B,hcat(A,AI))
        C = inv(B)
        nsv = C[:,end-n_null_space+1:end]        
        b = minimum(abs.(filter(e->e!=0,nsv)))
        nsv /= b
        return nsv'
    else
        B = rref(A)
        #check for rows containing only zero elements and augment with identity matrix 
        #to the right with size equal to the number of rows containing only zero as elements
        n_null_space = 0
        for i in 1:size(B,1)
            if count(e->e!=0, B(i,:)) == 0 
                B[i,i] = 1
                n_null_space += 1
            end
        end
        C = inv(B)
        nsv = C[:,end-n_null_space+1:end]        
        b = minimum(abs.(filter(e->e!=0,nsv)))
        nsv /= b
        return nsv'
    end

end

function print_reactions(stc::Matrix{Float64}, gasphase::Array{String})
    println("\nLinearly indepdendent reactions:\n")
    #Negative coefficients are assumed to be reactants     
    for i in 1:size(stc,1)
        rcts = ""
        prdts = ""
        for k in 1:length(gasphase)
            if stc[i,k] < 0 && stc[i,k] != 0
                rcts *= string(abs(stc[i,k]))*" "*gasphase[k]*" + "
            elseif stc[i,k] > 0 && stc[i,k] != 0
                prdts *= string(stc[i,k])*" "*gasphase[k]*" + "
            end
        end
        rcts = (strip(rcts))[1:end-1]
        prdts = (strip(prdts))[1:end-1]
        println(rcts * " > " * prdts)
    end
    
end

function minimize_G!(moles::Array{Float64},stc::Matrix, thermo_obj, Kp::Array{Float64} , T::Float64, p::Float64)        
    minimized = false
    nr = size(stc,1)
    Gmin = true 
    while Gmin
        mole_fracs = moles/sum(moles)                  
        G0 = IdealGas.Gmix(thermo_obj,T,p,mole_fracs)   
        for i in 1:nr
            equilibrate_reaction!(moles,stc[i,:],Kp[i],p)  
            #check if the Gibbs energy has decreased             
        end
        mole_fracs = moles/sum(moles)                  
        G = IdealGas.Gmix(thermo_obj,T,p,mole_fracs) 
        # println(G, "\t", G0, "\t", G-G0)     
        Gmin = G-G0 < 1e-15 ? false : true
    end
    # mole_fracs = moles/sum(moles)                  
    # println(mole_fracs)
    # readline()    
end

    
function equilibrate_reaction!(moles::Array{Float64}, rxn_stc::Array{Float64}, Kp::Float64, p::Float64 )
    # println("Equilibrating rxn: ", rxn )
    dn_min = -1e10
    dn_max = 1e10    
    #The maximum possible extent of a reaction depends on the limiting concentration that is present 
    for i in eachindex(rxn_stc)        
        if rxn_stc[i] < 0
            dn_max = min(dn_max,-moles[i]/rxn_stc[i])
        end
        if rxn_stc[i] > 0
            dn_min = max(dn_min,-moles[i]/rxn_stc[i])
        end
    end    
    
    local moles_updated = copy(moles)    

    function root_func(x)                
        n1 = (moles_updated .+ rxn_stc*x)                 
        mole_fracs = n1/sum(n1)
        n_pp = mole_fracs * (p/p_std)
        prdt_p = 1.0
        prdt_r = 1.0
        for k in eachindex(rxn_stc)
            if rxn_stc[k] < 0
                prdt_r *= n_pp[k]^abs(rxn_stc[k])
            elseif rxn_stc[k] > 0                
                prdt_p *= n_pp[k]^rxn_stc[k]
            end            
        end
        return (Kp*prdt_r - prdt_p)                
    end    
    dn = find_zero(root_func,(dn_max,dn_min),Bisection())
    #moles = moles_updated .+ rxn*dn
    moles .+= rxn_stc*dn
    
end

end