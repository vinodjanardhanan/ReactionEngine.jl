module Sens 
using LightXML, Printf
using DifferentialEquations, Statistics, GlobalSensitivity
using QuasiMonteCarlo

using ..Reactions
using ..SurfaceReactions
using ..Utils
using ..Cstr
using ..Batch
using ..Plug


struct  gsa_parameters 
    monitor::Array{String,1}   
    lower_pc::Float64
    upper_pc::Float64
    N::Int64
    gsa_srxn_ids::Array{Int64,1}
    gsa_srxn_constraint_ids::Array{Int64,1}    
end
struct Mechanisms 
    smd::MechanismDefinition
end

function gsa_xml(input_file::AbstractString, lib_dir::AbstractString)
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    # list of gasphase species 
    gsa_species = get_collection_from_xml(xmlroot,"gsa_species")
    # lower bound for parameters (in percentages)
    lower_pc = get_value_from_xml(xmlroot,"gsa_lb")
    # upper bound for parameters (in percentages)
    upper_pc = get_value_from_xml(xmlroot,"gsa_ub")
    
    # Size of the sample
    sample_N = Int(get_value_from_xml(xmlroot,"gsa_N"))

    # get the reactions id which forms the parameter 
    gsa_srxn_ids,gsa_srxn_constraint_ids = get_rxn_ids(xmlroot,"p_smech")  
        
    
    gsa_p = gsa_parameters(gsa_species,lower_pc,upper_pc,sample_N,gsa_srxn_ids,gsa_srxn_constraint_ids)
    

    # get the reactor model 
    model = get_text_from_xml(xmlroot,"gsa_model")
    
    model_root = root(parse_file(model))
    gsa_model = name(model_root)

    

    if lowercase(strip(gsa_model)) == "batch"
        params, prob, t_span = Batch.batch(model,lib_dir,true)        
        gsa_problem(params, prob, t_span, gsa_p)       
    elseif lowercase(strip(gsa_model)) == "plug"
        params, prob, t_span = Plug.plug(model,lib_dir, true)
        gsa_problem(params, prob, t_span, gsa_p)
    elseif lowercase(strip(gsa_model)) == "cstr"
        params, prob, t_span = Cstr.cstr(model,lib_dir,true)
        gsa_problem(params, prob, t_span, gsa_p)
    else
        println("Models yet to be implemented\n")
    end
end
    

"""
A common function for all reactor models to define the sensitivity problem 
#   Usage 
    gsa_problem(params, prob, gsa_p)
-   params : a tuple of parameters returned by the respective reactor model invoked 
-   prob : the problem definition for the reactor model 
-   gsa_p : struct of the type gsa_parameters
"""
function gsa_problem(params, prob, t_span, gsa_p::gsa_parameters)
    gsa_smech_params = Array{Float64,1}()
    state, thermo_obj, md, cp, o_stream = params        
    n_surface_species = length(md.sm.species) # number of surface species 
    n_gaspahse_species = length(thermo_obj.thermo_all) # number of gasphase species 
    n_species = n_gaspahse_species + n_surface_species # total number of species
    # get the parameters from the surface reaction mechaism  
    gsa_smech_parameter_map!(gsa_smech_params, gsa_p.gsa_srxn_ids, gsa_p.gsa_srxn_constraint_ids, md)
    println("Parameter size: ", length(gsa_smech_params))
    # create the lower limit for the parameters 
    lb = gsa_smech_params .* (1-0.01*gsa_p.lower_pc)
    # create the upper limit for the parameters 
    ub = gsa_smech_params .* (1+0.01*gsa_p.upper_pc)    
    parameter_ratio = Array{Float64,1}()  # Array for storing the original parameter ratio (forward to reverse)
    if count(x->x>0, gsa_p.gsa_srxn_constraint_ids) > 0
        # get the initial parameter ratios if the reverse reaction constraints are specified
        initial_parameter_ratio!(parameter_ratio, gsa_p.gsa_srxn_ids, gsa_p.gsa_srxn_constraint_ids,md)        
    end    
    
    perform_gsa(params,prob, t_span, gsa_p, lb, ub, parameter_ratio)
end


function perform_gsa(params, prob, t_span, gsa_p::gsa_parameters,lb::Array{Float64,1}, ub::Array{Float64,1}, parameter_ratio::Array{Float64,1})
    state, thermo_obj, md, cp, o_stream = params 
    time_points = collect(range(t_span[1],stop=t_span[2],length=50))
    #Get the name of all gasphase species 
    species_list = Array{String,1}()
    for sp_thermo in thermo_obj.thermo_all
        push!(species_list,sp_thermo.name)
    end
    append!(species_list,md.sm.species)
    #get the id of species to be monitored from the gsa_species list
    monitor_species_ids = map(x->get_index(x,species_list),gsa_p.monitor)    
    gsa_out = zeros(length(monitor_species_ids))
    func_eval_count = 0
    sens = function (gsa_params::Array{Float64})                
        func_eval_count += 1
        println("Sample @ ", func_eval_count)
        SurfaceReactions.update_params!(md,gsa_params,gsa_p.gsa_srxn_ids,gsa_p.gsa_srxn_constraint_ids,parameter_ratio)
        params_updated = (state, thermo_obj, md, cp, o_stream)
        s_prob = remake(prob,p=params_updated)                
        sol = solve(s_prob,alg_hints=[:stiff] , reltol=1e-6, abstol = 1e-8, saveat=time_points)  
        for i in eachindex(monitor_species_ids)
            gsa_out[i] = mean(sol[monitor_species_ids[i],:])
        end
        gsa_out
    end
    sampler = SobolSample()    
    A,B = QuasiMonteCarlo.generate_design_matrices(gsa_p.N,lb,ub,sampler)
    sobol_indices = gsa(sens,Sobol(order=[0,1,2]),A,B)

    # output file
    sobol_total = open("sobol_total.dat","w")     
    sobol_first = open("sobol_first.dat","w")
    sobol_second = open("sobol_second.dat","w")
    header_data = [ "p"*string(k) for k in gsa_p.gsa_srxn_ids]
    create_header(sobol_total,["species"],header_data)
    create_header(sobol_first,["species"],header_data)
    create_header(sobol_second,["species"],header_data)

    for i in eachindex(monitor_species_ids)
        write_to_file(sobol_total,gsa_p.monitor[i],sobol_indices.ST[i,:])
        write_to_file(sobol_first,gsa_p.monitor[i],sobol_indices.ST[i,:])
        write_to_file(sobol_second,gsa_p.monitor[i],sobol_indices.ST[i,:])
    end    

    close(sobol_total)
    close(sobol_first)
    close(sobol_second)
    
    
end

"""
A function to extract the parameters for sensitivity analysis from the xml input 
    The function returnds two integer arrays. forward_rxn_ids contains the ids of the 
    reactions whoes parameters will be considered as sensitivity parameters.
    reverse_rxn_ids contains the ids of the reactions aginst with the forward reaction 
    parameters are correlated. 
"""
function get_rxn_ids(xmlroot::XMLElement, xmltag::AbstractString)
    N = 10000
    p_mech = get_elements_by_tagname(xmlroot,xmltag)
    p_mech_content = ""
    if is_elementnode(p_mech[1])
        p_mech_content = content(p_mech[1])        
    end    
    p_mech_content_split = split(p_mech_content,",")           
    forward_rxn_ids = Array{Int64,1}()
    reverse_rxn_ids = Array{Int64,1}()
    if strip(p_mech_content_split[1]) == ":"        
        return [i for i in 1:N], [0 for i in 1:N]
    end
    for items in p_mech_content_split
        if occursin("=",items)
            rxn_sets = split(items,"=")        
            forward_sets = rxn_sets[1]
            reverse_sets = rxn_sets[2]
            forward_sets = map(x->parse(Int64,x),split(forward_sets,":"))
            reverse_sets = map(x->parse(Int64,x),split(reverse_sets,":"))
            if length(forward_sets) > 1
                for ids in forward_sets[1]:forward_sets[2]
                    push!(forward_rxn_ids,ids)
                end
                for ids in reverse_sets[1]:reverse_sets[2]
                    push!(reverse_rxn_ids,ids)
                end
            else                
                push!(forward_rxn_ids,forward_sets[1])
                push!(reverse_rxn_ids,reverse_sets[1])
            end       
        else            
            forward_sets = map(x->parse(Int64,x),split(items,":"))            
            if length(forward_sets) > 1
                for ids in forward_sets[1]:forward_sets[2]
                    push!(forward_rxn_ids,ids)
                    push!(reverse_rxn_ids,0)
                end
            else
                push!(forward_rxn_ids,forward_sets[1])
                push!(reverse_rxn_ids,0)
            end
        end
    end
    return (forward_rxn_ids, reverse_rxn_ids)
end

"""
This function returns the sticking coefficient or the pre-exponent factor 
    corresponding to the reaction ids present in gsa_rxn_ids 
#   Usage
    gsa_smech_parameter_map!(gsa_smech_params, gsa_srxn_ids, gsa_srxn_constraint_ids, md )    
-   gsa_smech_params : vector containing the parameters extracted from the inoput mechanism (output)
-   gsa_srxn_ids : reaction ids whoes parameters needs to be extracted from the mechanism (input)
-   gsa_srxn_constraint_ids : id of the reactions which needs to be considered for reverse reaction 
-   md: MechanismDefinition
"""
function gsa_smech_parameter_map!(gsa_smech_params::Array{Float64,1}, gsa_srxn_ids::Array{Int64,1}, gsa_srxn_constraint_ids::Array{Int64,1}, md )    
    
    function populate(rxn)
        if rxn.type == Reactions.stick
            push!(gsa_smech_params,rxn.params.s)
        else
            push!(gsa_smech_params,rxn.params.k0)
        end
    end
    # if the length of gsa_srxn_ids is very high then it means all reactions parameters need to be considered
    if length(gsa_srxn_ids) >= length(md.sm.reactions)
        for rxn in md.sm.reactions
            populate(rxn)
        end
        gsa_srxn_ids_copy = copy(gsa_srxn_ids)                
        empty!(gsa_srxn_ids)
        gsa_srxn_constraint_ids_copy = copy(gsa_srxn_constraint_ids)        
        empty!(gsa_srxn_constraint_ids)
        for i in 1:length(md.sm.reactions)
            push!(gsa_srxn_ids,gsa_srxn_ids_copy[i])
            push!(gsa_srxn_constraint_ids,gsa_srxn_constraint_ids_copy[i])
        end
        empty!(gsa_srxn_ids_copy)
        empty!(gsa_srxn_constraint_ids_copy)
    else
        for rxn in md.sm.reactions
            if rxn.id in gsa_srxn_ids
                populate(rxn)
            end        
        end
    end

    return gsa_smech_params
end


"""
A function to calculate the intial ratio of pre-exponential factors
"""
function initial_parameter_ratio!(p_ratio::Array{Float64,1},rxn_ids::Array{Int64,1},constraints_ids::Array{Int64,1}, md)
    for i in eachindex(rxn_ids)        
        if constraints_ids[i] != 0
            if md.sm.reactions[rxn_ids[i]].type == Reactions.stick
                ratio = md.sm.reactions[rxn_ids[i]].params.k0/md.sm.reactions[constraints_ids[i]].params.k0
                push!(p_ratio,ratio)
            else                
                ratio = md.sm.reactions[rxn_ids[i]].params.k0/md.sm.reactions[constraints_ids[i]].params.k0
                push!(p_ratio,ratio)
            end
        else
            push!(p_ratio,0.0)
        end
    end
end

#end of module
end