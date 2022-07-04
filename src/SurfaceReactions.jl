module SurfaceReactions
using LightXML
using DifferentialEquations
using ..IdealGas
using ..Reactions
using ..Utils
include("Constants.jl")

export compile_mech
export calculate_molar_production_rates!, calculate_ss_molar_production_rates!, covg_integration!
export update_params!, update_params

"""
Function to create composite type of SurfaceMechanism and
SpeciesRxnMap. In general the application programs shall only
call the compile_mech function and not directly the create_mech
and species_rxn_map function
#   Usage
    compile_mech(file_path,gas_species)
-   file_path::String : path to the mechanism xml file
-   gas_species::Array{String,1} : list of gasphase species
"""
function compile_mech(file_path::T, thermo_obj::IdealGas.SpeciesThermoObj, gas_species::Array{T,1}) where T <: AbstractString
    sm = create_mech(file_path,gas_species,thermo_obj.molwt)
    rxn_map = species_rxn_map(gas_species,sm)
    return SurfaceMechDefinition(sm,rxn_map)
end

"""
Function to read a surface reaction mechanism file
#   Usage:
    create_mech(file_path,gas_species)
-   file_path::String : file name including the relative path 
-   gas_species::Array{String,1} : gasphase species list 
"""
function create_mech(file_path::String, gas_species::Array{String,1},molwt::Array{Float64,1})
    xmldoc = parse_file(file_path)
    mech_root = root(xmldoc)    
    mech_unit = attribute(mech_root,"unit")    
    energy_factor = convert2si(mech_unit)

    #get the surface species
    surface_species = map(x->uppercase(x) ,split(content(get_elements_by_tagname(mech_root,"species")[1])))

    #default initial values
    site_coordination_vector = ones(length(surface_species))
    covg_vec = zeros(length(surface_species))
    site_density = 1.0
    site_name = ""
    #read the site node
    site_node = get_elements_by_tagname(mech_root,"site")
    for elements in site_node
        site_name = attribute(elements,"name")
        #get the site coordination
        coordination = get_elements_by_tagname(elements,"coordination")[1]
        if is_elementnode(coordination)
            site_coord_species = split(content(coordination),",")
            for items in site_coord_species
                sp,val = split(items,"=")
                sp_id = get_index(String(sp),surface_species)
                site_coordination_vector[sp_id] = parse(Float64,String(val))
            end
        end        
        #get the site density
        density = get_elements_by_tagname(elements,"density")[1]
        if is_elementnode(density)
            unit = strip(attribute(density,"unit"))
            site_density = parse(Float64, content(density))
            if uppercase(unit) == "MOL/M2"
                site_density *= 1e4
            end
        else
            throw(error("Surface site density not specified\n"))
        end
    
        #get the intitial coverages
        covg_node = get_elements_by_tagname(elements,"initial")[1]
        if is_elementnode(covg_node)
            ini_covg = split(content(covg_node),",")
            for items in ini_covg
                sp,val = split(items,"=")
                sp_id = get_index(String(sp),surface_species)
                covg_vec[sp_id] = parse(Float64,String(val))
            end
        else
            throw(error("Initial coverage not specified\n"))
        end
        
    end

    site_info = SiteInfo(site_name,site_density,site_coordination_vector,covg_vec)

    #Read the coverage dependencies if any. This must be read before reading the reactions
    covg_nodes = get_elements_by_tagname(mech_root,"coverage")
    cov_dep_rxns = Dict{Int64,Dict{Int,Float64}}() #rxn_id=>{species_id=>cov_act_energy}
    for covg in covg_nodes
        rxn_ids = split(attribute(covg,"id"))
        sp,val = split(content(covg),"=")
        sp_id = get_index(String(sp),surface_species)
        cov_dict = Dict(sp_id=>parse(Float64,val))

        for rxn in rxn_ids
           cov_dep_rxns[parse(Int64,rxn)] = cov_dict
        end        
    end
    #Read the order dependencies if any. This must be read before reading the reactions
    order_nodes = get_elements_by_tagname(mech_root,"order")
    order_dep_rxns = Dict{Int64,Dict{Int,Float64}}() #rxn_id=>{species_id=>order}
    for order in order_nodes
        rxn_ids = split(attribute(order,"id"))
        sp,val = split(content(order),"=")
        sp_id = get_index(String(sp),surface_species)
        order_dict = Dict(sp_id=>parse(Float64,val))
        for rxn in rxn_ids
            order_dep_rxns[parse(Int64,rxn)] = order_dict
        end
    end
    

    #Read Motz-Wise Correction
    mwc_node = get_elements_by_tagname(mech_root,"mwc")
    mwc_rxns = Array{Int64,1}()
    if length(mwc_node) > 0
        if is_elementnode(mwc_node[1])
            mwc_rxns = [ parse(Int64,rxn) for rxn in split(content(mwc_node[1]))]
        end    
    end
    
    #Define the reaction array to be returned
    rxn_array= Array{SurfaceRxns,1}()

    rct_ids = Array{Int64,1}()
    prdt_ids = Array{Int64,1}()
    all_species = Array{String,1}()
    all_species = copy(gas_species)
    append!(all_species,surface_species)
    
    #Read sticking reactions
    site_index = get_index(site_name,all_species)

    stick_node = get_elements_by_tagname(mech_root,"stick")[1]
    if is_elementnode(stick_node)
        stick_rxns = get_elements_by_tagname(stick_node,"rxn")
        for rxn in stick_rxns
            rxn_id = parse(Int64,strip(attribute(rxn,"id")))
            stoic, params = split(content(rxn),"@")
            reversible, r_species, p_species = parase_rxn_string(String(stoic))  
            rct_ids = [get_index(sp,all_species) for sp in r_species]          
            prdt_ids = [get_index(sp,all_species) for sp in p_species]
            s,β,E = parse_rxn_params(String(params))            
            #Check for coverage, order and mwc dependency
            cov,ord,mwc = get_dependencies(cov_dep_rxns,order_dep_rxns,mwc_rxns,rxn_id)
            #convert sticking coefficient to rate constant
            n_site = count(x->x==site_index, rct_ids)
            gas_species_id = filter(x->x != site_index, rct_ids)[1]
            sqrt_term = 100*sqrt(R/(2*π*molwt[gas_species_id]))
            s_by_gamma = s/site_density^n_site            
            k = sqrt_term*s_by_gamma #rate constant
            #Apply the MWC
            if mwc
                k /= (1-0.5*s)
            end
            p = SurfaceRxnParameters(s,k,β,E)

            stoic = SurfaceRxnStoichiometry(reversible,rct_ids,prdt_ids,cov,ord)
            elm_rxn = SurfaceRxns(rxn_id,Reactions.stick,p,mwc,stoic)
            push!(rxn_array,elm_rxn)
            
        end
    end
    
    #Read Arrhenius reactions
    arrhenius_node = get_elements_by_tagname(mech_root,"arrhenius")[1]
    if is_elementnode(arrhenius_node)
        arrhenius_rxns = get_elements_by_tagname(arrhenius_node,"rxn")
        for rxn in arrhenius_rxns
            rxn_id = parse(Int64,strip(attribute(rxn,"id")))
            stoic,params = split(content(rxn),"@")
            reversible, r_species, p_species = parase_rxn_string(String(stoic))
            rct_ids = [get_index(sp,all_species) for sp in r_species]
            prdt_ids = [get_index(sp,all_species) for sp in p_species]
            k,β,E = parse_rxn_params(String(params))
            p = SurfaceRxnParameters(0.0,k,β,E)
            #Check for coverage, order and mwc dependency. MWC not applicable to Arrhenius type
            cov,ord,mwc = get_dependencies(cov_dep_rxns,order_dep_rxns,mwc_rxns,rxn_id)   
            stoic = SurfaceRxnStoichiometry(reversible,rct_ids,prdt_ids,cov,ord)
            elm_rxn = SurfaceRxns(rxn_id,Reactions.arrhenius,p,mwc,stoic)  
            push!(rxn_array,elm_rxn)       
        end
    end
    
    return SurfaceMechanism(energy_factor,site_info,surface_species,rxn_array)
end


"""
Create array of structures of SpeciesRxnMap. The SpeciesRxnMap is a 
composit structure of species_id in participating reactions and the
stoichiometric coefficients
#   Usage
    species_rxn_map(gas_species,sm)
-   gas_species::Array{String,1} :  list of gas phase species
-   sm::SurfaceMechanism : SurfaceMechanism composite type
"""
function species_rxn_map(gas_species::Array{String,1},sm::SurfaceMechanism)
    #Array of SpeciesRxnMap
    srm_array = Array{SpeciesRxnMap,1}()
    #create list of all species
    all_species = Array{String,1}()
    all_species = copy(gas_species)
    append!(all_species,sm.species)
    #loop over all species
    for sp in all_species
        sp_id = get_index(sp,all_species)
        #check if the species is present in the reactant_ids of a reaction
        rxns = Array{Int64,1}()
        stc = Array{Int64,1}()
        for rxn in sm.reactions
            #count the number of occurances of species in reactants list (using ids)
            n = count(x->x==sp_id,rxn.rxn_stoic.reactant_ids)
            if n > 0
                push!(rxns,rxn.id)
                push!(stc,-n)
            end            
            #count the number of occurances of species in the products list (using ids)
            n = count(x->x==sp_id,rxn.rxn_stoic.product_ids)
            if n > 0
                push!(rxns,rxn.id)
                push!(stc,n)
            end
        end
        push!(srm_array,SpeciesRxnMap(sp_id,rxns,stc))
    end
    return srm_array
end



"""
calculate_molar_production_rates!(state::ReactionState, thermo_obj::IdealGas.SpeciesThermoObj,md::SurfaceMechDefinition)

Function to calculate the molar production rate of gas-phase and 
surface species based on the surface reaction mechanism. The function
will alter the field s_rate after the calculation.
#   Usage:
    calculate_molar_production_rate!(state,md)
-   state::State : state of the reaction, which is of the composite type State
-   thermo_obj::SpeciesThermoObj : Composite type SpeciesThermoObj. Required in case of reversible reactions
-   md::MechanismDefinition : the mechanism, which is of the composite type MechanismDefinition
"""
function calculate_molar_production_rates!(state::ReactionState, thermo_obj::IdealGas.SpeciesThermoObj,md::SurfaceMechDefinition)
    #get the mole fractions from state
    mf = state.mole_frac
    #convert mole fractions to concentrations
    all_conc = (state.p/R/state.T) .* mf    
    all_conc *= 1e-6 # convert to mol/cm3
    #get the surface coverages from state
    coverage = state.covg
    #convert coverages to surface concentrations mol/cm2    
    surf_conc =  (md.sm.si.density .* coverage) ./ md.sm.si.site_coordination
    #merge the concentrations
    append!(all_conc,surf_conc)    
    rxn_rate = Array{Float64,1}()
    
    #calcuate the rate of individual reactions
    for rxn in md.sm.reactions        
        local prdt_conc = prod(all_conc[rxn.rxn_stoic.reactant_ids])        
        #check for coverage dependency of the reaction
        local cov_act_E = 0
        #cov_act_E = sum(coverage[collect(keys(rxn.rxn_stoic.cov_dep))] .* collect(values(rxn.rxn_stoic.cov_dep)))
        if !isempty(rxn.rxn_stoic.cov_dep)
            #calculate coverage dependent activation energy
            for (sp_id,cov_e) in rxn.rxn_stoic.cov_dep
                cov_act_E += coverage[sp_id]*cov_e
            end                       
        end
        #check for order dependency
        local m_order = 1.0
        if !isempty(rxn.rxn_stoic.order_dep)
            #calculate order change
            for (sp_id,order) in rxn.rxn_stoic.order_dep
                m_order *= coverage[sp_id]^order
            end
        end
        #calculate the forward reaction rate                
        E = (rxn.params.E + cov_act_E)*md.sm.energy_factor        
        k = rxn.params.k0 * state.T^rxn.params.β * exp(-E/R/state.T)
        if rxn.type == Reactions.stick
            k *= sqrt(state.T)
        end        
        push!(rxn_rate,k*prdt_conc*m_order*1e4)#e4 is to convert the units to mol/m2

        if rxn.rxn_stoic.reversible
            throw(error("Reverse reaction not implemented\n"))
        end
        
    end
    #calculate the species molar production rates
    for srm in md.srm_array                            
        state.s_rate[srm.species_id] = sum(rxn_rate[srm.rxns] .* srm.stoic_coeff)
    end
end

"""
calculate_ss_molar_production_rates!(state::ReactionState, thermo_obj::IdealGas.SpeciesThermoObj,md::SurfaceMechDefinition, time =1.0)

Function to calculate the molar production rate of gas-phase species.
This function integrates the rate equations so that \frac{d\theta}{dt}=0
#   Usage:
    calculate_molar_production_rate!(state,md)
-   state::State : state of the reaction, which is of the composite type State
-   thermo_obj::SpeciesThermoObj : Composite type SpeciesThermoObj. Required in case of reversible reactions
-   md::MechanismDefinition : the mechanism, which is of the composite type MechanismDefinition
"""
function calculate_ss_molar_production_rates!(state::ReactionState, thermo_obj::IdealGas.SpeciesThermoObj,md::SurfaceMechDefinition, time =1.0)
    t_span = (0,time) 
    coverage = state.covg
    p = (state,thermo_obj,md)
    prob = ODEProblem(covg_integration!,coverage,t_span,p)
    sol = solve(prob, alg_hints=[:stiff] , reltol=1e-6, abstol=1e-10, save_everystep=false,save_start=false)        
    state.covg = sol.u[lastindex(sol.u)]
    return sol.t, state
end


function covg_integration!(du, u, p, t)
    state, thermo_obj,md = p
    state.covg = u    
    ng = length(state.mole_frac)
    calculate_molar_production_rates!(state,thermo_obj,md)         
    @inbounds for i in eachindex(u)
        du[i] = (md.sm.si.site_coordination[i]*state.s_rate[ng+i])/(md.sm.si.density*1e4)
    end
end


"""
Function to update the rate parameters. This function is only used for
sensitivity analysis
#   Usage
    update_params(rate_params,md)
-   rate_params::Array  : One dimentional array of floating point variables
-   md:SurfaceMechDefinition : SurfaceMechDefinition struct
"""
function update_params(rate_params::Array{Float64,1},md::SurfaceMechDefinition)
    #assign the new parameters to SurfaceMechDefinition object 
    for i in eachindex(rate_params)
        if md.sm.reactions[i].type == Reactions.stick
            k = md.sm.reactions[i].params.k0 #This is the current value of converted rate constant including MWC
            if md.sm.reactions[i].mwc == true
                #take out the MWC factor
                k *= (1-0.5*md.sm.reactions[i].s)
            end
            #new k value for the purturbed sticking coefficient
            k = (k/md.sm.reactions[i].params.s)*rate_params[i]
            md.sm.reactions[i].params.s = rate_params[i]
            if md.sm.reactions[i].mwc ==  true
                #Apply the MWC back
                k /= (1-0.5*rate_params[i])
            end
            md.sm.reactions[i].params.k0 = k
        else
            md.sm.reactions[i].params.k0 = rate_params[i]
        end                
    end

end

"""
A function for use in global sensitivity analysis. The function updates the parameters with the 
    supplied rate_params 
#   Usage
    update_params!(md, rate_params,rxn_ids, constraint_ids)
-   md : SurfaceMechDefinition structure
-   rate_params : parameters to that needs to be set to the mechanism object 
-   rxn_ids : id of the reactions considered in the sensitivity analysis 
-   constraint_ids : id of the reverse reaction ids considered
"""
function update_params!(md, rate_params::Array{Float64,1},rxn_ids::Array{Int64,1}, constraint_ids::Array{Int64,1}, pratio::Array{Float64,1})
    for id in eachindex(rxn_ids)
        if md.sm.reactions[rxn_ids[id]].type == Reactions.stick
            # The current value of sticking coefficient converted to rate constant including MWC
            k = md.sm.reactions[rxn_ids[id]].params.k0                       
            if md.sm.reactions[rxn_ids[id]].mwc == true
                #take out the MWC factor
                k *= (1-0.5*md.sm.reactions[rxn_ids[id]].params.s)
            end
            # multiplication factor used to convert sticking coefficient to rate constant 
            m_factor = k/md.sm.reactions[rxn_ids[id]].params.s
            # new sticking coefficient 
            md.sm.reactions[rxn_ids[id]].params.s = rate_params[id]
            #new k value for the purturbed sticking coefficient
            k  = md.sm.reactions[rxn_ids[id]].params.s * m_factor
            #Apply the MWC back
            if md.sm.reactions[id].mwc == true                
                k /= (1-0.5*rate_params[id])
            end
            md.sm.reactions[rxn_ids[id]].params.k0 = k            
            
            # Now if there is a constraint, i.e. if a reverse reaction present, then the 
            # parameter for that reaction needs to be changed accordingly. 
            if constraint_ids[id] != 0                
                md.sm.reactions[constraint_ids[id]].params.k0 = k/pratio[id]            
            end
        else            
            md.sm.reactions[rxn_ids[id]].params.k0 = rate_params[id]            
            if constraint_ids[id] != 0                
                md.sm.reactions[constraint_ids[id]].params.k0 = rate_params[id]/pratio[id]                
            end
        end
    end      
end



#end of module
end