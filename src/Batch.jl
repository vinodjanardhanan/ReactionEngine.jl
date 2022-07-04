module Batch
using LightXML, Printf
using DifferentialEquations

using ..Utils
using ..SurfaceReactions

include("Constants.jl")

struct ConstantParams{T1}
    Asv::T1
    T::T1
end

"""
A common function is written for the reading input data
for all methods that utilizes batch reactor model. The input
    struct defines the common Parameters    
"""
struct InputData
    T::Float64
    p_initial::Float64
    Asv::Float64    
    tf::Float64
    gasphase::Array{String,1}
    mole_fracs::Array{Float64,1}
    thermo_obj::SurfaceReactions.IdealGas.SpeciesThermoObj
    md::SurfaceReactions.MechanismDefinition
end


"""
This is the calling function for executing the batch reactor 
#   Usage
    run_batch(input_file, lib_dir)    
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
"""
function batch(input_file::AbstractString, lib_dir::AbstractString, sens=false)

    
    
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)
    id = input_data(xmlroot, lib_dir)    
    #create solution vector 
    soln = get_solution_vector(id)
    #total number of species
    n_species = length(id.gasphase) + length(id.md.sm.species)
    #storage for reaction rates
    rates = zeros(n_species)
    #create output files for vasing the data
    g_stream = open("gas_profile.dat","w")
    s_stream = open("surface_covg.dat","w")  
    o_stream = (g_stream, s_stream)  
    create_header(g_stream,["t","T","p","rho"],id.gasphase)
    create_header(s_stream,"t", "T" ,id.md.sm.species)

    #define the Parameters
    cp = ConstantParams(id.Asv,id.T)    
    state = SurfaceReactions.SurfaceRxnState(id.mole_fracs,id.md.sm.si.ini_covg,id.T,id.p_initial,rates)
    t_span = (0,id.tf)
    params = (state, id.thermo_obj, id.md, cp, o_stream)
    prob = ODEProblem(residual!,soln,t_span,params)
    if sens == true
        return (params, prob,t_span)
    end
    cb = FunctionCallingCallback(save_data)
    sol = solve(prob,alg_hints=[:stiff] , reltol=1e-6, abstol = 1e-8, save_everystep=false, callback=cb)
    
    close(g_stream)
    close(s_stream)    

    return sol.retcode
end


"""
Function for creating the solition vector
"""
function get_solution_vector(id::InputData)
    soln = zeros(length(id.gasphase))
    molefrac_to_massfrac!(soln,id.mole_fracs,id.thermo_obj.molwt)
    soln .*= density(id.mole_fracs,id.thermo_obj.molwt,id.T,id.p_initial)
    append!(soln,id.md.sm.si.ini_covg)    
    return soln
end


"""
Function for reading the common input parameters
"""
function input_data(xmlroot::XMLElement, lib_dir::AbstractString)

    # locate the lib directory    
    thermo_file = lib_dir*"therm.dat"

    #get the gasphase present
    gasphase = Array{String,1}
    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    
    thermo_obj = SurfaceReactions.IdealGas.create_thermo(gasphase, thermo_file)

    #Get the mole fractions from xml
    mole_fracs = get_molefraction_from_xml(xmlroot, thermo_obj.molwt, gasphase)
    
    #Get the temperature
    local T = get_value_from_xml(xmlroot,"T")

    #get the initial pressure 
    local p_initial = get_value_from_xml(xmlroot,"p")

    #get the surface area to volume ratio 
    local Asv = get_value_from_xml(xmlroot,"Asv")

    #get the integration time 
    local tf = get_value_from_xml(xmlroot,"time")

    #get the mechanism file
    mech_file = get_text_from_xml(xmlroot,"surface_mech")
    mech_file = lib_dir*"/"*mech_file
    
    #create the mechanism definition
    md = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
    

    id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,md)
    
    return id

end


"""
Residual function definition for the batch reactor
"""
function residual!(du,u,p,t)        
    state, thermo_obj, md, cp  = p          

    ng = length(state.mole_frac)
    ns = length(md.sm.species)    

    #density 
    ρ = sum(u[1:ng])
    #mass fractions 
    mass_fracs = u[1:ng]/ρ
    #mole fractions
    massfrac_to_molefrac!(state.mole_frac,mass_fracs,thermo_obj.molwt)
    #average molecular weight
    mlwt_avg = average_molwt(state.mole_frac,thermo_obj.molwt)
    #pressure update
    state.p = ρ*R*cp.T/mlwt_avg
    #coverage update
    state.covg = u[ng+1:ng+ns]
    
    #calculate the molar production rates
    SurfaceReactions.calculate_molar_production_rates!(state,thermo_obj,md)
    #gasphase residual
    du[1:ng] = (state.s_rate[1:ng] .* thermo_obj.molwt) * cp.Asv
    #coverage residuals
    du[ng+1:ng+ns] = (state.s_rate[ng+1:ng+ns] .* md.sm.si.site_coordination)/(md.sm.si.density*1e4)
    
end



"""
Function to save the variables into output file
"""
function save_data(u,t,integrator)
    state = integrator.p[1]
    g_stream, s_stream = integrator.p[5]
    #density 
    ρ = sum(u[1:length(state.mole_frac)])
    write_to_file(g_stream,t,state.T,state.p,ρ,state.mole_frac)
    write_to_file(s_stream,t,state.T,state.covg)   
    @printf("%4e\n",t) 
end


end