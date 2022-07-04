module Plug
using LightXML, Printf
using DifferentialEquations, Sundials
using ..SurfaceReactions
using ..Reactions
using ..Utils

"""
Plug flow
"""
function plug(input_file::AbstractString, lib_dir::AbstractString, sens= false)

    #locate the lib directory
    
    thermo_file = lib_dir*"therm.dat"

    local gasphase = Array{String,1}

    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    #Get gasphase from the xml input
    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    thermo_obj = SurfaceReactions.IdealGas.create_thermo(gasphase,thermo_file)        

    #Get the molefractions
    mole_fracs = get_molefraction_from_xml(xmlroot, thermo_obj.molwt, gasphase)
    
    
    #Convert mole fractions to mass fractions
    mass_fracs = zeros( length(gasphase))
    molefrac_to_massfrac!(mass_fracs,mole_fracs,thermo_obj.molwt)
    

    #Read the inlet temperature
    local T = get_value_from_xml(xmlroot,"T")

    #Read the velocity 
    local u = get_value_from_xml(xmlroot,"u")

    #Read the pressure
    local p = get_value_from_xml(xmlroot,"p")

    #Read the length of the reactor
    local l = get_value_from_xml(xmlroot,"length")

    #Read the diameter of the reactor 
    local dia = get_value_from_xml(xmlroot,"dia")
    local Ac = π*dia^2/4.0
    local As_per_unit_length = π*dia

    #Read the wall temperature
    local Tw = get_value_from_xml(xmlroot,"Tw")

    #Read the surface area per unit length
    local cat_geom = get_value_from_xml(xmlroot,"cat-geom-factor")

    #Read the temperature condition     
    isothermal = lowercase(get_text_from_xml(xmlroot,"isothermal")) == "true" ? true : false
    

    #Get the mechanism file from xml
    local mech_file = get_text_from_xml(xmlroot,"surface_mech")
    mech_file = lib_dir*"/"*mech_file

    #Create the mechanism definition
    md = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
    

    #create the soln vector
    soln = copy(mass_fracs) #mass fractions
    push!(soln,density(mole_fracs,thermo_obj.molwt,T,p)*u) #mass flux ρu
    if !isothermal
        push!(soln,T)
    end

    #file output streams for saving the data
    g_stream = open("gas_profile.dat","w")
    s_stream = open("surf_profile.dat","w")
    os_streams = (g_stream,s_stream)
    geometry = (Ac,As_per_unit_length,cat_geom)    

    create_header(g_stream,["z","T","p","u","rho"],gasphase)
    create_header(s_stream,"z","T",md.sm.species)
    #Set up the problem
    t_span = (0,l)
    covg = md.sm.si.ini_covg
    n_species = length(gasphase)+length(md.sm.species)
    rate = zeros(n_species)      
    #Create the state object    
    state = SurfaceReactions.SurfaceRxnState(mole_fracs,covg,T,p,rate)
    #this step is to get the steady state coverage before starting the plug flow integration
    t, state = SurfaceReactions.calculate_ss_molar_production_rates!(state,thermo_obj,md,1.0)
    params = (state,thermo_obj,md, geometry, os_streams)

    prob = ODEProblem(residual!,soln,t_span,params)
    if sens == true
        return (params, prob, t_span)
    end
    cb = FunctionCallingCallback(write_out_put)    
    sol = solve(prob, CVODE_BDF(linear_solver=:GMRES), reltol=1e-6, abstol=1e-10, save_everystep=false,callback=cb);        
    
    close(s_stream)
    close(g_stream)

    return sol.retcode
end

"""
residual function that defined the governing equations    
"""
function residual!(du,u,p,t)    
    #unpack the parameters
    state,thermo_obj,md,geom,os_stream = p
    Ac,Aspul,cat_geom = geom
    
    #Convert massfractions to molefractions
    ng = length(state.mole_frac)    
    massfracs = u[1:ng]

    massfrac_to_molefrac!(state.mole_frac,massfracs,thermo_obj.molwt)
    
    
    #Calculate the molar production rates
    SurfaceReactions.calculate_ss_molar_production_rates!(state,thermo_obj,md,1.0)
    
    #species residual    
    du[1:ng] = (state.s_rate[1:ng] .* thermo_obj.molwt) * cat_geom*(Aspul/Ac)/u[ng+1]
  
    #mass flux
    du[ng+1] = sum(state.s_rate[1:ng] .* thermo_obj.molwt)*cat_geom*Aspul/Ac
end



function write_out_put(u,t,integrator)    
    state = integrator.p[1]
    thermo_obj = integrator.p[2]
    g_stream, s_stream = integrator.p[5]
    
    d = density(state.mole_frac,thermo_obj.molwt,state.T,state.p)
    vel = u[length(state.mole_frac)+1]/d
    write_to_file(g_stream,t,state.T,state.p,vel,d,state.mole_frac)
    write_to_file(s_stream,t,state.T,state.covg)
    @printf("%.4e\n", t)   
end

end