module Cstr
using LightXML, Printf
using DifferentialEquations, Sundials
include("Utils.jl")
include("Constants.jl")
include("SurfaceReactions.jl")

struct ConstParams{T1}
    ρ_in # inlet density
    q_in    #inlet flow rate 
    avg_molwt_in   #average molecular weight at the inlet
    As  #surface area
    V   #reactor volume
    mass_fracs_in::Array{T1,1}  #inlet mass fractions
    T   # temperature 
    p   # pressure
    ρq_yk::Array{T1,1} #inlet mass flow rate ρ x q x y_k
end

function cstr(input_file::AbstractString, lib_dir::AbstractString, sens= false)

    
    local gasphase = Array{String,1}

    xmldoc = parse_file(input_file)
    
    xmlroot = root(xmldoc)

    #get the gasphase from xml input
    gasphase = get_collection_from_xml(xmlroot,"gasphase")

    #create the thermo object
    thermo_file = lib_dir*"therm.dat"
    thermo_obj = SurfaceReactions.IdealGas.create_thermo(gasphase,thermo_file)

    #get the molefractions from xml input
    mole_fracs = get_molefraction_from_xml(xmlroot,thermo_obj.molwt,gasphase)
    mass_fracs_in = zeros(length(gasphase))
    molefrac_to_massfrac!(mass_fracs_in,mole_fracs,thermo_obj.molwt)
    avg_molwt = average_molwt(mole_fracs,thermo_obj.molwt)
    #get the operating temperature
    local T = get_value_from_xml(xmlroot,"T")

    #get the initial pressure 
    local p = get_value_from_xml(xmlroot,"p")

    #get the reactor volume 
    local V = get_value_from_xml(xmlroot,"volume")

    #get the flow rate. The flow rate is expected in m3/s
    local q = get_value_from_xml(xmlroot,"flow-rate")

    #get the surface area
    local As = get_value_from_xml(xmlroot,"As")

    #get the mechanism input file
    local mech_file = get_text_from_xml(xmlroot,"surface_mech")
    mech_file = lib_dir*"/"*mech_file
    md = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)

    #get the integration time
    it = get_value_from_xml(xmlroot,"time")
    
    #density at the inlet 
    ρ = density(mole_fracs,thermo_obj.molwt,T,p)
    #toal number of species
    n_species = length(gasphase) + length(md.sm.species)
    #definition of array for the storage of reaction rates
    rates = zeros(n_species)
    #solution vector
    sol = copy(mass_fracs_in)
    #initial coverage as specified in the mechanism
    covg = md.sm.si.ini_covg
    #append the coverage to the solution vector
    append!(sol,md.sm.si.ini_covg)


    #file output stream for saving the data
    g_stream = open("gas_profile.dat","w")
    s_stream = open("surf_profile.dat","w")
    o_streams = (g_stream, s_stream)
    create_header(g_stream,["t","T","p","rho"],gasphase)    
    create_header(s_stream,"t","T",md.sm.species)

    #create state object 
    state = SurfaceReactions.State(mole_fracs,covg,T,p,rates)
    #inlet conditions    
    t_span = (0,it)    
    #the following is the first term in the right hand side of gasphase species residual
    ρq_mass_fracs = mass_fracs_in * (ρ*q)
    cp = ConstParams(ρ,q,avg_molwt,As,V,mass_fracs_in, T, p,ρq_mass_fracs)
    params = (state, thermo_obj, md, cp, o_streams)
    prob = ODEProblem(residual!,sol,t_span,params)
    if sens == true
        return (params, prob, t_span)
    end
    cb = FunctionCallingCallback(save_data)
    soln = solve(prob,alg_hints=[:stiff] , reltol=1e-6, abstol = 1e-8, save_everystep=false, callback=cb)
    println("Solver Integration: ", soln.retcode, "\t")
    #close the files
    close(g_stream)
    close(s_stream)
    
end

function residual!(du,u,p,t)
    state, thermo_obj, md, cp, o_streams = p
    ng = length(state.mole_frac)
    ns = length(md.sm.species)
    mass_fracs = u[1:ng]
    #update the state with latest mole fractions
    massfrac_to_molefrac!(state.mole_frac,mass_fracs,thermo_obj.molwt)
    #update the state with latest coverages
    state.covg = u[ng+1:ng+ns]
    #calculate the density
    ρ = density(state.mole_frac,thermo_obj.molwt,cp.T,cp.p)
    #calculate the molar production rates
    SurfaceReactions.calculate_molar_production_rates!(state,thermo_obj,md)
    #outlet volumetric flow rate 
    q_out = cp.q_in*cp.avg_molwt_in/average_molwt(state.mole_frac,thermo_obj.molwt)
    #gasphase species residual
    spIn = cp.ρq_yk/cp.V
    spOut = q_out * ρ * mass_fracs /cp.V
    rgVec = (state.s_rate[1:ng] .* thermo_obj.molwt)*cp.As/cp.V    
    du[1:ng] = (spIn - spOut + rgVec)/ρ
    #surface species residual 
    du[ng+1:ng+ns] = (state.s_rate[ng+1:ng+ns] .* md.sm.si.site_coordination)/(md.sm.si.density*1e4)            
end
    

function save_data(u,t,integrator)
    state = integrator.p[1]
    thermo_obj = integrator.p[2]
    g_stream, s_stream = integrator.p[5]
    d = density(state.mole_frac,thermo_obj.molwt,state.T,state.p)
    write_to_file(g_stream,t,state.T,state.p,d,state.mole_frac)
    write_to_file(s_stream,t,state.T,state.covg)    
    @printf("%.4e\n",t)
end


end