@enum RxnType stick=1 arrhenius=2

"""
An anstract type for surface reaction mechanism parameter definitions
"""
abstract type SurfaceMechDefinition end
export SurfaceMechDefinition
abstract type ReactionState end
export ReactionState
"""
composite type for the definition of site definition of a surface reaction mechanism
"""
struct SiteInfo <: SurfaceMechDefinition
    name::String    # name of the surface reaction mechanism
    density::Float64    # surface site density in mol/cm2    
    site_coordination::Array{Float64,1} # site coordination of adsorbed species by default 1 for all
    ini_covg::Array{Float64,1} # intitial surface coverage as read from the mechanism
end
export  SiteInfo
mutable struct Parameters <: SurfaceMechDefinition
    s::Float64  # stickig coefficient
    k0::Float64 # pre-exponential factor cm-mol-s
    β::Float64  # temperature exponent
    E::Float64  # Activation energy in J/mol
end
export  Parameters
struct Stoichiometry    <: SurfaceMechDefinition
    reversible::Bool    #whether the reaction is reversible or not
    reactant_ids::Array{Int64,1} # id of reactant species participating in an elementary reaction
    product_ids::Array{Int64,1} #id of product species prarticipating in an elementary reaction
    cov_dep::Dict{Int64,Float64} #coverage dependency of the reaction: species_id=>cov dependent activation energy   
    order_dep::Dict{Int64,Float64} #order dependency of the reaction: species_id=>order w.r.t the species
end
export Stoichiometry
struct ElementaryReactions <: SurfaceMechDefinition
    id::Int64 #reaction number
    type::RxnType #reaction type stick or arrhenius
    params::Parameters #composite type Parameters
    mwc::Bool #whether MWC is applied or not      
    rxn_stoic::Stoichiometry #composite type Stoichiometry
end
export ElementaryReactions

struct SurfaceMechanism <: SurfaceMechDefinition
    energy_factor::Float64 # Conversion factor for activation energy to SI units
    si::SiteInfo # Composite type SiteInfo
    species::Array{String,1} #list of surface species present in the mechanism
    reactions::Array{ElementaryReactions,1} #Array of composite type ElementaryReactions
end
export  SurfaceMechanism

struct SpeciesRxnMap <: SurfaceMechDefinition
    species_id::Int64 #id of the species
    rxns::Array{Int64,1}    #id of reactions in which the species participates
    stoic_coeff::Array{Int64,1} #stoichiometric coefficient of the species in the reaction
end
export SpeciesRxnMap

struct MechanismDefinition <: SurfaceMechDefinition
    sm::SurfaceMechanism
    srm_array::Array{SpeciesRxnMap,1}
end

mutable struct State <: ReactionState
    mole_frac::Array{Float64,1}
    covg::Array{Float64,1}
    T::Float64
    p::Float64 
    s_rate::Array{Float64,1}
end 
export State

"""
Separate the reactants and products from the stoichiometric equation
#   Usage:
    parase_rxn_string(rxn_str)
-   rxn_str::String:    String representing the reaction equation
"""
function parase_rxn_string(rxn_str::String)
    #check whether the reaction is reversible
    reversible = false
    m = match(r"<=>|=>|=|>",rxn_str)
    if m.match == "=" || m.match == "<=>"
        reversible = true
    end
    species = collect(split(rxn_str,m.match))    
    r_species = [String(strip(uppercase(i))) for i in split(species[1],"+")]
    p_species = [String(strip(uppercase(i))) for i in split(species[2],"+")]    
    return (reversible,r_species,p_species)
end


function parse_rxn_params(params::String)
    k = β = E = 0.0
    all_params = split(params)
    k = parse(Float64,all_params[1])
    if length(all_params) == 3
        β = parse(Float64,all_params[2])
        E = parse(Float64,all_params[3])
    elseif length(all_params)==2
        E = parse(Float64,all_params[3])
    end
    return k, β, E
end

"""
Function to determine the coverage dependency, order dependency and MWC
for a reaction
#   Usage
    get_dependencies(cov_dep_rxns,order_dep_rxns,mwc_rxns,rxn_id)
-   cov_dep_rxns::Dict{Int64=>Dict{Int64=>Float64}}
-   order_dep_rxns::Dict{Int64=>Dict{Int64=>Float64}}
-   mwc_rxns::Array{Int64}
-   rxn_id::Int64
The function returns a tuple of cov,order and mwc
"""
function get_dependencies(cov_dep_rxns,order_dep_rxns,mwc_rxns,rxn_id)
    #Check for coverage dependency
    cov = Dict{Int64,Float64}()
    ord = Dict{Int64,Float64}()
    if haskey(cov_dep_rxns,rxn_id)
        cov = cov_dep_rxns[rxn_id]
    end
    #Check for order dependency
    if haskey(order_dep_rxns,rxn_id)
        ord = order_dep_rxns[rxn_id]
    end
    #Check for mwc
    mwc = in(rxn_id,mwc_rxns)
    return cov,ord,mwc

end