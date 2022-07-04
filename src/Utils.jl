
module Utils

include("Constants.jl")
using LightXML, Printf

export convert2si, get_index, parse_composition
export massfrac_to_molefrac!, massfrac_to_molefrac, molefrac_to_massfrac!, average_molwt, density
export get_collection_from_xml, get_value_from_xml, get_text_from_xml, get_molefraction_from_xml
export create_header, write_to_file

cf = Dict("M"=>1,"CM"=>0.01,"KJ/MOL"=>1000)
convert2si(u::String) = cf[strip(uppercase(u))]


"""
Get the index of species in the gasphase limits
#   Usage
    get_index(sp,ig)
-   'sp::String' : species name
-   'ig::Gasphase' : gasphase object
"""
function get_index(sp::String,splist::Array{String,1})
    if uppercase(strip(sp)) in splist
        return findfirst(x-> x == uppercase(strip(sp)),splist)
    else
        throw(ErrorException("species $sp not present in $splist"))
    end
end


"""
Function to parse the composition string
#   Usage:
    parse_composition(composition_data,splist)
-   composition_data::String : composition specification String
-   sp_list:Array{String,1}: Array of species names
"""
function parse_composition(composition_data::String, sp_list::Array{String,1})
    species_composition = split(composition_data,",")
    composition = zeros(length(sp_list))
    for sc in species_composition
        species,val = split(sc,"=")
        composition[get_index(String(species), sp_list)] = parse(Float64,val)
    end
    return composition
end

"""
Function to convert massfractions to mole fractions
#   Usage:
    massfrac_to_molefrac!(massfracs,mol_wt,molefracs)
-   massfracs::Array{Float64,1} : Vector of massfractions
-   mol_wt::Array{Float64,1} : Vector of molecular weights
-   molefracs::Array{Float64,1} : Vector of molefractions 
"""
function massfrac_to_molefrac!(molefracs::T, massfracs::T, mol_wt::T) where T <: Array{Float64,1}
    moles = massfracs ./ mol_wt        
    molefracs .= moles/sum(moles)
end


"""
Function to convert massfractions to mole fractions
#   Usage:
    massfrac_to_molefrac(massfracs,mol_wt,molefracs)
-   massfracs::Array{Float64,1} : Vector of massfractions
-   mol_wt::Array{Float64,1} : Vector of molecular weights
"""
function massfrac_to_molefrac(massfracs::T,mol_wt::T) where T <: Array{Float64,1}
    moles = massfracs ./ mol_wt        
    return moles/sum(moles)
end


"""
Function to convert molefractions to massfractions
#   Usage:
    molefrac_to_massfrac!(massfracs,molefracs,mol_wt)
-   molefracs::Array{Float64,1} : Vector of molefractions
-   mol_wt::Array{Float64,1} : Vector of molecular weights
"""
function molefrac_to_massfrac!(massfracs::T, molefracs::T,mol_wt::T) where T <: Array{Float64,1}    
    wt = molefracs .* mol_wt    
    massfracs .= wt/sum(wt)    
end


"""
Average molecular weight calculation
"""
average_molwt(molfracs::T,molwt::T) where T = sum(molfracs .* molwt)
"""
density calculation
"""
density(avg_molwt::T1, T::T1, p::T1) where T1 = p*avg_molwt/(R*T)
density(molefrac::Array{T1,1},molwt::Array{T1,1}, T::T1, p::T1) where T1 = p*average_molwt(molefrac,molwt)/(R*T)


"""
Function to get floating point value from an  xml tag
#   Usage:
    get_value_from_xml(xmlroot,tag)
-   xmlroot::XMLElement : The xml root node
-   tag::String : name of the tag
The function returns the floating point value contained in the tag
"""
function get_value_from_xml(xmlroot::XMLElement,tag::String)
    tag_node = get_elements_by_tagname(xmlroot,tag)
    if is_elementnode(tag_node[1])
        tag_content = String(content(tag_node[1]))
        return parse(Float64,tag_content)
    else
        throw(error("Unable to find the node $tag in the xmlroot\n"))
    end    
end

"""
Function to get the file path from an xml tag
#   Usage:
    get_text_from_xml(xmlroot,tag)
-   xmlroot::XMLElement : The xml root node
-   tag::String : name of the tag
The function returns string of file path
"""
function get_text_from_xml(xmlroot::XMLElement,tag::String)
    tag_node = get_elements_by_tagname(xmlroot,tag)
    if is_elementnode(tag_node[1])
        return String(strip(content(tag_node[1])))
    else
        throw(error("Unable to find the node $tag in the xmlroot\n"))
    end    
end

"""
Function to get vector of strings contained in a tag.
Ideal for getting the species names 
#   Usage:
    get_collection_from_xml(xmlroot,tag)
-   xmlroot::XMLElement : The xml root node
-   tag::String : name of the tag
The function returns vector of strings. 

"""    
function get_collection_from_xml(xmlroot::XMLElement,tag::String)
    tag_node = get_elements_by_tagname(xmlroot,tag)   
    if is_elementnode(tag_node[1])
        return map(x->String(uppercase(x)),split(content(tag_node[1])))
    else
        throw(error("Unable to find the node $tag in the xmlroot\n"))
    end
end
    
"""
Function to get the mole fractions of species
#   Usage:
    get_molefraction_from_xml(xmlroot,molwt,gasphase)
-   xmlroot::XMLElement : The xml root node
-   molwt::Array{Float64} :  Vector of molecular weights
"""
function get_molefraction_from_xml(xmlroot::XMLElement, molwt::Array{Float64},gasphase::Array{String,1})
    tag_molefraction = get_elements_by_tagname(xmlroot,"molefractions")
    tag_massfraction = get_elements_by_tagname(xmlroot,"massfractions")
    local mass_fracs = Array{Float64,1}()
    local mole_fracs = Array{Float64,1}()
    if length(tag_massfraction) > 0
        massfraction = content(tag_massfraction[1])
        mass_fracs = parse_composition(massfraction,gasphase)
        mole_fracs = massfrac_to_molefrac(mass_fracs,molwt)
        return mole_fracs
    end    

    if length(tag_molefraction) > 0
        molefraction = content(tag_molefraction[1])
        mole_fracs = parse_composition(molefraction,gasphase)
        return mole_fracs
    end
end


function create_header(file_stream, args...)
    for i in eachindex(args)
        for k in args[i]
            @printf(file_stream,"%10s\t",k)
        end
    end
    @printf(file_stream,"\n")
end

function write_to_file(file_stream, args...)
    for i in eachindex(args)
        for k in args[i]
            @printf(file_stream,"%.4e\t",k)
        end
    end
    @printf(file_stream,"\n")
end


end # end of module