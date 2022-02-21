module ReactionEngine
# Write your package code here.
include("Inspect.jl")
include("ThermoProbe.jl")
include("Plug.jl")
include("Cstr.jl")
include("Batch.jl")
include("Equil.jl")
include("Transport.jl")
include("Sens.jl")
export inspect, thermoprobe, cstr, plug, batch, equilibrate, transport_properties, global_sensitivity

"""
inspect(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
inspect(input_file::AbstractString, lib_dir::AbstractString) = Inspect.run(input_file, lib_dir)

"""
thermoprobe(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
thermoprobe(input_file::AbstractString, lib_dir::AbstractString ) = ThermoProbe.run(input_file,lib_dir)

"""
cstr(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
cstr(input_file::AbstractString, lib_dir::AbstractString) = Cstr.cstr(input_file,lib_dir)

"""
plug(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
plug(input_file::AbstractString, lib_dir::AbstractString) = Plug.plug(input_file,lib_dir)


"""
batch(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
batch(input_file::AbstractString, lib_dir::AbstractString) = Batch.batch(input_file, lib_dir)



"""
equilibrate(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
equilibrate(input_file::AbstractString, lib_dir::AbstractString) = Equil.equil(input_file, lib_dir)

"""
transport_properties(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The  directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
transport_properties(input_file::AbstractString, lib_dir::AbstractString) = Transport.run(input_file,lib_dir)


"""
global_sensitivity(input_file::AbstractString, lib_dir::AbstractString)
- input_file : Input xml file (String)
- lib_dir : The directory where the data files are present. 
            In the present version, the data files are present in the lib directory within the test folder. 
            The user must specify the relative path to the lib directory. The name of the lib directory need not be
            "lib"
"""
global_sensitivity(input_file::AbstractString, lib_dir::AbstractString) = Sens.gsa_xml(input_file,lib_dir)


end
