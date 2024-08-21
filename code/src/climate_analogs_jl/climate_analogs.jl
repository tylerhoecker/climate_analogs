module ClimateAnalogs

# Function to load all the *fn.jl files
function load_functions()
    # Get the path to the directory containing the *fn.jl files
    script_dir = joinpath(dirname(@__FILE__), "..", "climate_analogs_jl")
    fn_files = filter(x -> endswith(x, r"fns?\.jl"), readdir(script_dir))
    
    # Load each *fn.jl file
    for file in fn_files
        include(joinpath(script_dir, file))
    end
end

# Load all the functions when the module is imported
load_functions()

export all_functions

end # module