# Function to load all the *fn.jl files
load_functions <- function() {
    # Get the path to the directory containing the *fn.jl files
    script_dir <- file.path(getwd(), "code/src/climate_analogs_R")
    fn_files <- list.files(script_dir, pattern = "fns?\\.R$", full.names = TRUE)
    fn_files <- c(fn_files, file.path(script_dir, "as.data.table.R"))

    # Load each *fn.jl file
    for (file in fn_files) {
        source(file)
    }
}

# Load all the functions when the module is imported
load_functions()

# Export all the functions
all_functions <- ls()
