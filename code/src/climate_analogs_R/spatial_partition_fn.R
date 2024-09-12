library(dplyr)

NoMemory <- function(msg) {
    structure(list(message = msg), class = "NoMemory")
}

check_memory <- function(focal_data_cov, focal_data_mean, analog_pool, n_analog_use) {
    if (Sys.getenv("SLURM_CPUS_ON_NODE") != "") {
        if (Sys.getenv("SLURM_MEM_PER_NODE") != "") {
            total_memory <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")) * 1e-3
        } else {
            total_memory <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")) * as.numeric(Sys.getenv("SLURM_MEM_PER_CPU")) * 1e-3
        }
    } else {
        total_memory <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)) * 1e-6
    }

    used_memory <- object.size(focal_data_cov) + object.size(focal_data_mean) + object.size(analog_pool)
    sampling_memory <- 4 * 7 * nrow(focal_data_mean)
    nprocesses <- parallel::detectCores() * 2
    environment_memory <- (used_memory + sampling_memory) * nprocesses
    output_memory <- 4 * 7 * nrow(focal_data_mean) * n_analog_use
    necessary_memory <- (environment_memory + output_memory) * 1e-9 + 30

    if (total_memory <= 30) {
        stop(NoMemory("Does not contain the minimal memory needed for overhead (30 GB)!"))
    }

    return(list(necessary_memory = necessary_memory, total_memory = total_memory))
}

create_outer_extents <- function(focal_data_coords, min_n_tiles) {
    min_x <- min(focal_data_coords$x)
    max_x <- max(focal_data_coords$x)
    min_y <- min(focal_data_coords$y)
    max_y <- max(focal_data_coords$y)

    width <- abs(max_x - min_x)
    height <- abs(max_y - min_y)

    n_tiles_x <- ceiling(sqrt(min_n_tiles * width / height))
    n_tiles_y <- ceiling(sqrt(min_n_tiles * height / width))

    n_tiles <- n_tiles_x * n_tiles_y
    tile_width <- width / n_tiles_x
    tile_height <- height / n_tiles_y

    tile_extents <- data.frame(
        tile_id = 1:n_tiles,
        min_x = numeric(n_tiles),
        max_x = numeric(n_tiles),
        min_y = numeric(n_tiles),
        max_y = numeric(n_tiles)
    )

    for (i in 1:n_tiles_x) {
        for (j in 1:n_tiles_y) {
            tile_id <- (i - 1) * n_tiles_y + j
            tile_extents[tile_id, "min_x"] <- min_x + (i - 1) * tile_width
            tile_extents[tile_id, "max_x"] <- min_x + i * tile_width
            tile_extents[tile_id, "min_y"] <- min_y + (j - 1) * tile_height
            tile_extents[tile_id, "max_y"] <- min_y + j * tile_height
        }
    }

    return(tile_extents)
}

find_tile_id <- function(x, y, tile_extents) {
    for (i in 1:nrow(tile_extents)) {
        if (x >= tile_extents[i, "min_x"] && x <= tile_extents[i, "max_x"] &&
            y >= tile_extents[i, "min_y"] && y <= tile_extents[i, "max_y"]) {
            return(as.integer(tile_extents[i, "tile_id"]))
        }
    }
}

spatial_partition <- function(focal_data_cov, focal_data_mean, analog_pool, n_analog_use) {
    memory_info <- check_memory(focal_data_cov, focal_data_mean, analog_pool, n_analog_use)
    necessary_memory <- memory_info$necessary_memory
    total_memory <- memory_info$total_memory

    if (necessary_memory > total_memory) {
        cat("Data does not fit in memory: necessary:", necessary_memory, ", total_memory:", total_memory, "\n")

        min_n_tiles <- ceiling(necessary_memory / total_memory)
        focal_data_coords <- data.frame(x = focal_data_mean$x, y = focal_data_mean$y)
        tile_extents <- create_outer_extents(focal_data_coords, min_n_tiles)

        tile_id <- data.frame(
            x = focal_data_coords$x,
            y = focal_data_coords$y,
            tile_id = sapply(1:nrow(focal_data_coords), function(i) {
                find_tile_id(focal_data_coords[i, "x"], focal_data_coords[i, "y"], tile_extents)
            })
        )

        n_tiles <- max(tile_extents$tile_id)
        cat("Splitting data into", n_tiles, "tiles\n")

        return(list(tile_id = tile_id, tile_extents = tile_extents))
    } else {
        cat("Data fits in memory: necessary:", necessary_memory, ", total_memory:", total_memory, "\n")
        return(list(tile_id = NULL, tile_extents = NULL))
    }
}

filter_analog_pool <- function(analog_pool, tile_extents, tile, max_dist) {
    min_x <- max_distance_coordinates(tile_extents[tile, "min_y"], tile_extents[tile, "min_x"], max_dist)$west[2]
    max_x <- max_distance_coordinates(tile_extents[tile, "max_y"], tile_extents[tile, "max_x"], max_dist)$east[2]
    min_y <- max_distance_coordinates(tile_extents[tile, "min_y"], tile_extents[tile, "min_x"], max_dist)$south[1]
    max_y <- max_distance_coordinates(tile_extents[tile, "max_y"], tile_extents[tile, "max_x"], max_dist)$north[1]

    analog_extents <- data.frame(
        min_x = min_x,
        max_x = max_x,
        min_y = min_y,
        max_y = max_y
    )

    filtered_analog_pool <- analog_pool %>%
        filter(y >= min_y & y <= max_y & x >= min_x & x <= max_x)

    return(filtered_analog_pool)
}

filter_cov_data <- function(focal_data_cov, tile_id, tile) {
    focal_data_cov_tile <- focal_data_cov
    for (i in 1:length(focal_data_cov)) {
        focal_data_cov_tile[[i]] <- focal_data_cov[[i]] %>%
            filter(x >= tile_id[tile, "min_x"] & x <= tile_id[tile, "max_x"] &
                y >= tile_id[tile, "min_y"] & y <= tile_id[tile, "max_y"])
    }

    return(focal_data_cov_tile)
}

filter_mean_data <- function(focal_data_mean, tile_id, tile) {
    focal_data_mean_tile <- focal_data_mean %>%
        filter(x >= tile_id[tile, "min_x"] & x <= tile_id[tile, "max_x"] &
            y >= tile_id[tile, "min_y"] & y <= tile_id[tile, "max_y"])

    return(focal_data_mean_tile)
}
