check_memory <- function(focal_data_cov, focal_data_mean, analog_pool, n_analog_use) {
    # get memory needed for all inputs in GB
    total_memory <- as.numeric(Sys.info()["mem_size"]) * 1e-9

    used_memory <- object.size(focal_data_cov) * 1e-9 +
        object.size(focal_data_mean) * 1e-9 +
        object.size(analog_pool) * 1e-9

    nprocesses <- parallel::detectCores()

    necessary_memory <- used_memory * nprocesses +
        sizeof(4) * 7 * nrow(focal_data_mean) * n_analog_use * 1e-9

    return(list(necessary_memory = necessary_memory, total_memory = total_memory))
}

create_outer_extents <- function(focal_data_coords, min_n_tiles) {
    # Get the minimum and maximum x and y coordinates
    min_x <- min(focal_data_coords$x)
    max_x <- max(focal_data_coords$x)
    min_y <- min(focal_data_coords$y)
    max_y <- max(focal_data_coords$y)

    # Calculate the width and height of the grid
    width <- abs(max_x - min_x)
    height <- abs(max_y - min_y)

    # Calculate the number of tiles in the x and y directions
    n_tiles_x <- ceiling(sqrt(min_n_tiles * width / height))
    n_tiles_y <- ceiling(sqrt(min_n_tiles * height / width))

    n_tiles <- n_tiles_x * n_tiles_y

    # Calculate the width and height of each tile
    tile_width <- width / n_tiles_x
    tile_height <- height / n_tiles_y

    # Create a data frame to store the extents of each tile
    tile_extents <- data.frame(
        tile_id = 1:n_tiles,
        min_x = rep(0, n_tiles),
        max_x = rep(0, n_tiles),
        min_y = rep(0, n_tiles),
        max_y = rep(0, n_tiles)
    )

    # Calculate the extents of each tile
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
    # Find the tile that contains the point (x, y)
    for (i in 1:nrow(tile_extents)) {
        if (x >= tile_extents[i, "min_x"] && x <= tile_extents[i, "max_x"] &&
            y >= tile_extents[i, "min_y"] && y <= tile_extents[i, "max_y"]) {
            return(as.integer(tile_extents[i, "tile_id"]))
        }
    }
}

spatial_partition <- function(focal_data_cov, focal_data_mean, analog_pool, n_analog_use) {
    # get memory needed for all inputs in GB
    memory <- check_memory(focal_data_cov, focal_data_mean, analog_pool, n_analog_use)
    necessary_memory <- memory$necessary_memory
    total_memory <- memory$total_memory

    if (necessary_memory > total_memory) {
        # calculate how many tiles we need to split the data into
        min_n_tiles <- ceiling(necessary_memory / total_memory)

        # Split the data into n_tiles using a spatial partitioning algorithm
        # Get the coordinates of the focal data
        focal_data_coords <- data.frame(
            x = focal_data_mean$x,
            y = focal_data_mean$y
        )
        tile_extents <- create_outer_extents(focal_data_coords, min_n_tiles)

        # Create a data frame to store the tile ID for each point
        tile_id <- data.frame(
            x = focal_data_coords$x,
            y = focal_data_coords$y,
            tile_id = rep(0, nrow(focal_data_coords))
        )

        # Assign each point to a tile
        for (i in 1:nrow(focal_data_coords)) {
            x <- focal_data_coords[i, "x"]
            y <- focal_data_coords[i, "y"]
            tile_id[i, "tile_id"] <- find_tile_id(x, y, tile_extents)
        }

        n_tiles <- max(tile_extents$tile_id)
        print(paste("Splitting data into", n_tiles, "tiles"))

        return(list(tile_id = tile_id, tile_extents = tile_extents))
    } else {
        print("Data fits in memory")
        return(list(tile_id = NULL, tile_extents = NULL))
    }
}

filter_analog_pool <- function(analog_pool, tile_extents, tile, max_dist) {
    # Filter the analog pool to only include points within max_dist of the focal points in the tile
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

    filtered_analog_pool <- subset(
        analog_pool,
        y >= min_y &
            y <= max_y &
            x >= min_x &
            x <= max_x
    )

    return(filtered_analog_pool)
}

filter_cov_data <- function(focal_data_cov, tile_id, tile) {
    # Filter the covariance data to only include points in the tile
    focal_data_cov_tile <- lapply(focal_data_cov, function(cov) {
        subset(
            cov,
            x >= tile_id[tile, "min_x"] &
                x <= tile_id[tile, "max_x"] &
                y >= tile_id[tile, "min_y"] &
                y <= tile_id[tile, "max_y"]
        )
    })

    return(focal_data_cov_tile)
}

filter_mean_data <- function(focal_data_mean, tile_id, tile) {
    # Filter the mean data to only include points in the tile
    focal_data_mean_tile <- subset(
        focal_data_mean,
        x >= tile_id[tile, "min_x"] &
            x <= tile_id[tile, "max_x"] &
            y >= tile_id[tile, "min_y"] &
            y <= tile_id[tile, "max_y"]
    )

    return(focal_data_mean_tile)
}
