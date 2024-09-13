# great circle distance between two points
# great circle distance between two points
haversine_formula(angle::AbstractFloat) = sin(angle / 2.0f0)^2.0f0
function great_circle_distance(
	lat1::T,
	lon1::T,
	lat2::T,
	lon2::T,
) where T <: AbstractFloat
	# Convert latitude and longitude from degrees to radians
	lat1 = lat1 * pi / 180.0
	lon1 = lon1 * pi / 180.0
	lat2 = lat2 * pi / 180.0
	lon2 = lon2 * pi / 180.0

	# Calculate the change in latitude and longitude
	dlat = lat2 - lat1
	dlon = lon2 - lon1



	# Calculate the great circle distance
	a = haversine_formula(dlat) + cos(lat1) * cos(lat2) * haversine_formula(dlon)
	c = 2 * asin(sqrt(a))
	R = 6378.137  # Earth's radius in km
	return R * c
end

function great_circle_distance(
	lat1::Vector{T},
	lon1::Vector{T},
	lat2::T,
	lon2::T,
) where T <: AbstractFloat
	# Convert latitude and longitude from degrees to radians
	pi32 = Float32(pi)
	lat1 = lat1 .* pi32 / 180.0f0
	lon1 = lon1 .* pi32 / 180.0f0
	lat2 = lat2 * pi32 / 180.0f0
	lon2 = lon2 * pi32 / 180.0f0

	# Calculate the change in latitude and longitude
	dlat = lat2 .- lat1
	dlon = lon2 .- lon1



	# Calculate the great circle distance
	a = haversine_formula.(dlat) .+ cos.(lat1) .* cos(lat2) .* haversine_formula.(dlon)
	c = 2.0f0 .* asin.(sqrt.(a))
	R = 6371.0f0  # Earth's radius in km
	return R .* c
end

function max_distance_coordinates(
	lat::AbstractFloat,
	lon::AbstractFloat,
	max_dist::Union{AbstractFloat, Integer})
	# Earth's radius in meters
	R = 6378.137

	# Calculate the change in latitude
	delta_lat = max_dist / R
	new_lat_north::Float32 = lat + (delta_lat * (180 / pi))
	new_lat_south::Float32 = lat - (delta_lat * (180 / pi))

	# Calculate the change in longitude
	delta_lon = max_dist / (R * cos(pi * lat / 180))
	new_lon_east::Float32 = lon + (delta_lon * (180 / pi))
	new_lon_west::Float32 = lon - (delta_lon * (180 / pi))

	# Return the new coordinates as a dictionary
	return Dict("north" => [new_lat_north],
		"south" => [new_lat_south],
		"east" => [new_lon_east],
		"west" => [new_lon_west])
end
