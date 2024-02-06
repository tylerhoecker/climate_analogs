# Veg extraction and vote fn ---------------------------------------------------
veg_vote_fn <- function(analog_data, impact_path, n_analog_keep){
  
  impact_data <- rast(impact_path) 

  result <- analog_data |> 
    mutate(focal_id = as.numeric(focal_id)) |>
    # Extract Landfire 2020 BPS code
    mutate(bps = terra::extract(impact_data, 
                                as.matrix(analog_data[,c('a_x','a_y')]))[[names(impact_data)]],
           bps_true = terra::extract(impact_data, 
                                as.matrix(analog_data[,c('f_x','f_y')]))[[names(impact_data)]]) |> 
    # Turn distances into weights
    group_by(focal_id) |>
    mutate(weight = scales::rescale(dist_m, to = c(1,0))) |>
    # For each unique BPS code...
    group_by(focal_id, bps) |> 
    summarize(# Preserve coordinates
              x = first(f_x), y = first(f_y),
              # For validation, 'true' BPS code
              bps_true = first(bps_true),
              # Weighted vote: sum of weights for BPS code
              wt_vote = sum(weight),
              # Raw vote: count of analogs w/ this BPS code
              raw_vote = n(),
              raw_prop = raw_vote/n_analog_keep,
              # Statistics for each BPS code
              mean_sigma = mean(sigma),
              min_dist = min(dist_m),
              mean_dist = mean(dist_m)) |> 
    # This sequence saves the rows for each focal point, one for each BPS class,
    # that account for 90% of the (unweighted) analogs 
    arrange(desc(raw_prop), .by_group = T) |> 
    filter(
      cumsum(
        raw_prop == accumulate(
          raw_prop, ~ ifelse(.x <= .90, .x + .y, .y))) == 1) |> 
    mutate(n_bps_90 = length(bps)) 
    # |>
    # # Sort by the weighted votes for BPS classes
    # arrange(desc(wt_vote), .by_group = T) |> 
    # mutate(bps_rank = paste0('bps_',1:length(bps))) |>
    # ungroup() |>
    # tidyr::pivot_wider(names_from = bps_rank, values_from = bps)
    # And save the 
    # slice_head(n = 3) %>%  
    # summarise(x = first(x), y = first(y),
    #           bps1 = bps[1],
    #           bps2 = bps[2],
    #           bps3 = bps[3],
    #           win_prop = first(wt_vote),
    #           win_vote_raw = first(win_vote_raw),
    #           n_bps_90 = first(n_bps),
    #           mean_sigma = first(mean_sigma),
    #           min_dist = first(min_dist),
    #           mean_dist = first(mean_dist))  
      
  return(result)
}

  
saveRDS(veg_votes_df, 'wa_veg_proj_fullBPS.Rds')
veg_votes_df <- readRDS('wa_veg_proj_fullBPS.Rds')
#-------------------------------------------------------------------------------
# Rasterize and write out the TIFFs
bps_key <- read_csv('../data/bps_simplification.csv')

list('current_bps','bps1','bps2','bps3','mean_sigma_100','dist_bps_km','win_prop_raw','n_bps') %>%  #'current_bps','bps1','bps2','bps3','win_prop_raw','n_analogs','n_bps','mean_sigma_100') %>%  
  walk(function(band) {
    
    print(paste0('Rasterizing ',band,'...'))
    
    r <- veg_votes_df %>% 
      select(x, y, band) %>% 
      rast(., type = 'xyz', crs = crs(lf_bps)) 
    
    # if(band %in% c('current_bps','bps1','bps2','bps3')){
    #   r <- classify(r, cbind(bps_key$VALUE,bps_key$veg_num))
    # }
    # 
    r <- project(r, 'EPSG:4326', method = 'near')
    
    writeRaster(r, filename = paste0('../output/gee_InfN/fullBPS/wa_',band,'.tiff'), overwrite = T)
    
    return(r)
  }) 




  