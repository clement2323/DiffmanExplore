
get_centroid_carreau <- function(data, var_carreau){
  
  taille_carreau <- readr::parse_number(stringr::str_extract(data[[var_carreau]][1], "RES[0-9]*m"))
  centroides <- data %>%
    select(carreau = all_of(var_carreau)) %>%
    mutate(
      ll_coord_x = readr::parse_number(stringr::str_extract(carreau, "E[0-9]*$")),
      ll_coord_y = readr::parse_number(stringr::str_extract(carreau, "N[0-9]*")),
      centroid_x = ll_coord_x + taille_carreau/2,
      centroid_y = ll_coord_y + taille_carreau/2,
      crs = readr::parse_number(stringr::str_extract(carreau, "CRS[0-9]*"))
    )
  return(list(df = centroides, epsg=centroides$crs[1], taille=taille_carreau))
}

carreaux_to_polygon <- function(data, var_carreau){
  
  require(btb)
  
  centroides_l <- get_centroid_carreau(data, var_carreau)
  
  centroides_l$df %>%
    select(carreau, x=centroid_x, y=centroid_y) %>%
    btb::dfToGrid(sEPSG = centroides_l$epsg, iCellSize = centroides_l$taille)
  
}

# pour préparer les géométries
preparer_geom <- function(input_dt){
  
  dt <- copy(input_dt)
  liste_carreau <- unique(dt$z2)
  
  polygone_carreau <-
    carreaux_to_polygon(data.frame(carreau = liste_carreau), var_carreau = "carreau") %>%
    sf::st_transform(crs = 4326)
  
  geom_z1 <- communes %>%
    filter(code %in% unique(c(dt$z1))) %>%
    select(code) %>%
    rename(z1 = code)
  
  geom_z2 <- polygone_carreau %>%
    rename(z2 = carreau) %>%
    select(-x,-y)
  
  out <- list(geom_z1 = geom_z1, geom_z2 = geom_z2)
  
  return(out)
}


search_diff_pb_one_resolution <- function(resolution,input_dt_grid){
  
  input_dt_resolution <- input_dt_grid %>% 
    select(z1, all_of(resolution),nb_obs) %>% 
    rename(z2 = all_of(resolution))
  
  input_dt_resolution <- input_dt_resolution[,.(nb_obs=sum(nb_obs)),by =.(z1,z2)]
  
  # controle simple, nombre d'ibnter carreaux x communes avec - de 11 menages
  
  if(nrow(input_dt_resolution[nb_obs < 11]) == 0) stop("Aucune intersection en dessous du seuil")
  
  save_name <- unlist(strsplit(resolution,"_"))[4]
  
  all_component_risk_extraction(
    input_dt_resolution,
    threshold = 11,
    max_agregate_size = 15,
    save_dir = paste0("diff_info_",save_name)
  )
}
