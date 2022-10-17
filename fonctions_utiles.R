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

nettoyer_base_init <- function(input_df){
  
  df <- copy(input_df)
  df <- df[,.(nb_obs=sum(nb_obs)),by = .(z1,z2)]# doublon
  ## Suppression des lignes hors format mettre en place des controles + suppression des zéros
  df <- df[z2 != "FR_unallocated"]
  df <- df[nb_obs != 0]
  
  df
}

# input df z1xz2 table areas with nb_obs
# return table with link between z1 (undirected graph )
build_link_table <- function(input_df){
  # donnees_rp <- readRDS("../data_rp.rds")
  # df <- setDT(donnees_rp)
  # input_df <- df[substr(z1,1,2) == "40"]
  
  df <- copy(input_df)
  
  df <- nettoyer_base_init(df)
  ## Je mets de côté les carreaux (z2) intersectant une et une seule commune (z1)
  df_z2 <- df[,.(nb_z1=length(z1)),by=.(z2)]
  df_z2_mono_z1 <- df[z2 %in% df_z2[nb_z1==1]$z2]
  df_z2_multi_z1 <- df[z2 %in% df_z2[nb_z1>1]$z2]
  
  
  link_table <- merge( 
                       df_z2_multi_z1,
                       df_z2_multi_z1,
                       by=c('z2'),
                       allow.cartesian = TRUE,
                       nomatch = 0,
                       suffix = c("_from","_to")
                       ) 
  
  colnames(link_table)[colnames(link_table)=="z1_from"] <- "from"
  colnames(link_table)[colnames(link_table)=="z1_to"] <- "to"
  
  link_table <- link_table[from != to]
  
  link_table$from_to_z2 <- apply(link_table[,c("from","to","z2")],1,function(x) paste0(sort(x),collapse= "-"))
  link_table <- link_table[!duplicated(from_to_z2)][,c("from","to","z2","nb_obs_from","nb_obs_to")]
  
  link_table[,n_z2 := length(z2) ,by =.(from,to)]
  
  retourner_graph(link_table)
}

retourner_graph<-function(link_table){
  # Construction du graph
  # https://igraph.org/r/#docs  joli !
  
  graph <- graph_from_data_frame(link_table,directed = FALSE)
  clust <- clusters(graph)
  nodes <- data.table(from = names(clust$membership), id_comp = clust$membership)
  
  # je regarde les clusters
  out <- nodes[link_table, on = "from"]
  
  out
}

# remettrela link table au niveau carreau, commune
long_table <- function(link_table){
  
  tab_from <- link_table[,c("from","z2","nb_obs_from","id_comp")] 
  colnames(tab_from) <- c("z1","z2","nb_obs","id_comp")
  
  tab_to <-  link_table[,c("to","z2","nb_obs_to","id_comp")] 
  colnames(tab_to) <- c("z1","z2","nb_obs","id_comp")
  
  out <- rbind(tab_from,tab_to)
  
  return(out)
  
}


# dessiner une situation donnée (situation = sous ensembe de la link table
# trace la situation à partir des géométries ciomplètes de z1 zet z2 en entrée : réalise l'intersection géométrique et &affiche le nb_obs

draw_situation <- function(situation_table,z2_to_nb_obs,geom_z1,geom_z2,liste_z1_to_color = NULL ,threshold = 11,save_name = NULL){
  
  # récupération du total par z2
  # z2_to_nb_obs <- long_table(link_table)[ ,.(nb_obs_z2 = sum(nb_obs)) ,by = .(z2)]
  
  # situation_table <- link_table[id_comp == 6] ;unique(c(situation_table$from,situation_table$to))
  # liste_z1_to_color <- c("01014")
  # liste_carreau <- situation_table$z2
  # polygone_carreau <-
  #   carreaux_to_polygon(data.frame(carreau = liste_carreau), var_carreau = "carreau") %>%
  #   st_transform(crs = 4326)
  # 
  # geom_z1 <- communes %>%
  #   filter(code %in% unique(c(situation_table$from,situation_table$to))) %>%
  #   select(code) %>%
  #   rename(z1 = code)
  # 
  # geom_z2 <- polygone_carreau %>%
  #   rename(z2 = carreau) %>%
  #   select(-x,-y)

  
  # recupération des intersections
  st_agr(geom_z1) = "constant"
  st_agr(geom_z2) = "constant"
  
  inter_carreau_commune <-st_intersection(
    geom_z1 %>% select(z1),
    geom_z2 %>% select(z2)
  )
  
  # je déconstruis la table pour avoir 
  # 
  # tab_from <- situation_table[,c("from","z2","nb_obs_from")] 
  # colnames(tab_from) <- c("z1","z2","nb_obs")
  # 
  # tab_to <-  situation_table[,c("to","z2","nb_obs_to")] 
  # colnames(tab_to) <- c("z1","z2","nb_obs")
  # 
  # tab_croisement <- rbind(tab_from,tab_to)
  
  tab_croisement<- long_table(situation_table)
  
  inter_carreau_commune <- merge(tab_croisement,inter_carreau_commune,by = c("z1","z2"))
  
  geom_z2 <- merge(geom_z2,z2_to_nb_obs,by ="z2",nomatch = 0)
  
  z1_fillColor <- with(geom_z1,ifelse(z1 %in% liste_z1_to_color,"orange","#3FC8FC"))
  
  highlightOptions_defaut <- highlightOptions(
    stroke = TRUE,
    weight = 6,
    color = "black",
    fillColor = "black",
    bringToFront = TRUE
  )
  
  
  carte <- 
    leaflet() %>% 
    addProviderTiles("GeoportailFrance.orthos") %>%  
    addPolygons(
      data = geom_z1,
      color = "#3FC8FC",
      fillColor = z1_fillColor,
      weight = 2,
      fillOpacity = 0.25,
      opacity = 1,
      label = geom_z1$z1
    ) %>% 
    addPolygons(
      data = geom_z2,
      color = "red",
      label = with(geom_z2, 
           sprintf(
             "<b> id z2 : </b> %s  <br/> <b> Number of observations : </b>  %s", ### c'est une d�finition de format qui vient du C
            z2, round(nb_obs_z2,1)
           ) %>% lapply(htmltools::HTML)
      ),
      weight = 2,
      fillOpacity = 0,
      opacity = 1,
      group ="z2 on two sides of one z1 area"
    ) %>% 
    addPolygons(
      data = inter_carreau_commune$geometry,
      color = ifelse(inter_carreau_commune$nb_obs < threshold,"red","#6E3DFF"),
      weight = 2,
      fillOpacity = 0.5,
      group = "intersections",
      highlightOptions = highlightOptions_defaut,
      label  =  with(inter_carreau_commune, 
                     sprintf(
                       "<b> id z1 : </b> %s  <br/> <b> id z2 : </b>  %s <br/> <b> Number of observations : </b>  %s", ### c'est une d�finition de format qui vient du C
                       z1, z2, round(nb_obs,1)
                     ) %>% lapply(htmltools::HTML)
      )
    ) %>% 
    addScaleBar(position="bottomright") %>% 
    hideGroup(c("z2 on two sides of one z1 area","intersections")) %>% 
    addLayersControl(
      overlayGroups = c("z2 on two sides of one z1 area","intersections"),
      options = layersControlOptions(collapsed = FALSE)
    ) %>% 
    addScaleBar(position="bottomright")
  
  if(!is.null(save_name)) htmlwidgets::saveWidget(carte, file=paste0(save_name,".html"),selfcontained = TRUE)
  
  carte
}
