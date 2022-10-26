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
  #input_df <- df
  df <- copy(input_df)
  df <- df[,.(nb_obs = sum(nb_obs)),by = .(z1,z2)]# doublon
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
  # input_df <- df[substr(z1,1,2) == "22"]
  
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


# construire la matrice de croisement à partir de la link table
construire_m_crois <- function(link_table){
  
  df <- link_table[ ,
                    .(
                      nb_obs_from = sum(nb_obs_from),
                      nb_obs_to = sum(nb_obs_to),
                      id_comp = unique(id_comp),
                      z2 = paste0(z2,collapse = "-"),
                      n_z2 = unique(n_z2)
                    ),
                    by = .(from,to)
  ]
  
  ltable <-long_table(df)
  
  ltable$z1  <- as.factor(ltable$z1)
  ltable$z2  <- as.factor(ltable$z2)
  
  # variable" when running R CMD check
  # Matrice z1Xz2
  m_crois <- Matrix::sparseMatrix(i=as.numeric(ltable$z1),
                                  j=as.numeric(ltable$z2),
                                  x=ltable$nb_obs,
                                  dimnames=list(levels(ltable$z1),levels(ltable$z2)))
  
  m_crois 
}
# a partir d'une liste de commune (en sortie ou non de diffman) sort toutes les informations permettant d'évaluer si un problèùe de différenciaition existe sur cette zone
# en sortie, 2 booleen diff_intern issue (resp external), alertant sur le type de pb de différenciation
# puis les tables avec les intersections carreaux x commune pour les carreaux à cheval sur la zone en question (définie par l'union des communes en paramètre dans liste_z1)
# TO DO ajouter l'ensemble des carreaux inclus totalement ou partiellement dans la zone (nécessite d'avoir les données en entrée ou de sortir lors de la construction de la linktable une liste avec link table d'un côté et de l'autre la table des carreaux inclus totalement )
return_info <- function(liste_z1,link_table, threshold){
  
  #liste_z1
  
  # on veut pour la zone en question sortir tous les carreaux problématiques (dire si on a différenciation externe ou interne)
  dataf <- long_table(link_table) 
  z2_target <- dataf[z1 %in% liste_z1,]$z2 # les carreaux impactés, par def de link table
  
  # sortir les z2 à cheval  sur 1 commune de la zone et une commune externe
  at_risk_crossing <- unique(dataf[z2 %in% z2_target &  ! z1 %in% liste_z1,"z2"])
  
  # je les récupère dans la table initiame
  diff_table <- dataf[z2 %in% at_risk_crossing$z2 ]
  external_diff_table <-diff_table[ !z1 %in% liste_z1]
  internal_diff_table <-diff_table[ z1 %in% liste_z1]
  
  external_diff_issue <- internal_diff_issue <- FALSE
  
  # Ok + qu'à sommer sur nb_obs pour savoir si le seuil est respecté ou non dans la diff interne ou externe
  if (sum(internal_diff_table$nb_obs) < threshold) internal_diff_issue <- TRUE
  if (sum(external_diff_table$nb_obs) < threshold) external_diff_issue <- TRUE
  
  # Il faut sortir les carreaux au bord de la zone
  
  out <- list(internal_diff_table = internal_diff_table,
              external_diff_table = external_diff_table,
              internal_diff_issue = internal_diff_issue, 
              external_diff_issue = external_diff_issue
  )
  return(out)
  # il faut maintenant sortir la partie des carreaux non incluse dans la zone de commune pour checker la diff externe
  # et isoler celle incluse pour la diffinterne (easy)
  
}

# fonction qui prend en entrée m_Crois, qui opère les focnctions de réduction de graph dessus et qui retourne in fine la liste des zones, (union d'élément de z1 ) à risque 
# TO DO refaire une fonction qui récupère la matrice m_Crois après pour mesurer la force simplificatrice (nb composantes ) checker)
# penser à regarder les codes c++

sortir_liste_zones_a_risque <- function(
    link_table,
    max_agregate_size = 15,
    save_file = NULL,
    simplify = TRUE,
    verbose = TRUE,
    threshold = 11
  ){
    
    m_crois <- construire_m_crois(link_table)
    # threshold = 7; max_agregate_size = 15;save_file = NULL; simplify = TRUE; verbose = TRUE
    if(simplify){ #one can choose to skip these steps of graph reduction if desired
      if(verbose) message("< --- Merging method 1 --- > ")
      m_crois <- diffman:::agregate(m_crois, threshold, methode = "m1", verbose = verbose)
      
      #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m1_",save_file,".RDS"))
      
      if(verbose) message("< --- Merging methods 1 and 2 --- >")
      m_crois <- diffman:::agregate(m_crois, threshold, methode = "both", verbose = verbose)
      
      # warnings()
      #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m2_",save_file,".RDS"))
      
      if(sum(dim(m_crois)==0)>0) {
        message("No differentiation problems detected !")
        return(NULL)
      }
      
      if(verbose) message("< --- Splitting the graph --- >") # à regarder
      l_decomp <- diffman:::decompose_m_crois(m_crois, max_agregate_size)
      
    }else{
      l_decomp <- diffman:::comp_connexe_list(m_crois)
    }
    
    if(verbose) message("< --- Exhaustive search of differentiation problems --- >")
    # sauvegarder l_decomp pourrait être intéressant
    l_ag <- diffman:::search_diff_agregate(l_decomp, threshold, max_agregate_size) 
    # dans fct_search_diff test composante connexe par composante connexe
    l_ag <- diffman:::desagregate_list(l_ag) 
    
    return(l_ag)
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
