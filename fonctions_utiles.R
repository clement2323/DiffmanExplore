
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
  # input_dt <- situation_table
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
  # resolution <- "id_carreau_niv_64km"
  # resolution <- "id_carreau_niv_32km"
  # resolution <- "id_carreau_niv_16km"
  # resolution <- "id_carreau_niv_8km"
  # resolution <- "id_carreau_niv_4km"
  # resolution <- "id_carreau_niv_2km"
  input_dt_grid <- copy(input_dt_grid)
  
  input_dt_resolution <- input_dt_grid[,c("z1",resolution,"nb_obs"),with = FALSE]
  colnames(input_dt_resolution)[colnames(input_dt_resolution) == resolution]<- "z2"
  
  
  input_dt_resolution <- input_dt_resolution[,.(nb_obs=sum(nb_obs)),by =.(z1,z2)]
  
  # j'échange les rôles de z1 et z2 si il y a un petit nombre de z2
  if(length(unique(input_dt_resolution$z1))>length(unique(input_dt_resolution$z2))){
    print("switch des rôles entre z1 et z2")
    input_dt_resolution[,z3:=z2];input_dt_resolution[,z2:=z1];input_dt_resolution[,`:=`(z1=z3,z3=NULL)]
  }
    
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
# readRDS("diff_info_16km/res_1.RDS")

build_compl <- function(checked_area,list_z1_compo){
  paste0(setdiff(list_z1_compo,unlist(strsplit(checked_area,"-"))),collapse="-")
}



# A partir de la liste des z1 dans une composante connexe et des problemes de différenciation extrait je construis une liste exhaustive 
# de zone interne (desfois une intersectionn n'apparait qu'en ionternal ou en external)
build_complete_internal_table <- function(comp_diff_info,list_z1_compo){
  
  external_area <- unique(comp_diff_info[type_diff == "external"]$checked_area)
  internal_area <- unique(comp_diff_info[type_diff == "internal"]$checked_area)

  # ça se passe bien sur la grosse compo -> 1 min # library(pbapply)
  internal_area_expected <- sapply(external_area,build_compl,list_z1_compo = list_z1_compo)
  
  missing_internal_area <- setdiff(internal_area_expected,internal_area)
  
  if (length(missing_internal_area) !=0){
    
    missing_internal_area <- internal_area_expected[internal_area_expected %in% missing_internal_area]
    
    # fait tout cracher relativement long pour du data table sur la grosse compo
    
    # je crée  la table de passage  external to missing internal
    external_to_internal <- data.table(external_checked_area = as.factor(names(missing_internal_area)), 
                                       internal_checked_area = as.factor(unname(missing_internal_area)),
                                       type_diff = "internal"
                                       )
    
    

    # avec les facteurs -> ultra rapide
    tmp<-comp_diff_info[,c("z1","z2","nb_obs","checked_area","id_comp")] 
    tmp[,checked_area := as.factor(checked_area)]
    
    # setkey(external_to_internal,"external_checked_area") # crash faire du as.factor ? les chaines sont trop longue ?
    # setkey(comp_diff_info,"checked_area")
    
    # controle 
    new_internal <- merge(tmp,external_to_internal,by.x = "checked_area" ,by.y = "external_checked_area",nomatch = 0)
    new_internal$checked_area <- NULL
    colnames(new_internal)[colnames(new_internal)== "internal_checked_area"] <- "checked_area"
    
    complete_internal_diff_info <-rbind(comp_diff_info[type_diff == "internal"],new_internal)
    
  }else{
    complete_internal_diff_info <- comp_diff_info[type_diff == "internal"]
  }
  return(complete_internal_diff_info)
}


# diffman ne sort pas touts les types de differenciation par exemple la commune 26267 ne sort pas mais son complémentaire oui donc pas grave..
# fonction proteger compo qui balai toutes les cheqck area blanchi de aprt et autres de la ckeked area (dans son complémentaire aussi) on blanchit jusqu'à ce que la différence dépasse 11

protect_component <- function(num_comp,global_diff_info,data_rp,z2_to_tag,threshold =11){ 
  # comp_to_nb_com
  # num_comp <- 997
  # chevauchant => pas de diff_interne possible mais une diff externe possible ici on traite ce sous cas ici
  
  dt <- copy(data_rp)
  list_z1_compo<- compo$z1[compo$id_comp == num_comp]
  input_dt <- clean_init_dt(dt[z1 %in% list_z1_compo])
  
  fully_included_z2<- diffman:::prepare_data(input_dt)$fully_included_z2
  z2_to_tag$full_incl <- z2_to_tag$z2 %in% unique(fully_included_z2$z2)
  comp_diff_info <- global_diff_info[id_comp == num_comp]
  comp_diff_info$checked_area_size <- NULL
  
  if(nrow(comp_diff_info)==0) return(NULL)  # Je build ici la table avectoutes les internal différences
  complete_internal_diff_info <- build_complete_internal_table(comp_diff_info,list_z1_compo)
  
  #### Et voici l'algorithme !!
  l <- split(complete_internal_diff_info,complete_internal_diff_info$checked_area)
  i <-1
  
  for(area_issue in l){
    # area_issue <- l[[1]]
    # area_issue <-   l$`15200`
    
    if (i%%200 == 0) print(i)
    i <- i+1
    
    z2_to_tag <- actualiser_z2_to_tag(area_issue,list_z1_compo,input_dt,z2_to_tag) 
     
  } # fin de la boucle 
  return(z2_to_tag)
}

# Je protège ici une zone à risque de différenciation obtenue via une checked area donnée
actualiser_z2_to_tag <- function(area_issue,list_z1_compo,input_dt,z2_to_tag){
    # threshold = 11
    z2_to_tag <- copy(z2_to_tag)
    nb_obs_at_risk <- sum(area_issue$nb_obs)
    nb_to_add <- threshold - nb_obs_at_risk
    
    # en amont récupérer 
    l <- calculer_zone_interne_externe(area_issue,list_z1_compo,z2_to_tag)
    z2_full_incl <- l$z2_full_incl
    z2_full_excl <- l$z2_full_excl
    
    # je réadapte le nb en fonction de ce qui a été blanchit déjà et je supprimeai les lignes dans la table après
    nb_to_add_full_incl <- nb_to_add - sum(z2_full_incl$nb_obs[z2_full_incl$tag])
    nb_to_add_full_excl <- nb_to_add - sum(z2_full_excl$nb_obs[z2_full_excl$tag])
    
    if(nb_to_add_full_incl >0 & nrow(z2_full_incl) > 0 ){
      
      z2_full_incl <- z2_full_incl[tag == FALSE]
      z2_full_incl[,cum_nb_obs := cumsum(nb_obs)]
      z2_full_incl[,higher := cum_nb_obs>=nb_to_add_full_incl]
      
      ind_sup <- first(which(z2_full_incl$higher))
      
      if (is.na(ind_sup)){ # pas assez d'obs pour proteger on blanchit tout même l'iontersection et on actualise z2 tot ag direct
        z2_to_mask_incl <- NULL
        z2_to_tag[z2 %in% input_dt[z1 %in% z1_in_area]$z2]$tag <- TRUE #on blanchit tout même l'intersection qu'opn récupère (la commune n'a pas assez de ménages )
      }else{ # fonctionnement normal
      ind_z2 <- 1:ind_sup
      z2_to_mask_incl <- z2_full_incl$z2[ind_z2]
      }
    }
    # RAJOUTER CONDITION  SUR LA TAILLE DE LA ZONE EXTERNE ICI POUR SAVOIR SI ONLA PROTEGE AUSSI
    if(nb_to_add_full_excl >0 & nrow(z2_full_excl)> 0 ){
      
      z2_full_excl <- z2_full_excl[tag == FALSE]
      z2_full_excl[,cum_nb_obs := cumsum(nb_obs)]
      z2_full_excl[,higher := cum_nb_obs>nb_to_add_full_excl]
      
      ind_sup <- first(which(z2_full_excl$higher)) 
      
      # ne pas oublier que les inférieurs au seuil seront nécessairement blanchis il faut juste en ajouter si ils ne suffisent pas (si il n'yen a pas  1 seul carreau suffira)
      
      if (is.na(ind_sup)){ # pas assez d'obs pour proteger on blanchit tout même l'iontersection et on actualise z2 tot ag direct
        z2_to_mask_excl <- NULL
        z2_to_tag[z2 %in% z2_full_excl$z2]$tag <- TRUE #on blanchit tout aussi
      }else{
      ind_z2 <- 1:ind_sup
      z2_to_mask_excl <- z2_full_excl$z2[ind_z2]
      }
    }
    # carreaux dans la zone (il faut blanchir jusqu'à 11-nbobs atrisk)
    
    z2_to_tag[z2 %in% c(z2_to_mask_incl,z2_to_mask_excl)]$tag <- TRUE
    
    return(z2_to_tag)
}



calculer_zone_interne_externe <- function(area_issue,list_z1_compo,z2_to_tag){
  
  # je sors l'ensemble des carreaux de la zone interne
  z1_in_area <-unlist(strsplit(unique(as.character(area_issue$checked_area)),"-"))
  z1_out_area <- setdiff(list_z1_compo,z1_in_area)  
  
  #On se limite au blanchiment des full_incl quand on regarde la differenciation interne
  z2_full_incl <-input_dt[z1 %in% z1_in_area & z2 %in% z2_to_tag[full_incl == TRUE]$z2] 
  # Rq : on a bien 1 commune pour 1 carreau ici par definition des z2 totalement inclus
  

  # On définit  z2_full_excl par le complémentaire et on fait passer la table au niveau carreau (cf intersection prises en compte)
  z2_full_excl <- 
    input_dt[!paste0(z1,z2) %in% paste0(z2_full_incl$z1,z2_full_incl$z2)]# la zone a risque
  
  # il faut récupérer le tag !!
  z2_full_incl <- merge(z2_full_incl,z2_to_tag,by ="z2")[,.(nb_obs = sum(nb_obs)), by = "z2"][order(nb_obs)]
  z2_full_excl <- merge(z2_full_excl,z2_to_tag,by ="z2")[,.(nb_obs = sum(nb_obs)), by = "z2"][order(nb_obs)]
  
  return(list(z2_full_incl = z2_full_incl, z2_full_excl = z2_full_excl))
}



# comp_to_nb_com
# num_comp <- 263

# z2_to_tag <- data.table(z2  = unique(input_dt$z2), tag  = FALSE, tag_sup = FALSE)
# z2_to_z2_sup <- emboitement[,c("id_carreau_petit","id_carreau_niv6")]
# colnames(z2_to_z2_sup)<- c("z2","z2_sup")
# z2_to_tag<- merge(z2_to_tag,z2_to_z2_sup,by = "z2")


# verifier que l'on ne peut reconstruire les zones interne avec les carreaux de niveau supérieur
# il suffit de faire une fusion entre les deux zones, si elle est non vide on a gagné sinon il faut blanchir au niveau au dessus 

# merge_z2_sup <-merge(z2_full_excl,z2_full_incl,by ="z2_sup")
# arrive beaucoup trop souvent..
# if(nrow(merge_z2_sup==0)){print(paste0("level sup issue for area ",i))}

