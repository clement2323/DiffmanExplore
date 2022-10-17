
install.packages("igraph")
install.packages("visNetwork")
installed.packages("tidygraph")
install.packages("ggraph")
install.packages("readr")
install.packages("btb")
install.packages("pbapply")


library(pbapply)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggplot2)
library(stringr)
library(sf)
library(leaflet)

source("Z:/carreaux_rp/R/fonctions_geo.R", encoding = "UTF-8")

com_2021 <- st_read("S:/PROD/creacartes/pd/fichiers-ihm/2021/franceentiere/commune_franceentiere_2021.shp")
data <- readRDS("X:/HAB-Traitement-Confidentialite/RP/Census_2021/inputs/pop_carreau_dmrg.rds")

# on a des doublons mais ok on verra dans un second temps
doublons <-data %>%
  select(depcom,carreau,nb_log) %>%
  filter(duplicated(paste0(depcom,carreau))) %>%
  left_join(data %>% select(depcom,carreau,nb_log), by = c("depcom","carreau"))

# Constitution des graphs :
data_graph <- data %>% 
  select(depcom,carreau,nb_log) %>% 
  filter(carreau != "FR_unallocated") %>% # on enlève les carreaux "solde"
  group_by(carreau) %>% 
  mutate(ncommune = sum(nb_log != 0)) %>% 
  ungroup()



# comptage du nombre de commune intersectée par carreau
table(data_graph$ncommune)
#0      1      2      3      4      5      6 
#14 312673 138574  19962   1360    120      6 



# On peut dégager les carreaux totalement inclus dans une commune de l'analyse égalemencarreau_full_included <- data %>%
fully_included <-data_graph %>% 
  filter(ncommune <= 1) 


# Je travaille avec la popuplation pour le moment comme nb log continet les résidences principales pour le moment
data_graph <- data_graph %>% 
  anti_join(fully_included, by = "carreau") %>% 
  select(depcom,carreau,nb_log,ncommune)

communes_isolees <- fully_included %>%
  filter(!depcom %in% data_graph$depcom) %>% 
  select(depcom) %>% 
  unique() 
# Ok ce sont les communes qui ont des carreaux inclus en leur sein 
# et aucun carreaux intersectant (au sens du nb log inside)

inner_join(fully_included, by = "carreau") %>% 
  select(depcom) %>%
  unique()


# contrôle de disparition de communes
setdiff(unique(data$depcom),unique(data_graph$depcom))


# existe-t'il des communes qui n'apparaissent que dans les fully included ?
# ie des composantes connexes à 1 élément

commune_carreaux_fu
nombre_de_commune_carreau_fully_include <- 
  sum(!unique(fully_included$depcom) %in% unique(data_graph$depcom)) # aucune !!

# je vais dégager les 0 count car ils ne contribuent pas non plus aux liens 
data_graph_2<-data_graph %>% 
  group_by(carreau) %>% 
  summarise(communes_reliee = paste0(depcom,collapse ="_"), 
            ncommune = unique(ncommune),nb_log= paste0(nb_log,collapse ="_"))


carreaux_to_communes_liees<- pblapply(seq_along(data_graph_2$communes_reliee),function(i){
  # i = 952
  # i =1
  # i = 46899
  
  ncommune <- data_graph_2$ncommune[i]
  commune_relie <- data_graph_2$communes_reliee[i]
  nb_log <- data_graph_2$nb_log[i]
  carreau <-data_graph_2$carreau[i]
  
  vec_com <- unlist(str_split(commune_relie,"_"))
  vec_nb_log <- as.numeric(unlist(str_split(nb_log,"_")))
  
  df <- data.frame(depcom = vec_com,nb_log =vec_nb_log)
  
  out <- t(combn(vec_com,2)) %>% 
    `colnames<-`(c("from","to")) %>%
    as.data.frame()
  
  out <- out %>% 
    filter(from != to) %>% 
    inner_join(df, by = c("from" = "depcom")) %>% 
    rename( inter1 = nb_log) %>% 
    inner_join(df, by = c("to" = "depcom")) %>% 
    rename( inter2 = nb_log)
  if(dim(out)[1] !=0)  out$carreau <- carreau
  return(out)
})


# y a des doublons
# TO DO faire un graph non  orienté avec attribut
# combn pour 2 parmi n
# CARTO Faire une carte
# inter1 donne bien le nombre dans la commune de départ intersecté avec le carreau

#carreaux_to_communes_liees <- readRDS("carreaux_to_communes_liees.RDS")
# Construction du graph
table_for_graph <-do.call(rbind,carreaux_to_communes_liees)

saveRDS(table_for_graph,"table_for_graph.RDS")
table_for_graph <- readRDS("table_for_graph.RDS")

#y a t'il des doublons dans les couples de commune ? auquel cas j'e sommerai le contenu des carreaux correspondant

table_for_graph %>% 
  group_by(from,to) %>% 
  summarise(n_carreau = n()) %>% 
  group_by(n_carreau) %>% 
  count()


## Ok donc il y a des communes qui sont reliées par beacuoupe de carreaux on réagrege par couple(from_to)
## pour avoir les comptes exacts

table_for_graph <- table_for_graph %>% 
  group_by(from,to) %>%
  summarise(inter1 = sum(inter1),inter2 = sum(inter2), carreau = paste0(carreau,collapse = "*"))


graph <- graph_from_data_frame(table_for_graph,directed = TRUE)
graph_tidy <- as_tbl_graph(graph)


# je regarde les clusters
graph_tidy <-
  graph_tidy %>%
  activate(nodes) %>%
  mutate(
    id_composante = (graph_tidy %>% clusters())$membership
  )

# J'extrait les communes connectées entre elles (je dégage les carreaux)

nodes <- as.data.frame(graph_tidy %>%  activate(nodes))
compo_to_nb_com <-nodes %>% 
  group_by(id_composante) %>% 
  summarise(n_commune =n()) %>% 
  arrange(-n_commune)


## histogramme
compo_to_nb_com %>% 
  filter(id_composante !=1) %>% 
  ggplot() +
  geom_bar(aes(x = n_commune))

compo_to_nb_com %>% 
  group_by(n_commune) %>% 
  count()

max(nodes$id_composante) # 1895 composantes connexes ?

## filtre des liens potentiellement problématiques
edge <- graph_tidy %>% 
  activate(edges) %>% 
  as.data.frame()

# Faire des unions de polygone de communes et les représenter
# stats, nombre de clusters, nombre de carreaux sous le seuil  
# Représentation graphique des clusters en union
# cluster diff_interne externe des chiffres
# taille en nombre de communes



### Carte 1 : Exemples d"taillés de composantes connexes

num_compo <- 1251
# num_compo <- 9
nodes_compo <- nodes %>%
  mutate(id_node = row_number()) %>% 
  filter(id_composante == num_compo) 

carreaux_cluster <- edge %>% 
  filter(from %in% nodes_compo$id_node | to %in% nodes_compo$id_node) %>%
  select(carreau) %>%
  pull()

liste_carreau<- data.frame(carreau = str_split(carreaux_cluster,"[*]")  %>% unlist())
polygone_carreau <- carreaux_to_polygon(liste_carreau, var_carreau = "carreau") %>% st_transform(crs = 4326)

liste_carreau2<- data %>% ungroup() %>% filter(depcom %in% nodes_compo$name ) %>% filter( carreau != "FR_unallocated") %>%  select(carreau)
polygone_carreau_hors_inter <- carreaux_to_polygon(liste_carreau2, var_carreau = "carreau") %>% st_transform(crs = 4326)

polygone_commune<- com_2021 %>% 
  filter(code %in% nodes_compo$name)


# plot(polygone_commune$geometry)
# plot(polygone_carreau$geometry, col = "red", add = TRUE)

# recupération des intersections
info_croisement <- data %>% 
  select(depcom,carreau,nb_log) %>% 
  filter(depcom %in% nodes_compo$name & carreau %in% liste_carreau$carreau)

inter_carreau_commune <-st_intersection(
  polygone_commune %>% select(code),
  polygone_carreau %>% select(carreau)
)

inter_carreau_commune <- info_croisement %>% 
  left_join(inter_carreau_commune, by = c("depcom"="code","carreau"))


#### Carte
dep <- st_read("S:/PROD/creacartes/pd/fichiers-ihm/2021/francemetro/dep_francemetro_2021.gpkg") %>% 
  st_transform(crs = 4326)

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
    data = dep %>% filter(code %in% polygone_commune$dep ),
    color = "black",
    weight = 3,
    fillOpacity = 0,
    opacity = 1,
    group = "departement"
  ) %>% 
  addPolygons(
    data = polygone_commune,
    color = "#3FC8FC",
    weight = 2,
    fillOpacity = 0,
    opacity = 1
  ) %>% 
  addPolygons(
    data = polygone_carreau_hors_inter,
    color = "orange",
    weight = 2,
    fillOpacity = 0,
    group = "ensemble des carreaux",
    opacity = 1
  ) %>% 
  addPolygons(
    data = polygone_carreau,
    color = "red",
    weight = 2,
    fillOpacity = 0,
    opacity = 1,
    group ="carreaux à l'intersection des communes"
  ) %>% 
  addPolygons(
    data = inter_carreau_commune$geometry,
    color = ifelse(inter_carreau_commune$nb_log < 11,"red","#6E3DFF"),
    weight = 2,
    fillOpacity = 0.5,
    group = "intersections",
    highlightOptions = highlightOptions_defaut,
    label  =  with(inter_carreau_commune, 
                   sprintf(
                     "<b> id commune : </b> %s  <br/> <b> id carreau : </b>  %s <br/> <b> Nombre de logements au RIL median : </b>  %s", ### c'est une d�finition de format qui vient du C
                     depcom, carreau, round(nb_log,1)
                   ) %>% lapply(htmltools::HTML)
    )
  ) %>% 
  addScaleBar(position="bottomright") %>% 
  hideGroup(c("departement","ensemble des carreaux","carreaux à l'intersection des communes","intersections")) %>% 
  addLayersControl(
    overlayGroups = c("departement","ensemble des carreaux","carreaux à l'intersection des communes","intersections"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  addScaleBar(position="bottomright")


carte

htmlwidgets::saveWidget(carte, file=paste0("carte_compo_",num_compo,".html"),selfcontained = TRUE)



## carte 2 :  cartes des composantes connexes en tant qu'union de commune 

# juste à faire un group by
nodes_geom <-nodes %>%
  left_join(com_2021 %>% select(code,geometry),by = c("name"="code")) 

l <- split(nodes_geom,nodes_geom$id_composante)
sf_use_s2(FALSE)

# 


# la première composante fait tout crasher, on utilise reduce
n_com_comp1 <- dim(l[[1]])[1]
compo1 <- Reduce(st_union,l[[1]]$geometry)
saveRDS(compo1,"union_compo_1.rds")
compo1 <- readRDS("union_compo_1.rds")

autres_compo <- nodes_geom %>% 
  filter(id_composante !=1) %>% 
  group_by(id_composante) %>% 
  summarise(geometry_union = st_union(geometry),
            n_communes = length(id_composante))


# 
# noeuds_uniques <- data %>% 
#   filter(carreau!= "FR_unallocated") %>% 
#   group_by(depcom) %>% 
#   summarise(ncarreau = n()) %>% 
#   filter(ncarreau==1) %>% 
#   left_join(com_2021 %>% select(code),by = c("depcom"="code"))
# 


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
    data = dep,
    color = "black",
    weight = 3,
    fillOpacity = 0,
    opacity = 1,
    group = "departement"
  ) %>% 
  addPolygons(
    data = autres_compo$geometry_union,
    color = "#FFF821",
    fillColor = "red",
    weight = 1,
    fillOpacity = 0.5,
    opacity = 1,
    group  ="autres composantes",
    highlightOptions = highlightOptions_defaut,
    label = sprintf("<b> Nombre de communes : </b> %s",autres_compo$n_communes)%>% lapply(htmltools::HTML)
  ) %>% 
  addPolygons(
    data = compo1,
    color = "#47F5D3",
    weight = 2,
    fillOpacity = 0.3,
    opacity = 1,
    group = "enorme composante",
    label = sprintf("<b> Nombre de communes : </b> %s",n_com_comp1)%>% lapply(htmltools::HTML)
  )%>% 
  addScaleBar(position="bottomright") %>% 
  hideGroup(c("enorme composante","autres composantes")) %>% 
  addLayersControl(
    overlayGroups = c("departement","enorme composante","autres composantes"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  addScaleBar(position="bottomright")


carte


htmlwidgets::saveWidget(carte, file="carte_compo_gobales.html",selfcontained = TRUE)


