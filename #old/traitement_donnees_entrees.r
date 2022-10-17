
install.packages("btb")
install.packages("remotes")
install.packages("tidygraph")

library(data.table)
library(sf)
library(igraph)
library(Matrix)
library(tidygraph)
library(leaflet)
# installation diffman pour récupérer les donénes jouets
remotes::install_github("https://github.com/InseeFrLab/diffman-1")

# applishare
system("mc cp s3/cguillo/data_rp.rds data_rp.rds") # données avec la pop seulement

source("fonctions_utiles.R")

communes <- st_read("../commune.shp")
donnees_rp <- readRDS("data_rp.rds")

df <- setDT(donnees_rp)
df <- df[substr(z1,1,2) == "01"] 
# df <- setDT(diffman::t_ex)[,.(nb_obs=length(id)),by=.(z1,z2)]


# création de la table des liens

## graph en sortie non orienté => pas de doublon A-> B et B-> A
## on a pas agrégé par from,to ici encore 
link_table <- build_link_table(df) 
link_table %>% head()

## etat des lieux des composantes connexes
ltable <- long_table(link_table)
ltable[,.(n_z1 = length(unique(z1))),.(id_comp)][rev(order(n_z1))]


# représentation d'une composante connexe
situation_table <- link_table[id_comp %in% c(1,15)] ;unique(c(situation_table$from,situation_table$to))

## construction des géométries nécessaires à l'affichage (en dehirs du package) 
  liste_carreau <- situation_table$z2
  
  polygone_carreau <-
    carreaux_to_polygon(data.frame(carreau = liste_carreau), var_carreau = "carreau") %>%
    st_transform(crs = 4326)
  
  geom_z1 <- communes %>%
    filter(code %in% unique(c(situation_table$from,situation_table$to))) %>%
    select(code) %>%
    rename(z1 = code)
  
  geom_z2 <- polygone_carreau %>%
    rename(z2 = carreau) %>%
    select(-x,-y)

draw_situation(situation_table,geom_z1,geom_z2,threshold = 20)




# Controler les sorties de diffman sur des petits examples
# Ok cette ligne d'instruction ce sera pour après dans les recherches de pb etc..

# search pb_diff
# en rentrée je mets forcément une table de la forme link table
situation_table <- link_table[id_comp %in% c(15)] ;unique(c(situation_table$from,situation_table$to))
df <- situation_table[ ,
                          .(
                            nb_obs_from = sum(nb_obs_from),
                            nb_obs_to = sum(nb_obs_to),
                            id_comp = unique(id_comp),
                            z2 = paste0(z2,collapse = "-"),
                            n_z2 = unique(n_z2)
                          ),
                          by = .(from,to)
]




# TO DO :
# faire un rmd,
# maitriser diffman sur une petite comoposante




# Rq :

# Ok les carreaux liant plus que 3 communes ont des oobsvertions comptées en doublons :
# cf triangle de com A1 A2 A3 lavec un carreau liant ca
# ca inter A1 est compter dans la liaison A1-A2 et A1-A3

