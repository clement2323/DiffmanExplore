
install.packages("dplyr")
install.packages("datatable")
rm(list = ls())

library(dplyr)
library(Matrix) 
library(data.table)

system("mc cp s3/cguillo/pop_carreau_dmrg.rds pop_carreau_dmrg.rds") # invite de commande
install.packages("remotes")
remotes::install_github("https://github.com/InseeFrLab/diffman-1")

# Lancement naïf : ça marche !
# diffman::t_ex

res_diff <- diffman::find_pbm_diff(
  t_ind = diffman::t_ex,
  threshold = 8,
  max_agregate_size = 15
  )

# Ok ::: pour chercher les fonctions non visibles
# La fonction find_id_obs_risque de diffman en, fait partie , ici je la réutilise à partir des résultats présédents
# uil faut construire la liste de vecteur de zonage et il sort les informations à risque associé
# finalement nous on ne veut pas ces individus, mais plutôt l'intersection carreau x commune qui dépasse
# TO DO faire un test de diffman sur une composante connexe isolée, celle de l'exemple !
list_agregat <- strsplit(res_diff$agregat_z1,"-")

test <- diffman:::find_id_obs_risque( 
  list_agregat = list_agregat ,
  t_ind = diffman::t_ex,
  threshold = 11,
  verbose = TRUE
)


# " agregate etc.."
source("../R/fct_merging.R") 

# " test fus etc.."
source("../R/fct_merging_test.R") 

# " m_crois"
source("../R/fct_reshaping.R")

# " Decompose mcrois"
source("../R/fct_splitting.R")

# search_diff_agregat
source("../R/fct_search_diff.R")

# differencierRcpp
source("../R/RcppExports.R")

find_pbm_diff_tab_clem <- function(t_crois_init, threshold, max_agregate_size, save_file = NULL, simplify = TRUE, verbose = TRUE){
  # threshold = 7; max_agregate_size = 15;save_file = NULL; simplify = TRUE; verbose = TRUE
  # t_crois_init = diffman::t_ex %>% group_by(z1,z2) %>% summarize(nb_obs = n()) %>% ungroup()
  # define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = NULL
  
  if(verbose) message("******** Start of the process *********")
  
  to_save <- !is.null(save_file) #si save_file est renseignee, cela permet de sauvegarder les matrices croisees apres agregations
  
  t_crois_init <- as.data.table(t_crois_init) #conversion to data.table format
  
  if(verbose) message("< --- Creation of the crossing matrix --- >")
  # Ok les lignes qui utilisaient simplify_z2_rm et tabl_crois 
  # servait  à :
  # i) filtrer les zonages z2 qui n'étaient pas à cheval sur au moins 2 zonages z1
  # ii) constituer la tabulation Z1 x Z2 à partir de la table individuelle sans les zonages Z2 évoqués pprécedemment
  # on à déjà tut ce qu'il faut si ce n'est qu'il faut aussi dégager les Z2 -> à traiter quand mêm après -> diffman le fait-il ?
  
  t_z2 <- t_crois_init[,.(nb_z1=length(z1)),by=.(z2)]
  t_crois <- t_crois_init[z2 %in% t_z2[nb_z1>=2]$z2]
  t_crois_final <- t_crois[,z2_b:=paste0(sort(z2),collapse="-"),by=.(z1)][ #z2_b est une "étiquette" qui donne les zones de z1 recouvertes par chaque zone de z2
    ,.(nb_obs=sum(nb_obs)),by=.(z1,z2_b)]
  
  colnames(t_crois_final)[colnames(t_crois_final)=="z2_b"] <- "z2"
  
  # t_crois <- t_crois[z1 != "blanchi" & z2 != "blanchi"]
  m_crois <- matrix_crois(t_crois_final)
  
  #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_",save_file,".RDS"))
  
  if(simplify){ #one can choose to skip these steps of graph reduction if desired
    
    if(verbose) message("< --- Merging method 1 --- > ")
    m_crois <- agregate(m_crois, threshold, methode = "m1", verbose = verbose)
    
    #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m1_",save_file,".RDS"))
    
    if(verbose) message("< --- Merging methods 1 and 2 --- >")
    m_crois <- agregate(m_crois, threshold, methode = "both", verbose = verbose)
    
    # warnings()
    #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m2_",save_file,".RDS"))
    
    if(sum(dim(m_crois)==0)>0) {
      message("No differentiation problems detected !")
      return(NULL)
    }
    
    if(verbose) message("< --- Splitting the graph --- >")
    l_decomp <- decompose_m_crois(m_crois, max_agregate_size)
  }
  else{
    l_decomp <- comp_connexe_list(m_crois)
  }
  
  if(verbose) message("< --- Exhaustive search of differentiation problems --- >")
  # sauvegarder l_decomp pourrait être intéressant
  l_ag <- search_diff_agregate(l_decomp, threshold, max_agregate_size) 
  # dans fct_search_diff test composante connexe par composante connexe
  l_ag <- desagregate_list(l_ag)
  
  # Ok dans l_ag un élément est un vecteur d'élément de z1, 
  # l'union de ces éléments désigne l'agrégat d'elt de z1 différenciable
  # but du jeu : identifier les intersections concernées par le risque (diff interne ou externe) 
  
  return(l_ag)  
}




retrouver_zones <- function(l_ag){
out <- lapply(l_ag,function(agregat_z1){
  # agregat_z1 <- l_ag[[9]]
  
  #je récupère tous les croisements qui font intervenir l'agregat dans la table avec tous lkes croisements
  t_crois_ag <- t_crois_init[z1 %in% agregat_z1] 
  
  #j'extirpe les zonages z2 concernés
  z2_ag <- as.character(unique(t_crois_ag$z2))
  
  # je regarde où interviennet également les z2 en plus de l'agrégat initial
  zone_externe <- t_crois_init[z2 %in% z2_ag] # a comprendre comme ensemble des carreaux partiellementou totalment inclu
  
  # les z2 carreaux non totalement inclus sont ceux qui apparaissent dans une zone de z1 hiors agregat
  z2_not_fully_include <-as.character(unique(zone_externe[!z1 %in% agregat_z1]$z2))
  zone_interne<- t_crois_init[z1 %in% agregat_z1 & !z2 %in% z2_not_fully_include]# a comprendre comme ensemble des carreaux z2 totalment inclus dans l'agregat z1
  
  diff_externe <- sum(zone_externe$nb_obs) - sum(t_crois_ag$nb_obs)
  diff_interne <- sum(t_crois_ag$nb_obs) - sum(zone_interne$nb_obs)  
  
  liste_res = list(
    agregat_z1 = paste0(agregat_z1,collapse="-"),
    agregat_z2_inter_z1 = z2_ag,
    zone_externe = NULL,
    zone_interne = NULL
  )
  if(diff_externe < threshold){
    liste_res$zone_externe <- zone_externe
  }
  if(diff_interne < threshold){
    liste_res$zone_interne <- zone_interne
  }
  return(liste_res)
})
out
}


# la fonction retrouver zone est en stand by pour le moment
#" Application cas réel"
# Lancement avec la table pop 
# ça ne peut bien marcher car on est déjà tabulé alors que diffman attend de l'indiv
# test sur un département
# Décorticage de la fonction find_pbm_diff 
# 

donnees_rp <- readRDS("pop_carreau_dmrg.rds")
donnees_rp <- donnees_rp %>% select(depcom,carreau,nb_log) %>% `colnames<-`(c("z1","z2","nb_obs"))

# mc cp s3/cguillo/data_rp.rds data_rp.rds
donnees_rp <- readRDS("data_rp.rds")

Sys.time()
res <- find_pbm_diff_tab_clem(donnees_rp,threshold = 11, max_agregate_size =15)
Sys.time()


saveRDS(donnees_rp,"data_rp.rds")
readRDS("data_rp.rds")