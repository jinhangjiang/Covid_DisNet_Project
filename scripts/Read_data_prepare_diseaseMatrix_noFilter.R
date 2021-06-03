#### Read the AZDHS Hospital discharge data file and create disease co-occurrence network WITHOUT FILTER


library(tm)
library(topicmodels)
library(tidyr)
library(reshape)
library(igraph)
library(Hmisc)
library(plyr)
library(dplyr)
library(Rcpp)
library(parallel)
library(Matrix)
library(data.table)
library(caret)
library(readr)


setwd("C:\\Users\\k639s258\\OneDrive - The University of Kansas\\RootFolder\\Research\\Covid_disNet")

colwidth <- read.csv("data/colwidth.csv")

tr20192 <- read_fwf("raw_data\\puf_201902ip.txt", fwf_widths(colwidth$Length, colwidth$Varname) )
tr20201 <- read_fwf("raw_data\\puf_202001ip.txt", fwf_widths(colwidth$Length, colwidth$Varname) )
tr20202 <- read_fwf("raw_data\\puf_202002ip.txt", fwf_widths(colwidth$Length, colwidth$Varname) )

tr20192$Dur <- "2019Y2"
tr20201$Dur <- "2020Y1"
tr20202$Dur <- "2020Y2"

generate_network <- function(data_in)
{

  ## Take the Diag columns and create the DTM
  diag <- data_in[,c(2,seq(7,54,2))] 
  
  #diag = data.frame(diag)
  
  ## Single vector
  diag_vec <- unite(diag, col= "all",match(names(diag),names(diag)),sep=" ",remove=T, na.rm=T)
  
  diag_vec <- data.frame(diag_vec)
  diag_a <- as.vector(diag_vec[,1])
  
  system.time( corpus <- Corpus(VectorSource(diag_a)) )  ### This takes 10 minutes to run, for ~380k records 
  
  #processedCorpus <- tm_map(corpus, PlainTextDocument)
  #dtm <- DocumentTermMatrix(processedCorpus)
  dtm <- DocumentTermMatrix(corpus)
  
  Dense <- sparseMatrix(dtm$i,dtm$j,x=dtm$v)
  
  d1 <- Dense
  count_diag <- colSums(d1) # Prevalence of each diagnosis
  
  # rare_diseases <- which(count_diag < 10 )  ## Around 24196 out of 31525 are rare diagnoses in 2019-2
  
  ## Also remove codes belonging to categories: o-pregnancy, s,t-injury or poisoning, v-z - (External):
  
  names_diag0 <- dtm$dimnames[[2]]
  
  ### At this stage, you can remove unused objects to make space in your memory, to run other memory intensive code
  ##sort( sapply(ls(),function(x){object.size(get(x))}))
  ### rm(processedCorpus);rm(aritemp);rm(ari12);rm(ari14);rm(ari13);
  ### mem_used()
  ### memory.size()
  
  ### Exclusion index: 
  
  # excl_index = c(rare_diseases,startsWith(names_diag0, c('o','s','t','v','w','x','y','z')))
  excl_index = c(startsWith(names_diag0, c('o','s','t','v','w','x','y','z')))
  
  d2 <- d1[,-excl_index]  
  names_diag <- dtm$dimnames[[2]][-excl_index]
  
  count_diagnoses <- colSums(d2)
  
  icd9_list <- data.frame(Id = names_diag,weight = count_diagnoses)
  
  #icd9_list <- data.frame(Id = names_diag,weight = count_diagnoses) 
  #write.csv(icd9_list,"data/ICD9_weights_Arizona20192.csv",row.names=F)
  
  system.time(tot <- crossprod(d2)) 
  N <- nrow(d2)
  trm <- as.matrix(tot)
  
  trm_CC <- trm
  P <- diag(trm)   ## Diagonal elements
  
  system.time(
    for(i in 1:nrow(trm))
    {
      # print(i)
      PX <- P[i]
      for (j in 1:ncol(trm))
      {
        if(trm[i,j]!=0)
        {
          if(i!=j) 
          {
            C <- trm[i,j]
            PY <- P[j]
            trm_CC[i,j] <- sqrt(2)*C/sqrt(PX^2 + PY^2)   
            N_test <-  max(PX,PY)
            # ## Significant links only using t-test 
            # if(  (trm_CC[i,j] <= 0 ) || (2*pt(-trm_CC[i,j]/sqrt((1-trm_CC[i,j]^2)/(N_test -2)),df= (N_test-1)) >= 0.05) ) 
            # {
            #   trm_CC[i,j] <- NA  
            # }   
            # 
          }
        }
        
      }    
    }
  )
  
  create_undirected_graph <- function(adj_mat,names_diag,P)
  {
    t_g <- graph.adjacency(adj_mat,weighted=TRUE,mode = c("undirected"))
    simple_graph <- simplify(t_g,remove.loops=TRUE)
    V(simple_graph)$name <- names_diag
    V(simple_graph)$weight <- P
    return(simple_graph)
  }
  
  create_directed_graph <- function(adj_mat,names_diag,P)
  {
    t_g <- graph.adjacency(adj_mat,weighted=TRUE,mode = c("directed"))
    simple_graph <- simplify(t_g,remove.loops=TRUE)
    #simple_graph <- t_g
    V(simple_graph)$name <- names_diag
    V(simple_graph)$weight <- P
    return(simple_graph)
  }
  
  simple_graph <- create_undirected_graph(trm,names_diag,P)
  CC_graph <- create_undirected_graph(round(trm_CC,4),names_diag,P)
  
  # giant.component <- function(graph) { 
  #   cl <- clusters(graph) 
  #   induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))} 
  # 
  clean_edges <- function(graph)
  {
    graph <- subgraph.edges(graph, E(graph)[!is.na(E(graph)$weight)], delete.vertices = TRUE)
    return(graph)
  }
  
  CC_g <- clean_edges(CC_graph) ## 
  simple_g <- clean_edges(simple_graph) ## 

return(list(icd9_list, CC_g))
}

disnet20192_l <- generate_network(tr20192)
disnet20201_l <- generate_network(tr20201)
disnet20202_l <- generate_network(tr20202)

disnet20192 <- disnet20192_l[[2]]
disnet20201 <- disnet20201_l[[2]]
disnet20202 <- disnet20202_l[[2]]

count20192 <- disnet20192_l[[1]]
names(count20192)[2] <- "Prevalence_20192"
count20201 <- disnet20201_l[[1]]
names(count20201)[2] <- "Prevalence_20201"
count20202 <- disnet20202_l[[1]]
names(count20202)[2] <- "Prevalence_20202"

dis_categories <- read.csv("data/icd10cm_codes_2021.csv")

dis_categories_prev <- left_join(dis_categories, count20192,by =c("Id"="Id"))
dis_categories_prev <- left_join(dis_categories_prev, count20201,by =c("Id"="Id"))
dis_categories_prev <- left_join(dis_categories_prev, count20202,by =c("Id"="Id"))

rowsums <- rowSums(dis_categories_prev[,c(5,6,7)])
dis_categories_prev2 <- dis_categories_prev[which(!is.na(rowsums) ),] 

write.csv(dis_categories_prev2, "data/icd10cm_codes_prev_2021_unfiltered.csv", na="",row.names=F)

e_disnet20192 <- get.data.frame(disnet20192)
e_disnet20201 <- get.data.frame(disnet20201)
e_disnet20202 <- get.data.frame(disnet20202)

names(e_disnet20192) <- c("Source","Target","Weight")
names(e_disnet20201) <- c("Source","Target","Weight")
names(e_disnet20202) <- c("Source","Target","Weight")

write.csv(e_disnet20192, "data/Edgelist2019_2_unfiltered.csv",row.names=F)
write.csv(e_disnet20201, "data/Edgelist2020_1_unfiltered.csv",row.names=F)
write.csv(e_disnet20202, "data/Edgelist2020_2_unfiltered.csv",row.names=F)

save.image("Generate_unfiltered_disNets.RData")
