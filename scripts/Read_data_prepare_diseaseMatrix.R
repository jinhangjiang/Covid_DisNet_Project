#### Read the AZDHS Hospital discharge data file and create disease co-occurrence network 
## Phase 2: Use node2vec for embedding and visualize using T-SNE

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
  rare_diseases <- which(count_diag < 10 )  ## Around 24196 out of 31525 are rare diagnoses in 2019-2
  
  ## Also remove codes belonging to categories: o-pregnancy, s,t-injury or poisoning, v-z - (External):
  
  names_diag0 <- dtm$dimnames[[2]]
  
  ### At this stage, you can remove unused objects to make space in your memory, to run other memory intensive code
  ##sort( sapply(ls(),function(x){object.size(get(x))}))
  ### rm(processedCorpus);rm(aritemp);rm(ari12);rm(ari14);rm(ari13);
  ### mem_used()
  ### memory.size()
  
  ### Exclusion index: 
  
  excl_index = c(rare_diseases,startsWith(names_diag0, c('o','s','t','v','w','x','y','z')))
  
  d2 <- d1[,-excl_index]  
  names_diag <- dtm$dimnames[[2]][-excl_index]
  
  count_diagnoses <- colSums(d2)
  
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
            ## Significant links only using t-test 
            if(  (trm_CC[i,j] <= 0 ) || (2*pt(-trm_CC[i,j]/sqrt((1-trm_CC[i,j]^2)/(N_test -2)),df= (N_test-1)) >= 0.05) ) 
            {
              trm_CC[i,j] <- NA  
            }   
            
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
  
  giant.component <- function(graph) { 
    cl <- clusters(graph) 
    induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))} 
  
  clean_edges <- function(graph)
  {
    graph <- subgraph.edges(graph, E(graph)[!is.na(E(graph)$weight)], delete.vertices = TRUE)
    return(graph)
  }
  
  CC_g <- clean_edges(CC_graph) ## 
  simple_g <- clean_edges(simple_graph) ## 

return(CC_g)
}

disnet20192 <- generate_network(tr20192)
disnet20201 <- generate_network(tr20201)
disnet20202 <- generate_network(tr20202)

e_disnet20192 <- get.data.frame(disnet20192)
e_disnet20201 <- get.data.frame(disnet20201)
e_disnet20202 <- get.data.frame(disnet20202)

names(e_disnet20192) <- c("Source","Target","Weight")
names(e_disnet20201) <- c("Source","Target","Weight")
names(e_disnet20202) <- c("Source","Target","Weight")

write.csv(e_disnet20192, "data/Edgelist2019_2.csv",row.names=F)
write.csv(e_disnet20201, "data/Edgelist2020_1.csv",row.names=F)
write.csv(e_disnet20202, "data/Edgelist2020_2.csv",row.names=F)

dis_categories <- read.csv("data/icd10cm_codes_2021.csv")

d20192 <- disnet20192
d20201 <- disnet20201
d20202 <- disnet20202

temp = data.frame(Id = V(d20192)$name)
temp2 <- left_join(temp, dis_categories, by =c("Id"="Id"))
V(d20192)$Description = temp2$Description
V(d20192)$Category = temp2$Category


## Visualizing using igraph networkD3: https://christophergandrud.github.io/networkD3/
library(networkD3)

# Convert to object suitable for networkD3
d20192_d3 <- igraph_to_networkD3(d20192, group = V(d20192)$Category)

# Create force directed network plot
# forceNetwork(Links = d20192_d3$links, Nodes = d20192_d3$nodes, 
#              Source = 'source', Target = 'target', 
#              NodeID = 'name', Group = 'group')  %>% saveNetwork(file = "C:/Users/k639s258/OneDrive - The University of Kansas/RootFolder/Research/Covid_disNet/htmlout/Disnet20192.html",selfcontained = TRUE)

visNetwork(nodes, edges, height = "500px") %>%
  visIgraphLayout() %>%
  visNodes(size = 10)

simpleNetwork(d20192_d3)  %>% saveNetwork(file = "/htmlout/Disnet20192.html",selfcontained = TRUE)


## The ICD10 codes for COVID-19  ## Use U071
## https://www.who.int/classifications/icd/COVID-19-coding-icd10.pdf?ua=1
#### U07.1 COVID-19, virus identified
#### U07.2 COVID-19, virus not identified
  #o Clinically-epidemiologically diagnosed COVID-19
  #o Probable COVID-19
  #o Suspected COVID-19

## node2vec embedding in R: Does not work
library(node2vec)
embed_disnet20192 <- node2vecR(data = e_disnet20192[,c(1:2)], num_walks = 5, walk_length = 5, dim = 10)
embed_disnet20201 <- node2vecR(e_disnet20201, num_walks = 5, directed = FALSE, walk_length = 5, dim = 10)
embed_disnet20202 <- node2vecR(e_disnet20202, num_walks = 5, directed = FALSE, walk_length = 5, dim = 10)


