
##Loading Data


library(bibliometrix)
library(igraph)
library(dplyr)

#Reading the data
#This data was from WOS using crawling 

D<-readFiles("C://Users/admin/Desktop/final.ahrq.txt")
M <- convert2df(D, dbsource="isi",format="plaintext")


cr<-cocMatrix(M,Field = "CR",sep =";")
dim(cr)

#The function cocMatrix makes the data bipartite matrix(96*1939)
#The Field "CR" means Cited references
#There are 96 articles in AHRQ
#There are 1939 references in 96 articles


z<-biblioNetwork(M,analysis ="co-citation",network = "references")
g1<-graph.adjacency(z,mode="undirected",weighted = T)
#The function biblioNetwork makes the data co-citation matrix(1939*1939)






## Data cleasing


A<-as.matrix(cr)
name<-dimnames(A)[[2]]
name[1:10]


for(i in 1:length(name)){
  
  tmp<-strsplit(name[i],' ')[[1]]
  tmp[1]<-tolower(tmp[1])
  ttmp<-strsplit(tmp[1],'')[[1]]
  ttmp[1]<-toupper(ttmp[1])
  tmp[1]<-paste(ttmp,collapse = "")
  
  ind<-grep("\\d",tmp)
  
  tmp[ind]<-paste0(tmp[ind],",")
  tmp[ind-1]<-paste0(tmp[ind-1],",")
  name[i]<-paste(tmp,collapse=" ")
}

dimnames(A)[[2]]<-name

h<-data.frame(colSums(A))
h$name<-name
#write.csv(arrange(h,desc(h$colSums.A.)),"ahrq.totla.freq.author.csv")

name[1:10]



## Centrality and Igraph


co.rf<-t(A) %*% A
co.rf.g<-graph.adjacency(co.rf,mode = "undirected",weighted = T)
#t(A)*A is co-citation matrix


deg.rf<-degree(co.rf.g)
deg.rf2<-as.data.frame(deg.rf)
deg.rf2$name<-name
deg.rf3<-arrange(deg.rf2,desc(deg.rf2$deg))
deg.rf3[1:31,]

#Top30 nodes(references) by degree from the total data

e1<-eigen_centrality(co.rf.g)[1]
e11<-data.frame(e1)

b1<-betweenness(co.rf.g)
b11<-data.frame(b1)


d1<-degree(co.rf.g)
d11<-data.frame(d1)


dc1<-d11/(dim(d11)[1]-1) 

cen<-cbind.data.frame(e11,b11,d11)
cen$name<-name
names(cen)<-c("eigen_cen","betweenness","degree_cen","name")

head(cen)
#You can choose the Top30 nodes by eigen centrality or beteweenness or degree centrality

E(co.rf.g)$width = E(co.rf.g)$weight
V(co.rf.g)$size<-deg.rf

#Edge size is weighted edge
#Vertex size is degree

induce<-induced.subgraph(co.rf.g,V(co.rf.g)$size >=207)

#Inducing top 30 nodes

V(induce)$initial<-ifelse(V(induce)$size>=240,V(induce)$name,"")
g2<-simplify(induce,remove.loops=TRUE,remove.multiple = T)

V(g2)$color<-ifelse(V(g2)$size>=593,"black",
                    ifelse(V(g2)$size>=403,"Dim gray","Gray"))

V(g2)$initial2<-" "
set.seed(12)
plot(g2,layout=layout.fruchterman.reingold, vertex.color=V(g2)$color, vertex.size=V(g2)$size/70, vertex.label.color="black",
     vertex.label.cex=0.8, vertex.label.dist=1, edge.curved=0,edge.width=E(g2)$weight/3,edge.color="gray",vertex.label=V(g2)$initial)



##Community detection

eig<-leading.eigenvector.community(co.rf.g)
lou<-cluster_louvain(co.rf.g)

V(co.rf.g)$lou<-lou$membership
V(co.rf.g)$eig<-eig$membership

#You can check the cluster including 
