#arbol con ggtree
library(ggtree)
library(ape)
library(pwalign)
library("msa")
library("BSgenome")
library("ggmsa")
library(phangorn)
library(phytools)
library(seqinr)

arbol_1 <- phangorn::upgma(matriz_dist)
plot(arbol_1)


arbol <- ggtree(arbol_1, layout="slanted")+
  geom_treescale(linesize = 2, )+ #grosor de la linea/editar la linea
  geom_tiplab(size = 3, color = "purple")+
  geom_highlight(node = 1:2, color = "orange", fill="orange")+
  geom_point2(aes(subset=node==1), color='darkorange', size=3)+
  geom_point2(aes(subset=node==2), color='darkorange', size=3)+
  geom_highlight(node = 4:5, color = "green", fill="green")+
  geom_highlight(node = 3, color = "pink", fill="pink")+
  geom_point2(aes(subset=node==3), color='red', size=3)+
  geom_point2(aes(subset=node==4), color='darkgreen', size=3)+
  geom_point2(aes(subset=node==5), color='darkgreen', size=3)

plot(arbol)

ggsave("RESULTADOS/Arbol modificado.png", plot = arbol)
