library("msa")
library("BSgenome")
library("ggmsa")
browseVignettes("msa")
library(ape)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)
install.packages("seqinr")
library(seqinr)
install.packages("ggtree")
library(ggtree)
browseVignettes("ggtree")


#Utilizar varias secuencias descargadas y hacer alineamientos de dos tipos y un arbol
hemoglobinas <- readAAStringSet("Datos/DivergentGlobins.fasta")
print(hemoglobinas)

alineamiento_hemo <- msa (hemoglobinas, "Muscle")
print(alineamiento_hemo, show = "complete")

alineamiento_hemo2 <- msa(hemoglobinas, "ClustalW")
print(alineamiento_hemo2, show = "complete")

ggmsa(hemoglobinas)

#ggmsa graficas del alineamiento.


##arbol filogenetico con el metodo UPGMA
#primer calculo de los alineamientos.
matriz_dist <- ape::as.AAbin(alineamiento_hemo)
matriz_dist <- ape::dist.aa(matriz_dist)
matriz_dist

matriz_dist_2 <- ape::as.AAbin(alineamiento_hemo2)
matriz_dist_2 <- ape::dist.aa(matriz_dist_2)
matriz_dist_2

pdf("RESULTADOS/Arbol alineamiento 1. Metodo UPGMA.pdf")
arbol_1 <- phangorn::upgma(matriz_dist)
plot(arbol_1)
dev.off()

arbol_1

pdf("RESULTADOS/Arbol alineamiento 2 Metodo UPGMA.pdf")
arbol_2 <- phangorn::upgma(matriz_dist_2)
plot(arbol_2)
dev.off()

#######################################
#Intento de arbol utilizando el msa
class(alineamiento_hemo)
distancia_1 <- dist.alignment(alineamiento_hemo, "identity")
as.matrix(distancia_1)

hemoTree <- nj(distancia_1)
plot(hemoTree)

##### Arbol con ggtree: solo del primer alineamiento.
pdf("RESULTADOS/Arbol alineamiento 1 Metodo ggtree.pdf")
p1 <- ggtree(arbol_1)
plot(p1)
dev.off()
