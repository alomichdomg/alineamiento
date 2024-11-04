library("msa")
library("BSgenome")
library("ggmsa")
browseVignettes("ggmsa")
library(ape)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)

#Utilizar varias secuencias descargadas y hacer alineamientos de dos tipos y un arbol
hemoglobinas <- readAAStringSet("Datos/DivergentGlobins.fasta")
print(hemoglobinas)

alineamiento_hemo <- msa (hemoglobinas, "Muscle")
print(alineamiento_hemo, show = "complete")

alineamiento_hemo2 <- msa(hemoglobinas, "ClustalW")
print(alineamiento_hemo2, show = "complete")

ggmsa(hemoglobinas)

#ggmsa graficas del alineamiento.


##arbol filogenetico:
#primer calculo de los alineamientos.
matriz_dist <- ape::as.AAbin(alineamiento_hemo)
matriz_dist <- ape::dist.aa(matriz_dist)
matriz_dist

matriz_dist_2 <- ape::as.AAbin(alineamiento_hemo2)
matriz_dist_2 <- ape::dist.aa(matriz_dist_2)
matriz_dist_2

arbol_1 <- phangorn::upgma(matriz_dist)
plot(arbol_1)

arbol_2 <- phangorn::upgma(matriz_dist_2)
plot(arbol_2)
