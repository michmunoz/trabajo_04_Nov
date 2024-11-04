

### LEER EL ARCHIVO FASTA


library(BiocManager)
library (Biostrings)
library(msa)
library (ape)
library (seqinr)

sec <- readDNAStringSet("data/secuenciasCAA.txt")
sec


# hacer los alineamientos

alin1 <- msa (sec, method = "ClustalW")
alin1

alin2 <- msa (sec, method = "Muscle")
alin2

## se realizan las matrices de distancias

a1 <- msaConvert (alin1, type="seqinr::alignment")
a1

dist1 <- dist.alignment(a1, "identity")
dist1

a2 <- msaConvert (alin2, type="seqinr::alignment")
a2

dist2 <- dist.alignment(a2, "identity")
dist2


### hacer los arboles

arbol1 <- nj (dist1)

pdf ("results/arbol ClustalW.pdf")
plot (arbol1)
dev.off() 

arbol2 <- nj (dist2)

pdf ("results/arbol Muscle.pdf")
plot (arbol2)
dev.off() 










