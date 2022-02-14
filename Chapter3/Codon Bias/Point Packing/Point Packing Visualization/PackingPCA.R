library(factoextra)
packing <- read.csv(file.choose(), header = FALSE)
colnames(packing) <- c("AAA(K)","AAC(N)","AAG(K)","AAT(N)","ACA(T)","ACC(T)","ACG(T)","ACT(T)","AGA(R)","AGC(S)",
                        "AGG(R)","AGT(S)","ATA(I)","ATC(I)","ATT(I)","CAA(Q)","CAC(H)","CAG(Q)","CAT(H)","CCA(P)",
                        "CCC(P)","CCG(P)","CCT(P)","CGA(R)","CGC(R)","CGG(R)","CGT(R)","CTA(L)","CTC(L)","CTG(L)",
                        "CTT(L)","GAA(E)","GAC(D)","GAG(E)","GAT(D)","GCA(A)","GCC(A)","GCG(A)","GCT(A)","GGA(G)",
                        "GGC(G)","GGG(G)","GGT(G)","GTA(V)","GTC(V)","GTG(V)","GTT(V)","TAC(Y)","TAT(Y)","TCA(S)",
                        "TCC(S)","TCG(S)","TCT(S)","TGC(C)","TGT(C)","TTA(L)","TTC(F)","TTG(L)","TTT(F)")

p.pca <- prcomp(packing, scale = TRUE)

fviz_pca_ind(p.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping  
            )
p.eigvec <- p.pca$rotation[,1]
p.eigvec
