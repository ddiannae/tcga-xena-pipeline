nIndex =  which(normalized.cqn$Targets$Group %in% "N")
tIndex = which(normalized.cqn$Targets$Group %in% "T")
rownames(normalized.cqn$M) <- normalized.cqn$Annot$Symbol
nExpression = normalized.cqn$M[,nIndex]
tExpression = normalized.cqn$M[,tIndex]
write.table(nExpression, file="NormalExpression_geneNames.tsv", quote=FALSE, sep="\t")
write.table(tExpression, file="TumorExpression_geneNames.tsv", quote=FALSE, sep="\t")

