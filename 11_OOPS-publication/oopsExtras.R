
# GO enrichment
go.en = as.data.frame(enricher(gene = three.oops,TERM2GENE = univ.go[,c(2,1)],universe = u2os.univ))
goterms <- Term(GOTERM)
ont<- Ontology(GOTERM)
go.en$Description = goterms[go.en$ID]
go.en$Ontology = ont[go.en$ID]

go.mf = enricher(gene = three.oops,TERM2GENE = univ.go[,c(2,1)],universe = u2os.univ)
go.cc = enrichGO(gene = three.oops,OrgDb = org.Hs.eg.db,keytype = "UNIPROT",ont = "CC",qvalueCutoff = 0.05, universe = u2os.univ )
go.bp = enrichGO(gene = three.oops,OrgDb = org.Hs.eg.db,keytype = "UNIPROT",ont = "BP",qvalueCutoff = 0.05, universe = u2os.univ )

three.go = rbind(cbind(Ontology = "MF",data.frame(go.mf)),cbind(Ontology = "BP",data.frame(go.bp)),cbind(Ontology = "CC", data.frame(go.cc)))
one.go$Symbol = sapply(all.go$geneID, function(x) paste(bitr(unlist(strsplit(x,"/")),fromType="UNIPROT", toType=c("SYMBOL"), OrgDb=org.Hs.eg.db)$SYMBOL,collapse="/"))

# Plots showing enrichment
dotplot(go.cc,title = "GO:CC", showCategory = 20)
dotplot(go.bp, title = "GO:BP",showCategory = 20)

write.table(one.go,paste(outdir,"GO-annotation-for-proteins-found-in-a-single-OOPS-experiment.txt",sep="/"),sep="\t",row.names=F,quote=F)

# KEGG enrichment
one.entrez = bitr(geneID = one.oops ,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F)$ENTREZID
kg = enrichKEGG(one.entrez,qvalueCutoff = 0.05)
#browseKEGG(kg, 'hsa04144')