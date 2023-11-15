#cleanup
#rm(finalData, genes, GOfilter, names, all.features, all.genes, r, row, top10,
#   top10P, val, genescount, geneList, Object_Pooled.markers)
#Reference ARI and NMI
Ref <- clustComp(Object@meta.data$class_label, Object@meta.data$cluster_label)

#Comparing phase1 with phase2 in a table
P1_P2 <- clustComp(Idents(Object_Pooled), Idents(Object))

table <- matrix(ncol = 6)
table <- rbind(table, round(c(Obj_Class$ARI, Obj_Cluster$ARI, ObjPooled_Class$ARI,
                              ObjPooled_Cluster$ARI, P1_P2$ARI,Ref$ARI), digits = 3))
table <- rbind(table, round(c(Obj_Class$NMI, Obj_Cluster$NMI, ObjPooled_Class$NMI,
                              ObjPooled_Cluster$NMI, P1_P2$NMI, Ref$NMI), digits = 3))
table <- table[-1,]
rownames(table) <- c("ARI", "NMI")
colnames(table) <- c("P1-Class", "P1-Cluster", "P2-Class",
                     "P2-Cluster", "P1-P2", "Class-Cluster(Reference)")
table

#Comparison via alluvial plot
rm(alluvium)
alluvium <- data.frame(Idents(Object), Idents(Object_Pooled))
alluvium$P1 <- as.numeric(alluvium$Idents.Object.)
alluvium$P2 <- as.numeric(alluvium$Idents.Object_Pooled.)
alluvium$frequency <- alluvium$P1 / nrow(alluvium)
alluvium <- alluvium[,-1:-2]
#alluvium

ggplot(alluvium,
       aes(y = frequency,
           axis1 = P1, axis2 = P2)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_alluvium(aes(fill = P1), width = 1/12, segments = 450, curve_type = "sine") +
  guides(fill = "none") +
  geom_stratum(alpha = .25, reverse = TRUE, width = 1/12) +
  scale_fill_distiller(palette = "Set2") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("P1", "P2"), labels = c("P1", "P2"), expand = c(0.05, 0.05)) +
  ggtitle("Phase 1 to Phase 2 clustering")
