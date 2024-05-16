
grep("CCND3 Cyclin D3",golub.gnames[,2])
ccnd3_data <- golub[1042,]

#  Perform hierarchical clustering with single linkage & Ward linkage

gol.fac <- factor(golub.cl, levels = 0:1, labels = c("ALL", "AML"))
par(mfrow = c(1, 2))

# single linkage 
plot(ccnd3_data, pch = as.numeric(gol.fac))
legend("topright", legend = c("ALL", "AML"), pch = 1:2)
plot(hclust(dist(ccnd3_data, method = "euclidean"), method = "single"), labels = gol.fac)

# Ward linkage
plot(ccnd3_data, pch = as.numeric(gol.fac))
plot(hclust(dist(ccnd3_data, method = "euclidean"), method = "ward.D2"), labels = gol.fac)

cut_single <- cutree(hclust(dist(ccnd3_data, method = "euclidean"), method = "single"), k = 2)
cut_ward <- cutree(hclust(dist(ccnd3_data, method = "euclidean"), method = "ward.D2"), k = 2)
table(cut_single, gol.fac)
table(cut_ward, gol.fac)

#  k-means cluster analysis 
cl.mean <- kmeans(ccnd3_data, centers = 2, nstart = 10)
# Compare clusters with patient groups ALL/AML
table(cl.mean$cluster, gol.fac)

cl.mean <- kmeans(ccnd3_data, centers = 2, nstart = 10)
cl.mean

# 
initial <- cl.mean$centers
n <- dim(ccnd3_data)[1]
nboot <- 1000
boot.cl <- matrix(NA, nrow = nboot, ncol = 4)

for (i in 1:nboot) {
  dat.star <- ccnd3_data[sample(1:n+1, replace = TRUE), ]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot.cl[i, ] <- c(cl$centers[, 1], cl$centers[, 2])
}

mean_boot <- apply(boot.cl, 2, mean)
quantile(boot.cl[, 1], c(0.025, 0.975))
quantile(boot.cl[, 2], c(0.025, 0.975))

# 
k <- 1:30
sse <-rep(NA,length(K))

for(k in K){
  sse[k] <- kmeans(ccnd3_data,centers=k,nstart = 10)$tot.withinss
  
}
plot(K,sse,type='o',xaxt="n");axis(1,at=K,las = 2)

# Identify the elbow point
elbow_point <- which.min(diff(sse))
abline(v = elbow_point, col = "red", lty = 2)

cat("Suggested number of clusters:", elbow_point, "\n")




library(multtest)
data(golub)

# a) Select Oncogenes and Antigens:
sel1 <- grep("oncogene", golub.gnames[, 2])
sel2 <- grep("antigen", golub.gnames[, 2])

# b) Clustering Analysis with K-means and K-medoids:

clusdata <- rbind(golub[sel1, ], golub[sel2, ])
g.name <- rep(c("oncogene", "antigen"), c(length(sel1), length(sel2)))
kmeans_clusters <- kmeans(clusdata, centers = 2)
kmedoids_clusters <- pam(clusdata, k = 2)

# Compare clusters with gene groups
table(kmeans_clusters$cluster, g.name)
table(kmedoids_clusters$cluster, g.name)

# Test the Marginal Independence:

chi_square_kmeans <- chisq.test(table(kmeans_clusters$cluster, g.name))
chi_square_kmedoids <- chisq.test(table(kmedoids_clusters$cluster, g.name))

chi_square_kmeans
chi_square_kmedoids

# Plot Cluster Dendrograms

# Hierarchical clustering with single linkage
hc_single <- hclust(dist(clusdata), method = "single")
plot(hc_single, hang = -1, labels = g.name)

# Hierarchical clustering with complete linkage
hc_complete <- hclust(dist(clusdata), method = "complete")
plot(hc_complete,hang = -1, labels = g.name)



# K-means clustering and plot K versus SSE
K <-(1:30)
sse <- rep(NA,length(K))

for (k in K) {
  sse[k] <- kmeans(ncidata, centers = k,nstart = 10)$tot.withinss
  
}

# Plot K versus SSE
plot(K, sse, type = "o",xaxt='n',main = "K versus SSE")
axis(1,at=K,las=2)
# Identify the optimal K based on the plot
optimal_k <- which.min(sse)
abline(v = optimal_k, col = "blue", lty = 2)
text(optimal_k, sse[optimal_k], labels = paste("Optimal K =", optimal_k), pos = 4, col = "red")


# Perform K-medoids clustering with K=7 and 1-correlation as the dissimilarity measure
kmedoids_cluster <- pam(as.dist(1 - cor(t(ncidata))), k = 7)
clusters <- kmedoids_cluster$clustering
cluster_summary <- table(clusters, ncilabs)
cluster_summary

# Identify well-grouped and not-grouped cancers
well_grouped_cancer <- names(which.max(rowSums(cluster_summary)))
not_grouped_cancer <- names(which.max(colSums(cluster_summary)))
cat("Well-Grouped Cancer Type:", well_grouped_cancer, "\n")
cat("Not-Grouped Cancer Type:", not_grouped_cancer, "\n")

ovarian_cluster <- clusters[ncilabs=="OVARIAN"]
similar_cancers <- names(sort(table(ovarian_cluster),decreasing = TRUE)[2:3])
cat("Cancer types most similar to ovarian cancer:", paste(similar_cancers,collapse = ", "), "\n")
