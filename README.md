# AACluster
# Parallel Optimization of Categorical Clustering via an Adaptive Agglomerative Algorithm
# How to install  AAACluster:

# Go to the package directory
setwd("~/.../AACclust")

# Load devtools
library(devtools)

# Clean completely
devtools::clean_dll()

# Recompile
devtools::document()

#### Step 10: Use your package

### Close R and reopen it (to load the newly installed package):

library(AACclust)

# Load data
data(felinos)

# Run AAC with k=5
particion <- aac_cluster(felinos, alpha = 5)

print(particion)

# Calculate QI
qi <- qi_index(felinos, particion)

cat("QI =", qi, "\n")

# Search for optimal k automatically
resultado <- buscar_k_optimo(felinos)

cat("Optimal k =", resultado$k_optimo, "\n")

cat("Maximum QI =", resultado$QI_maximo, "\n")
