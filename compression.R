# Packages ----
library(EMD) # Lena's image
library(blockmatrix) # For creating block matrices
library(ggplot2) # Graphics
library(redR) # For MSE and PSNR measures
library(SPUTNIK) # For SSIM measure

# Loading data 
data(lena) # Storage the image
image(lena, col=gray(0:100/100), axes=FALSE) # See the image

# Functions ----
blockToObs <- function(block_matrix) {
  obs_matrix <- c()
  for (i in 1:64) {
    for (j in 1:64) {
      obs <- as.vector(t(block_matrix[i,j]))
      obs_matrix <- rbind(obs_matrix, obs)
    }  
  }
  return(obs_matrix)
}

obsToBlock <- function(Obs_matrix) {
  matrix <- matrix(0, nrow = 512, ncol = 512) 
  cont <- 1
  for (i in seq(1,512,8)) {
    for (j in seq(1,512,8)) {
      block_matrix <- matrix(Obs_matrix[cont,], nrow = 8, ncol = 8, byrow = T)
      matrix[i:(i+8-1), j:(j+8-1)] <- block_matrix
      cont <- cont + 1
    }
  }
  return(matrix)
}

Tmatrix <- function(m) {
  t_m <- diag(0, ncol = 64, nrow = 64)
  for (i in 1:m) {
    t_m[i,i] <- 1
  }
  return(t_m)
}

pca <- function(data) {
  lena_block <- as.blockmatrix(data, nrowe = 8, ncole = 8) # Create a block matrix for image data
  x <- blockToObs(lena_block) # Transform a block matrix into a "dataframe" (blocks by row)
  x_cent <- scale(x, center = T, scale = F) # Centering the matrix 
  cov_matrix <- t(x_cent) %*% x_cent # Calculating the Covariance matrix
  
  eigen_pairs <- eigen(cov_matrix) # Calculating eigenvalues and eigenvectors
  eigenvalues <- eigen_pairs$values # Accessing eigenvalues
  eigenvectors <- eigen_pairs$vectors # Accessing eigenvectors
  
  results <- list("x" = x, "eigenvalues" = eigenvalues, "eigenvectors" = eigenvectors)
  return(results)
}

compression <- function(x, eigenvectors, m) {
  y <- x %*% eigenvectors # Creating the new "dataframe" with new p principal components 
  t_m <- Tmatrix(m) # T matrix for define the number of p taking into account
  y_m <- y %*% t_m # Creating the new "dataframe" with new p = dim(t_m) principal components
  
  b <- y_m %*% t(eigenvectors) # Inverse of the transformation
  B <- obsToBlock(b) # The image generating matrix
  return(B) # Returning the image matrix
}

# Compression ----
PCA <- pca(lena)
image <- compression(PCA$x, PCA$eigenvectors, 2) # See the image with m PC

# Graphics ----

# Eigenvalues' bar graphic 
pca_prop_var <- PCA$eigenvalues/sum(PCA$eigenvalues)
prop_var <- as.data.frame(pca_prop_var[1:8])
colnames(prop_var) = "prop_var"

p <- ggplot(prop_var, aes(x=c(1:8), y=prop_var)) + 
  geom_bar(stat="identity", fill="#494d4d") + 
  labs(title="", x="Principal Component", 
       y="% of Variance") + geom_line() + geom_point() + 
  scale_x_continuous(breaks=c(1:8))
p

# Matrices of the eigenvectors
par(mfrow = c(3,3))
image(matrix(PCA$eigenvectors[, 1], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 4], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 8], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 25], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 30], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 35], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 55], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 60], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)
image(matrix(PCA$eigenvectors[, 64], 8, 8, byrow = T), col=gray(0:100/100), axes=FALSE)

# Lena's plots
lena_1 <- compression(PCA$x, PCA$eigenvectors, 1) 
lena_3 <- compression(PCA$x, PCA$eigenvectors, 3) 
lena_7 <- compression(PCA$x, PCA$eigenvectors, 7) 
lena_10 <- compression(PCA$x, PCA$eigenvectors, 10) 
lena_13 <- compression(PCA$x, PCA$eigenvectors, 13) 

image(lena_1, col=gray(0:100/100), axes=FALSE) # See the image compression for 1 PC
image(lena_3, col=gray(0:100/100), axes=FALSE) # See the image compression for 3 PC
image(lena_7, col=gray(0:100/100), axes=FALSE) # See the image compression for 7 PC
image(lena_10, col=gray(0:100/100), axes=FALSE) # See the image compression for 10 PC
image(lena_13, col=gray(0:100/100), axes=FALSE) # See the image compression for 13 PC
image(lena, col=gray(0:100/100), axes=FALSE) # See the original image

# Quality of reconstruction ----
# PSNR
psnr_lena_1 <- PSNR(lena, lena_1)
psnr_lena_3 <- PSNR(lena, lena_3)
psnr_lena_7 <- PSNR(lena, lena_7)
psnr_lena_10 <- PSNR(lena, lena_10)
psnr_lena_13 <- PSNR(lena, lena_13)

#SSIM
ssim <- c()
for (i in 1:64){
  ssim[i] <- SSIM(lena, compression(PCA$x, PCA$eigenvectors, i)) 
}

ssim_lena_1 <- ssim[1]
ssim_lena_3 <- ssim[3]
ssim_lena_7 <- ssim[7]
ssim_lena_10 <- ssim[10]
ssim_lena_13 <- ssim[13]

q <- ggplot(data=as.data.frame(ssim), aes(x=1:64, y=ssim, group=1)) +
     geom_line()


# MSE
mse <- c()
for (i in 1:64){
  mse[i] <- MSE(lena, compression(PCA$x, PCA$eigenvectors, i)) 
}

mse_lena_1 <- mse[1]
mse_lena_3 <- mse[3]
mse_lena_7 <- mse[7]
mse_lena_10 <- mse[10]
mse_lena_13 <- mse[13]

r <- ggplot(data=as.data.frame(mse), aes(x=1:64, y=mse, group=1)) +
  geom_line()
