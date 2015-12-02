#Smyth May
#INFO 370 Assignment 2
#2015/12/01

#rawMatrix for Assignment 2 example input
rawMatrix <- mat.or.vec(6, 6)
numArticles <- c(3, 2, 5, 1, 2, 1) #Article number per source
rowCount <- 6 #current number of rows
e <- 0.0001

#Insert data into matrices
rawMatrix[1,1] = 1
rawMatrix[1,2] = 0
rawMatrix[1,3] = 2
rawMatrix[1,4] = 0
rawMatrix[1,5] = 4
rawMatrix[1,6] = 3
rawMatrix[2,1] = 3
rawMatrix[2,2] = 0
rawMatrix[2,3] = 1
rawMatrix[2,4] = 1
rawMatrix[2,5] = 0
rawMatrix[2,6] = 0
rawMatrix[3,1] = 2
rawMatrix[3,2] = 0
rawMatrix[3,3] = 4
rawMatrix[3,4] = 0
rawMatrix[3,5] = 1
rawMatrix[3,6] = 0
rawMatrix[4,1] = 0
rawMatrix[4,2] = 0
rawMatrix[4,3] = 1
rawMatrix[4,4] = 0
rawMatrix[4,5] = 0
rawMatrix[4,6] = 1
rawMatrix[5,1] = 8
rawMatrix[5,2] = 0
rawMatrix[5,3] = 3
rawMatrix[5,4] = 0
rawMatrix[5,5] = 5
rawMatrix[5,6] = 2
rawMatrix[6,1] = 0
rawMatrix[6,2] = 0
rawMatrix[6,3] = 0
rawMatrix[6,4] = 0
rawMatrix[6,5] = 0
rawMatrix[6,6] = 0

#copy the raw matrix so we can keep it
diagonalAdjustMatrix <- rawMatrix

#make the diagonal 0
for (i in 1:rowCount) {
  for (k in 1:rowCount) {
    if (i == k) {
      diagonalAdjustMatrix[i,k] <- 0
    }
  }
}

sums <- colSums(diagonalAdjustMatrix) #sums of columns

#copy of sums to make dangling vector
dangle <- sums

#Check for dangling nodes and normalize the diagonal matrix
for (i in 1:rowCount) {
  for (k in 1:rowCount) {
    if (diagonalAdjustMatrix[i,k] != 0) {
      diagonalAdjustMatrix[i,k] <- (diagonalAdjustMatrix[i,k]/sums[k])
    } 
  }
  if (sums[i] != 0) { #dangling nodes
    dangle[i] <- 0
  } else {
    dangle[i] <- 1
  }
  #official normalized matrix
  normalizedMatrix <- diagonalAdjustMatrix
}

#H' for P equation, which we never use...
primeMatrix <- normalizedMatrix
#total number of articles
sumOfArticles <- sum(numArticles)

#normalize the article vector
for (i in 1:length(numArticles)) {
  numArticles[i] = numArticles[i] / sumOfArticles
}

initialVector <- mat.or.vec(6,1) #initial start vector
#Normalize the initial vector a
for (i in 1:rowCount) {
  initialVector[i] = 1 / rowCount
}

#Make H' by replacing dangling nodes.
for (i in 1:rowCount) {
  if (dangle[i] == 1) {
    primeMatrix[,i] <- numArticles
  }
}

#Define variables for iterative solving.
#Done in one place to keep things simple and readable
alpha <- 0.85 
epsilon <- 0.0001 
H <- normalizedMatrix 
pi_K <- initialVector
d <- dangle
a <- numArticles

# constructs the equation åHπ^k + [å.d.π^k + (1 - å)]a
pi_K1 <- (alpha*H%*%pi_K) + ((alpha*(d%*%pi_K) + 0.15)*a)

#first residual value, on the offchance the first iteration is correct...
residual <- sum(abs(pi_K1 - pi_K))

#continue with the iterative solving until Epsilon is met.
while (residual > epsilson) {
  pi_K <- pi_K1 #set pi^k+1 as the current initial start vector
  pi_K1 = (alpha*H%*%pi_K) + ((alpha*(d%*%pi_K) + 0.15)*a)
  residual <- sum(abs(pi_K1 - pi_K)) #calculate new residual
}

#final influence vector defined
influence <- pi_K1

#find the value we need to multiply by 100...
efEquation <- ((H%*%influence)/(colSums(H%*%influence)))
ef <- 100*efEquation
ef

#Answer for example data set...
#        [,1]
#[1,] 31.938649
#[2,] 20.031575
#[3,] 11.984812
#[4,]  3.340727
#[5,] 32.704237
#[6,]  0.000000