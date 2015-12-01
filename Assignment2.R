#A <- mat.or.vec(2,2)
#B <- mat.or.vec(2,3)

#A[1,1] = 1
#A[2,1] = 2
#A[1,2] = 0
#A[2.2] = 2

#B[1,1] = 0
#B[1,2] = 4
#B[1,3] = 2
#B[2,1] = 1
#B[2,2] = 1
#B[2,3] = 2

#A%*%B
rawMatrix <- mat.or.vec(6, 6)
numArticles <- c(3, 2, 5, 1, 2, 1)
rowCount <- 6
a <- 0.85
e <- 0.0001

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

diagonalAdjustMatrix <- rawMatrix

for (i in 1:rowCount) {
  for (k in 1:rowCount) {
    if (i == k) {
      diagonalAdjustMatrix[i,k] <- 0
    }
  }
}

sums <- colSums(diagonalAdjustMatrix) #sums of columns

dangle <- sums

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
  normalizedMatrix <- diagonalAdjustMatrix
}

primeMatrix <- normalizedMatrix
sumOfArticles <- sum(numArticles)


for (i in 1:length(numArticles)) {
  numArticles[i] = numArticles[i] / sumOfArticles
}
sum <- sum(numArticles)

initialVector <- c() #inital start vector
for (i in 1:rowCount) {
  initialVector[i] = 1 / rowCount
}

for (i in 1:rowCount) {
  if (dangle[i] == 1) {
    primeMatrix[,i] <- numArticles
  }
}

primeMatrix

#calculating influence vector

aeT <- mat.or.vec(6, 6)

for (i in 1:rowCount) {
  for (k in 1:rowCount) {
    aeT[i,k] = numArticles[i]
  }
}

P <- (0.85*primeMatrix+(0.15)*aeT)