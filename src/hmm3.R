
# Get the DNA data
source('~/Desktop/bgse/projects/github/smo/src/correct_dna.R')
res <- correct.dna(cut.matrix = FALSE)

# Data obtained
dna <- res[['dna']]
sequence <- res[['sequence']]
nucleotids <- res[['nucleotids']]

# Sequence of aminoacids
amino <- sapply(seq(1, length(nucleotids), 3), function(i) {
  paste(nucleotids[i:(i + 2)], collapse = '')
})

################################################################################
# HMM matrix
# Skeleton
hmm <- vector(mode = 'list', length = 5)
names(hmm) <- c('States', 'Symbols', 'startProbs', 'transProbs', 
                 'emissionProbs')

# Names of the states
hmm[[1]] <- c('E', 'I', 'N')

# Name of the aminoacids
hmm[[2]] <- unique(amino)

# Initial state probabilities (a priori)
tt <- table(substr(dna[, 1], 1, 1))
hmm[[3]] <- as.numeric(tt / sum(tt))
names(hmm[[3]]) <- c('E', 'I', 'N')

# Simulate transition probabilities
set.seed(666)
first <- sample(1:nrow(dna), 1)
labs <- substr(dna[, 1], 1, 1)[first]

# Conditional probabilities
pI <- tt[names(tt) %in% c('E', 'N')] / sum(tt[names(tt) %in% c('E', 'N')])
pE <- tt[names(tt) %in% c('I', 'N')] / sum(tt[names(tt) %in% c('I', 'N')])
pN <- tt[names(tt) %in% c('E', 'I')] / sum(tt[names(tt) %in% c('E', 'I')])

# Start simulation
for (i in 2:nrow(dna)) {
  if (labs[i - 1] == 'I') {
    draw <- ifelse(rbinom(1, 1, pI[1]) == 1, 'E', 'N')
  } else if (labs[i - 1] == 'E') {
    draw <- ifelse(rbinom(1, 1, pE[1]) == 1, 'I', 'N')
  } else if (labs[i - 1] == 'N') {
    draw <- ifelse(rbinom(1, 1, pN[1]) == 1, 'E', 'I')
  }
  labs <- c(labs, draw)
}

# Extend the sequence times 20 (one per aminoacid)
#ext.labs <- sapply(labs, function(x) { rep(x, 20) })
ext.labs <- c()
for (lab in labs) {
  ext.labs <- c(ext.labs, rep(lab, 20))
}

# Transition probabilities
trI <- table(ext.labs[which(ext.labs == 'I') + 1]) /
       sum(table(ext.labs[which(ext.labs == 'I') + 1]))
trE <- table(ext.labs[which(ext.labs == 'E') + 1]) /
       sum(table(ext.labs[which(ext.labs == 'E') + 1]))
trN <- table(ext.labs[which(ext.labs == 'N') + 1]) /
       sum(table(ext.labs[which(ext.labs == 'N') + 1]))
trN <- table(ext.labs[which(ext.labs == 'N') + 1]) /
       sum(table(ext.labs[which(ext.labs == 'N') + 1]))

# Fill the matrix
res <- matrix(nrow = 3, ncol = 3)
colnames(res) <- c('E', 'I', 'N')
rownames(res) <- c('E', 'I', 'N')
res[1, ] <- trE
res[2, ] <- trI
res[3, ] <- trN
# if (FALSE) {
#   res[1, ] <- c(0.95, 0.025, 0.025)
#   res[2, ] <- c(0.025, 0.95, 0.025)
#   res[3, ] <- c(0.025, 0.025, 0.95)
# }
hmm[[4]] <- res
names(attr(hmm$transProbs, "dimnames"))[1] <- 'from'
names(attr(hmm$transProbs, "dimnames"))[2] <- 'to'

# Emission probabilities
classes <- c()
for (i in 1:nrow(dna)) {
  classes <- c(classes, rep(dna[i, 1], 20))
}
amino2 <- cbind(amino, classes)
amino2[, 2] <- substr(amino2[, 2], 1, 1)
camino <- paste(amino2[, 2], amino2[, 1], sep = '')

# Fill the matrix
res <- matrix(nrow = 3, ncol = 64)
uamino <- unique(amino)
for (i in 1:64) {
  aux <- c()
  aux[1] <- length(which(camino == paste('E', uamino[i], sep = ''))) /
            length(which(substr(camino, 1, 1) == 'E'))
  aux[2] <- length(which(camino == paste('I', uamino[i], sep = ''))) /
            length(which(substr(camino, 1, 1) == 'I'))
  aux[3] <- length(which(camino == paste('N', uamino[i], sep = ''))) /
            length(which(substr(camino, 1, 1) == 'N'))
  res[, i] <- aux
}
colnames(res) <- uamino
rownames(res) <- c('E', 'I', 'N')
hmm[[5]] <- res
names(attr(hmm$emissionProbs, "dimnames"))[1] <- 'states'
names(attr(hmm$emissionProbs, "dimnames"))[2] <- 'symbols'
################################################################################

################################################################################
# Run it MOTHERFUCKER
total <- c()
for (i in 1:nrow(dna)) {
  j <- 20 * (i - 1) + 1
  obs2 <- amino[j:(j + 19)]
  trial <- HMM::viterbi(hmm, obs2)
  post1 <- rowMeans(HMM::posterior(hmm, obs2))
  post2 <- rowMeans(HMM::posterior(hmm, rev(obs2)))

  lab1 <- names(which.max(post1))
  lab2 <- names(which.max(post2))

  if (lab1 == lab2) {
    result <- lab1
  } else {
    result <- ifelse(which.max(c(max(post1), max(post2))) == 1, lab1, lab2)
  }

  total <- c(total, names(table(trial))[1])
}

table(dna[, 1], total)
sum(diag(table(dna[, 1], total))) / sum(table(dna[, 1], total))
(table(dna[, 1], total)[, 1] / rowSums(table(dna[, 1], total))[1])[1]
(table(dna[, 1], total)[, 2] / rowSums(table(dna[, 1], total))[2])[2]
(table(dna[, 1], total)[, 3] / rowSums(table(dna[, 1], total))[3])[3]



