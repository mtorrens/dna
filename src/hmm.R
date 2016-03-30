source('~/Desktop/bgse/projects/github/smo/src/correct_dna.R')

dna <- res[['dna']]
sequence <- res[['sequence']]
nucleotids <- res[['nucleotids']]

priors <- table(nucleotids) / sum(table(nucleotids))
poster <- as.data.frame(matrix(ncol = 4, nrow = 4))
colnames(poster) <- c('A', 'C', 'G', 'T')
rownames(poster) <- c('A', 'C', 'G', 'T')

tt <- table(nucleotids[which(nucleotids == 'A') + 1])
poster[1, ] <- tt / sum(tt)

tt <- table(nucleotids[which(nucleotids == 'C') + 1])
poster[2, ] <- tt / sum(tt)

tt <- table(nucleotids[which(nucleotids == 'G') + 1])
poster[3, ] <- tt / sum(tt)

tt <- table(nucleotids[which(nucleotids == 'T') + 1])
poster[4, ] <- tt / sum(tt)

library(tileHMM)
library(HMM)

priorsS <- table(dna[, 1]) / sum(table(dna[, 1]))
posterS <- as.data.frame(matrix(ncol = 3, nrow = 3))

tt <- table(dna[which(dna[, 1] == 'EI') + 1, 1])
posterS[1, ] <- tt / sum(tt)

tt <- table(dna[which(dna[, 1] == 'IE') + 1, 1])
tt <- table(dna[which(dna[, 1] == 'N') + 1, 1])


####
{
  tt <- table(dna[, 1])
  mhmm <- vector(mode = 'list', length = 5)
  names(mhmm) <- c('States', 'Symbols', 'startProbs', 'transProbs', 'emissionProbs')
  mhmm[[1]] <- c('E', 'I', 'N')
  mhmm[[2]] <- c('A', 'C', 'G', 'T')
  mhmm[[3]] <- as.numeric(tt / sum(tt))
  names(mhmm[[3]]) <- c('E', 'I', 'N')
  res <- matrix(nrow = 3, ncol = 3)
  colnames(res) <- c('E', 'I', 'N')
  rownames(res) <- c('E', 'I', 'N')
  res[, 1] <- c(0, tt[2], tt[3]) / sum(tt[2] + tt[3])
  res[, 2] <- c(tt[1], 0, tt[3]) / sum(tt[1] + tt[3])
  res[, 3] <- c(tt[1], tt[2], 0) / sum(tt[1] + tt[2])
  mhmm[[4]] <- res

  res2 <- res
  #attr(res, 'dimnames')[1] <- 'from'

  classes <- c()
  for (i in 1:nrow(dna)) {
    classes <- c(classes, rep(dna[i, 1], 60))
  }
  nucs2 <- cbind(nucleotids, classes)
  nucs2[, 2] <- substr(nucs2[, 2], 1, 1)
  cnucs <- paste(nucs2[, 2], nucs2[, 1], sep = '')

  res <- matrix(nrow = 3, ncol = 4)
  res[, 1] <- table(cnucs[grepl('A', cnucs)]) / sum(table(cnucs[grepl('A', cnucs)]))
  res[, 2] <- table(cnucs[grepl('C', cnucs)]) / sum(table(cnucs[grepl('C', cnucs)]))
  res[, 3] <- table(cnucs[grepl('G', cnucs)]) / sum(table(cnucs[grepl('G', cnucs)]))
  res[, 4] <- table(cnucs[grepl('T', cnucs)]) / sum(table(cnucs[grepl('T', cnucs)]))
  colnames(res) <- c('A', 'C', 'G', 'T')
  rownames(res) <- c('E', 'I', 'N')
  mhmm[[5]] <- res

  names(attr(mhmm$transProbs, "dimnames"))[1] <- 'from'
  names(attr(mhmm$transProbs, "dimnames"))[2] <- 'to'

  names(attr(mhmm$emissionProbs, "dimnames"))[1] <- 'states'
  names(attr(mhmm$emissionProbs, "dimnames"))[2] <- 'symbols'
}
#####

# hmm = initHMM(c('E', 'I', 'N'), c('A', 'C', 'G', 'T'), startProbs=as.numeric(tt / sum(tt)),
#      transProbs=res, emissionProbs=res2)

obs2 = nucleotids[1:60]

result <- viterbi(hmm, observations)
result2 <- viterbi(mhmm, obs2)

# for (ct in names(table(cnucs))) {
#   aux <- which(cnucs == ct) + 1

#   table(cnucs[aux])
# }


# ext.labs <- c()
# labs <- substr(dna[, 1], 1, 1)
# labs2 <- labs
# tt <- table(labs)
# set.seed(666)
# #ord <- sample(1:nrow(dna), 1)
# ch <- sample(1:length(labs2), 1)
# ext.labs <- c(ext.labs, rep(labs[ord], 60))
# labs2 <- labs2[-ord]
# if (ch == 'I') {
#   labs3 <- labs2[labs2 != 'I']
# } else if (ch == 'E') {
#   labs3 <- labs2[labs2 != 'E']
# } else if (ch == 'E') {
#   labs3 <- labs2[labs2 != 'I']
# }

or <- substr(dna[, 1], 1, 1)
set.seed(666)
first <- sample(1:nrow(dna), 1)
labs <- or[first]

pI <- tt[names(tt) %in% c('E', 'N')] / sum(tt[names(tt) %in% c('E', 'N')])
pE <- tt[names(tt) %in% c('I', 'N')] / sum(tt[names(tt) %in% c('I', 'N')])
#pN <- tt[names(tt) %in% c('E', 'N')] / sum(tt[names(tt) %in% c('E', 'N')])
#pN <- tt / sum(tt)
pN <- tt[names(tt) %in% c('E', 'I')] / sum(tt[names(tt) %in% c('E', 'I')])

for (i in 2:nrow(dna)) {
  if (labs[i - 1] == 'I') {
    draw <- rbinom(1, 1, pI[1])
    end <- ifelse(draw == 1, 'E', 'N')
  } else if (labs[i - 1] == 'E') {
    draw <- rbinom(1, 1, pE[1])
    end <- ifelse(draw == 1, 'I', 'N')
  } else if (labs[i - 1] == 'N') {
    draw <- rbinom(1, 1, pN[1])
    end <- ifelse(draw == 1, 'E', 'I')
  }
  labs <- c(labs, end)
}

# # Or...
# labs <- sample(or, length(or))

ext.labs <- c()
for (lab in labs) {
  ext.labs <- c(ext.labs, rep(lab, 60))
}

trI <- table(ext.labs[which(ext.labs == 'I') + 1]) / sum(table(ext.labs[which(ext.labs == 'I') + 1]))
trE <- table(ext.labs[which(ext.labs == 'E') + 1]) / sum(table(ext.labs[which(ext.labs == 'E') + 1]))
trN <- table(ext.labs[which(ext.labs == 'N') + 1]) / sum(table(ext.labs[which(ext.labs == 'N') + 1]))
#trN <- table(ext.labs) / sum(table(ext.labs))
#trN <- c(trN[1], 0, trN[2])
#names(trN)[2] <- 'I'

res <- matrix(nrow = 3, ncol = 3)
colnames(res) <- c('E', 'I', 'N')
rownames(res) <- c('E', 'I', 'N')
res[1, ] <- trE
res[2, ] <- trI
res[3, ] <- trN
mhmm[[4]] <- res

names(attr(mhmm$transProbs, "dimnames"))[1] <- 'from'
names(attr(mhmm$transProbs, "dimnames"))[2] <- 'to'

result2 <- HMM::viterbi(mhmm, obs2)

total <- c()
for (i in 1:nrow(dna)) {
  j <- 60 * (i - 1) + 1
  obs2 <- nucleotids[j:(j + 59)]
  trial <- HMM::viterbi(mhmm, obs2)
  total <- c(total, names(table(trial))[1])
}

table(dna[, 1], total)
sum(diag(table(dna[, 1], total))) / sum(table(dna[, 1], total))






















