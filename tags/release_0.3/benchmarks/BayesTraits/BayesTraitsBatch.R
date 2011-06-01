nul <- system('echo "6\n1\nkappa 0\nlambda 0\ndelta 0\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
fix <- system('echo "6\n1\nkappa 0.5\nlambda 0.5\ndelta 0.5\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
kld <- system('echo "6\n1\nkappa 1\nlambda 1\ndelta 1\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
Kld <- system('echo "6\n1\nkappa  \nlambda 1\ndelta 1\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
kLd <- system('echo "6\n1\nkappa 1\nlambda  \ndelta 1\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
klD <- system('echo "6\n1\nkappa 1\nlambda 1\ndelta  \nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
kLD <- system('echo "6\n1\nkappa 1\nlambda  \ndelta  \nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
KlD <- system('echo "6\n1\nkappa  \nlambda 1\ndelta  \nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
KLd <- system('echo "6\n1\nkappa  \nlambda  \ndelta 1\nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
KLD <- system('echo "6\n1\nkappa  \nlambda  \ndelta  \nrun\n" | ./BayesTraits test.nex test.dat', intern=TRUE)
                            
mods <- list(nul=nul[57], fix=fix[57], kld = kld[57], Kld = Kld[57], kLd = kLd[57],
	         klD = klD[57], kLD = kLD[57], KlD = KlD[57], KLd = KLd[57], KLD = KLD[57])

titles <- list(nul[56], fix[56], kld[56], Kld[56], kLd[56], klD[56], kLD[56], KlD[56], KLd[56], KLD[56])
all(sapply(titles, '==', titles[[1]]))

modsep <- sapply(mods, strsplit, split='\t')
modsep <- sapply(modsep, as.numeric)

# oddities in row names - 1 extra name and actually repeats Rsq
modsep <- t(modsep)[,-6]
cols <- strsplit(kld[56], split='\t')[[1]][-7:-8]

colnames(modsep) <- cols
BayesTraitsMods <- as.data.frame(modsep)
BayesTraitsMods$model <- names(mods)
BayesTraitsMods <- BayesTraitsMods[, c(10, 2:7, 9,8)]

save(BayesTraitsMods, file='../../data/benchBayesTraitsOutputs.rda')