covarCheck <- function(phy){
	
	# get the VCV array
	V <- VCV.array(phy)
	
	# get the variance components from the diagonal
	var <- diag(V)
	
	# sweep these values out of the covariances - which should be smaller
	Vsweep <- sweep(V, 1, var)	
	
	# which covariance - variance values equal to zero
	zeros <- which(Vsweep == 0, arr.ind=TRUE)
	zeros <- apply(Vsweep, 1, function(X) colnames(V)[X ==0])

	# ditch single zeros and find unique sets of names
	zeros <- unique(zeros[sapply(zeros, length) > 1])
	
	return(zeros)
	
}

good <-  read.tree(text="((((t9:0.127910607,((t2:0.06831648015,t3:0.06831648015):0.05324626053,t4:0.1215627407):0.006347866363):0.1737651682,(t1:0.2938401338,t5:0.2938401338):0.00783564141):0.127341787,(t12:0.08350726193,t8:0.08350726193):0.3455103003):0.4600415532,((t7:0.02071294976,t10:0.02071294976):0.3000792894,(t6:0.1940654776,t11:0.1940654776):0.1267267615):0.5682668764);")

bad  <-  read.tree(text="((((t9:0.1071976573,((t2:0.04760353039,t3:0.04760353039):0.05324626053,t4:0.1008497909):0.006347866363):0.1737651682,(t1:0.2731271841,t5:0.2731271841):0.00783564141):0.127341787,(t12:0.06279431217,t8:0.06279431217):0.3455103003):0.4600415532,((t7:0,t10:0):0.3000792894,(t6:0.1733525279,t11:0.1733525279):0.1267267615):0.5682668764);")

vbad <- read.tree(text="((((t9:0.1071976573,(t2:0,t3:0,t4:0):0.1071976573):0.1737651682,(t1:0.2731271841,t5:0.2731271841):0.00783564141):0.127341787,(t12:0.06279431217,t8:0.06279431217):0.3455103003):0.4600415532,((t7:0,t10:0):0.3000792894,(t6:0.1733525279,t11:0.1733525279):0.1267267615):0.5682668764);")

covarCheck(good)
covarCheck(bad)
covarCheck(vbad)