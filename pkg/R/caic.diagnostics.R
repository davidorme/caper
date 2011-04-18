caic.diagnostics <- function(caicObj, which.terms=NULL, which.plots=c("NV","SD","AGE"),
                             outlier.val=3, label=FALSE, ultrametric.tol=0.0001, 
                             ask=TRUE, test.signif=TRUE, plot.signif=TRUE, alpha=0.05, ...){
	
	    opar <- par(no.readonly=TRUE)
	    on.exit(par(opar))

	    if (ask) {
	        op <- par('ask')
			par(ask = TRUE)
	        on.exit(par(op))
	    }


	    which.plots <- match.arg(which.plots, c("NV","SD","AGE"), several.ok=TRUE)
    
	    if(is.null(which.terms)){
	        which.terms <- dimnames(attr(terms(formula(caicObj$mod)), "factors"))[[2]] # explanatory terms
	     } else {
	        if(! all(which.terms %in% dimnames(attr(terms(formula(caicObj$mod)), "factors"))[[2]])) stop("Not all specified terms were present in the model.")
	     }
     
 
     
	    # if(attr(caicObj, "contr.method") == "macrocaic"){
        
        
        
        
	    # } else {
        
	        tab <- caic.table(caicObj, nodalValues=TRUE, validNodes=TRUE, ultrametric.tol)
        
	        if(any(is.na(tab$nodeAge))) {
	            warning("Plots of absolute contrasts against node ages requested where node ages are not available.\n",
	                    "The tree may not be ultrametric, or the tolerance may need to be adjusted.")
	            which.plots <- which.plots[! which.plots=="AGE"]
	        }

	        outlier <- (tab$studResid >= outlier.val) + 20
        
	        tests <- list()
        
	        for(var in which.terms){
            
	            tests[[var]] <- list()
	            ylabExpr <-substitute(expression(abs(plain('Contrasts in') ~~ VAR)), env=list(VAR=as.name(var)))
            
	            # plot absolute contrasts against nodal values
	            if("NV" %in% which.plots){               
	                eval(substitute(plot(abs(CNT) ~ NV, data=tab, pch=outlier, ylab=ylabExpr, xlab=paste("Nodal values in", var), ...), 
	                     env=list(CNT= as.name(var), NV =as.name(paste("nodal.", as.character(var), sep="")))))
                     
	                if(test.signif){
	                    NVmod <- eval(substitute(lm(abs(CNT) ~ NV, data=tab), 
	                        env=list(CNT= as.name(var), NV =as.name(paste("nodal.", as.character(var), sep="")))))
                        
	                    tests[[var]][["NV"]] <- NVmod
                    
	                    if(plot.signif){
	                        if(summary.aov(NVmod)[[1]][1,5] <= alpha) abline(NVmod)
	                    }
	                 }
	             }
            
	            # plot absolute contrasts against sd
	           if("SD" %in% which.plots){
	               eval(substitute(plot(abs(CNT) ~ sqrt(contrVar), data=tab, pch=outlier, ylab=ylabExpr,
	                    xlab=expression(sqrt(plain('Variance at node'))), ...), env=list(CNT= as.name(var))))

	               if(test.signif){
	                    SDmod <- eval(substitute(lm(abs(CNT) ~ sqrt(contrVar), data=tab), env=list(CNT= as.name(var))))
                        
	                    tests[[var]][["SD"]] <- SDmod
                    
	                    if(plot.signif){
	                        if(summary.aov(SDmod)[[1]][1,5] <= alpha) abline(SDmod)
	                    }
	                 }
	            }
            
	           if("AGE" %in% which.plots){
	                 eval(substitute(plot(abs(CNT) ~ log(nodeAge), data=tab, pch=outlier, ylab=ylabExpr, xlab="Ln Node Age", ...), env=list(CNT= as.name(var))))
                 
	                 if(test.signif){
	                    AGEmod <- eval(substitute(lm(abs(CNT) ~ log(nodeAge), data=tab), env=list(CNT= as.name(var))))
                    
	                    tests[[var]][["AGE"]] <- AGEmod
                
	                    if(plot.signif){
	                        if(summary.aov(AGEmod)[[1]][1,5] <= alpha) abline(AGEmod)
	                    }
                 
	                }
	            }

	        }
        
	    # } 
    
	    class(tests) <- "caic.diagnostics"
	    return(tests)

}    
    
print.caic.diagnostics <- function(x, ...){
    
    for(vars in seq(along=x)){
        cat(names(x)[vars], ":\n")
        
        datf <- array(NA, dim=c(length(x[[vars]]), 4), dimnames=list(names(x[[vars]]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        
        for(tests in seq(along=x[[vars]])){
            datf[tests, ] <- coef(summary(x[[vars]][[tests]]))[2,]
        }
        
        print(datf)
        
    }
}

