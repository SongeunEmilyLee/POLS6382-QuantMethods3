## Heteroskedastic Normal Regression via ML in R
# Version 2. 
# Charles H. Franklin, Monday, June 30, 2003 at 19:37
# Modified: Thursday, July 1, 2004 at 20:19--Make y 1st argument
# Modified: Saturday, July 17, 2004 at 20:00

  MLhetreg<-function(y,X,Z,method='BFGS',Xnames=colnames(X),Znames=colnames(Z)){
    X<-cbind(1,X)
    colnames(X)[1]<-"Constant"
    nx<-ncol(X)
    Z<-cbind(1,Z)
    colnames(Z)[1]<-"ZConstant"
    nz<-ncol(Z)
    neglnl<-function(theta,X,Z,y){
        b<-theta[1:ncol(X)]
        g<-theta[ncol(X)+1:ncol(Z)]
        lnl<-as.vector(-.5*(Z%*%g)-(.5/exp(Z%*%g))*(y-X%*%b)^2)
        -sum(lnl)
        }        
    result<-c(optim(c(mean(y),rep(0,ncol(X)-1),log(var(y)),rep(0,ncol(Z)-1)), 
                neglnl, hessian=T, method=method, X=X, Z=Z, y=y),
                list(varnames=c(Xnames,Znames),nx=nx,nz=nz))
    class(result)<-"MLhetreg"        
    return(result)
    }
    
## Print method for MLhetreg

    print.MLhetreg <- function(object){
        coef<-object$par
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n hetreg converged\n')
        if (!object$convergence==0) cat('\n *** hetreg failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLhetreg

    summary.MLhetreg<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        nz<-object$nz
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\nHeteroskedastic Linear Regression\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        ghat<-coef[(nx+2):length(coef)]
        gvc<-vc[(nx+2):length(coef),(nx+2):length(coef)]
        wald<-t(ghat)%*%solve(gvc)%*%ghat
        pwald<-1-pchisq(wald,nz-1)
        cat("\nWald Test for Heteroskedasticity\n")
        cat("  Wald statistic: ",wald, "with ", nz-1," degrees of freedom\n")
        cat("                   p=",pwald,"\n")
        }
