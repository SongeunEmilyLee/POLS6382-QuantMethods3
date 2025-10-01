##MLLibrary--- ML Routines for chf
## Version 3.0
## (c) 2004-2006 by Charles H. Franklin
## Latest Revision: Monday, July 19, 2004 at 13:53
##                  Tuesday, July 4, 2006 at 13:57
# Contents:
#    MLcauchy3.r
#    MLcloglog3.r
#    MLefx.r
#    MLhetReg.r
#    MLlogit3.r
#    MLloglog3.r
#    MLlrtest.r
#    MLmnl3.r
#    MLnegbin3.r
#    MLologit3.r
#    MLoprobit3.r
#    MLpoisson3.r
#    MLProbit3.r
#    MLwald.r
#    MLzip3.r
#
# Cauchy Regression via ML in R
# Version 2. 
# Charles H. Franklin, Monday, July 19, 2004 at 10:47

  MLcauchy<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
   
    llf0<-log(mean(y))*sum(y) + log(1-mean(y))*sum(1-y)
    neglnl<-function(theta,X,y){
            p<-as.vector(atan(X%*%theta)/pi+.5)
            lnl<-(y*log(p) + (1-y)*log(1-p))
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
        part1<-as.vector(y*(1/(pi^(-1)*atan(X%*%theta)+.5))*(pi^(-1)/(1+(X%*%theta)^2)) +
                        (1-y)*(1/(.5-pi^(-1)*atan(X%*%theta)))*(-(pi^(-1)/(1+(X%*%theta)^2))))
        -apply(part1*X,2,sum)
   }
    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    se0<-sqrt(sum((y-X%*%b0)^2)/(length(y)-nx))
    theta<-b0/se0*1.6
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),nx=nx,llf0=llf0))
    result$fitted.values<-as.vector(atan(X%*%result$par)/pi+.5)
    result$xb<-X%*%result$par
    class(result)<-"MLcauchy"        
    return(result)
    }
    
## Print method for MLcauchy

    print.MLcauchy <- function(object){
        coef<-object$par
        print(coef)
        if (object$convergence==0) cat('\n MLcauchy converged\n')
        if (!object$convergence==0) cat('\n *** MLcauchy failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLcauchy

    summary.MLcauchy<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Cauchy\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*(object$llf0--object$value)
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
# Complementary Log-log Regression via ML in R
# Version 2. 
# Charles H. Franklin, Wednesday, July 7, 2004 at 14:22

  MLcloglog<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
   
    llf0<-(log(mean(y))*sum(y) + log(1-mean(y))*sum(1-y))
    neglnl<-function(theta,X,y){
            p<-as.vector(1-exp(-exp(X%*%theta)))
            lnl<-y*log(p) + (1-y)*log(1-p)
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
        emexb<-exp(-exp(X%*%theta))       
        part1<-as.vector(((y*(-emexb)/(1-emexb))+(1-y))*-exp(X%*%theta))   
        -apply(part1*X,2,sum)
   }
    

##Bad starting values, but seem to work, eventually. 
## Problems when p->1 or 0 with linear start values.

    theta<-rep(0,len=nx)
    
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),nx=nx,llf0=llf0))
    result$fitted.values<-as.vector(1-exp(-exp(X%*%result$par)))
    result$xb<-X%*%result$par
    class(result)<-"MLcloglog"        
    return(result)
    }
    
## Print method for MLcloglog

    print.MLcloglog <- function(object){
        coef<-object$par
        print(coef)
        if (object$convergence==0) cat('\n MLcloglog converged\n')
        if (!object$convergence==0) cat('\n *** MLcloglog failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLcloglog

    summary.MLcloglog<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Complementary Log-log model\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*(object$llf0--object$value)
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
##MLefx --- Marginal Effects Routines for MLlibrary
## by chf Saturday, July 17, 2004 at 13:47
##  added MLefx.poisson Tuesday, July 17, 2007 at 21:09
##

MLefx.oprobit<-function(object,atx=NULL,aty=NULL){
    coef<-object$par
    names(coef)<-object$varnames
    nx<-object$nx
    J<-object$J
    b<-coef[1:nx]
    tau<-c(0,coef[(nx+1):length(coef)])
    
    if (!is.null(aty)){
        pj<-table(aty)/length(aty)
        cpj<-cumsum(pj)
        fx<-rep(0,J)
        for (j in 1:J){
            if (j==1) 
                fx[j]<--dlogis(qlogis(cpj[j]))
            else {
                if (j==J) 
                    fx[j]<-dlogis(qlogis(cpj[j-1]))
                else
                    fx[j]<--(dlogis(qlogis(cpj[j]))-dlogis(qlogis(cpj[j-1])))
                }
            }
        efx<-outer(coef[2:nx],fx)
        }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        xb<-atx%*%b
        fx<-rep(0,J)
        for (j in 1:J){
            if (j==1) 
                fx[j]<--dlogis(-xb)
            else {
                if (j==J) 
                    fx[j]<-dlogis(tau[j-1]-xb)
                else
                    fx[j]<--(dlogis(tau[j]-xb)-dlogis(tau[j-1]-xb))
                }
            }   
        efx<-outer(coef[2:nx],fx)
        }
        else stop("Error in MLefx: must provide either atx or aty values")
    efx<-cbind(coef[2:nx],efx)    
    dimnames(efx)<-list(object$varnames[2:nx],c("Coef",paste("Y=",1:J,sep="")))
    cat("\n Marginal Effects for Ordered Probit\n")
    if (!is.null(aty)) cat("\n Evaluated at P(y) =",round(pj,3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X =",round(atx,3),"\n")
    print(round(efx,4))
    } 
    


MLefx.ologit<-function(object,atx=NULL,aty=NULL){
    coef<-object$par
    names(coef)<-object$varnames
    nx<-object$nx
    J<-object$J
    b<-coef[1:nx]
    tau<-c(0,coef[(nx+1):length(coef)])
    
    if (!is.null(aty)){
        pj<-table(aty)/length(aty)
        cpj<-cumsum(pj)
        fx<-rep(0,J)
        for (j in 1:J){
            if (j==1) 
                fx[j]<--dlogis(qlogis(cpj[j]))
            else {
                if (j==J) 
                    fx[j]<-dlogis(qlogis(cpj[j-1]))
                else
                    fx[j]<--(dlogis(qlogis(cpj[j]))-dlogis(qlogis(cpj[j-1])))
                }
            }
        efx<-outer(coef[2:nx],fx)
        }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        xb<-atx%*%b
        fx<-rep(0,J)
        for (j in 1:J){
            if (j==1) 
                fx[j]<--dlogis(-xb)
            else {
                if (j==J) 
                    fx[j]<-dlogis(tau[j-1]-xb)
                else
                    fx[j]<--(dlogis(tau[j]-xb)-dlogis(tau[j-1]-xb))
                }
            }   
        efx<-outer(coef[2:nx],fx)
        }
        else stop("Error in MLefx: must provide either atx or aty values")
    efx<-cbind(coef[2:nx],efx)    
    dimnames(efx)<-list(object$varnames[2:nx],c("Coef",paste("Y=",1:J,sep="")))
    cat("\n Marginal Effects for Ordered Logit\n")
    if (!is.null(aty)) cat("\n Evaluated at P(y) =",round(pj,3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X =",round(atx,3),"\n")
    print(round(efx,4))
    } 
    

MLefx.mnl<-function(object,atx=NULL,aty=NULL){
    nx<-object$nx
    J<-object$J
    b<-cbind(0,matrix(object$par,nx,J-1))

    if (!is.null(aty)){
        pj<-table(aty)/length(aty)
        efx<-matrix(pj,nx,J,byrow=TRUE)*(b-matrix(b%*%pj,nx,J))
    }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        exb<-exp(atx%*%b)
        denom<-sum(exb)
        pj<-exb/denom
        efx<-matrix(pj,nx,J,byrow=TRUE)*(b-matrix(b%*%t(pj),nx,J))        
    }
        else stop("Error in MLefx: must provide either atx or aty values")
    dimnames(efx)<-list(object$varnames[1:nx],paste("Y=",1:J,sep=""))    
    cat("\n Marginal Effects for Multinomial Logit\n")
    if (!is.null(aty)) cat("\n Evaluated at P(y) =",round(pj,3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X = ",round(atx,3),"\n")
    if (!is.null(atx)) cat(" Predicted P(Y=j)=",round(pj,3),"\n")
    print(round(efx,4))
    } 


MLefx.probit<-function(object,atx=NULL,aty=NULL){
    nx<-object$nx
    b<-object$par
    
    if (!is.null(aty)){
        pj<-table(aty)/length(aty)
        efx<-dnorm(qnorm(pj[2]))*b
        }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        xb<-as.vector(atx%*%b)
        pj<-c(pnorm(-xb),pnorm(xb))
        efx<-dnorm(xb)*b
        }
        else stop("Error in MLefx: must provide either atx or aty values")
    efx<-cbind(b,efx,b*dnorm(0))
    dimnames(efx)<-list(object$varnames,c("Coef","  d(Y=1)/dx","  max(dy/dx)"))    
    cat("\n Marginal Effects for Binary Probit\n")
    if (!is.null(aty)) cat("\n Evaluated at P(y) =",round(pj[2],3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X = ",round(atx,3),"\n")
    if (!is.null(atx)) cat(" Predicted P(Y=1)=",round(pj[2],3),"\n")
    print(round(efx,4))
    }         


MLefx.logit<-function(object,atx=NULL,aty=NULL){
    nx<-object$nx
    b<-object$par
    
    if (!is.null(aty)){
        pj<-table(aty)/length(aty)
        efx<-dlogis(qlogis(pj[2]))*b
        }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        xb<-as.vector(atx%*%b)
        pj<-c(plogis(-xb),plogis(xb))
        efx<-dlogis(xb)*b
        }
        else stop("Error in MLefx: must provide either atx or aty values")
    efx<-cbind(b,efx,b*.25)
    dimnames(efx)<-list(object$varnames,c("Coef","  d(Y=1)/dx","  max(dy/dx)"))    
    cat("\n Marginal Effects for Binary Logit\n")
    if (!is.null(aty)) cat("\n Evaluated at P(y) =",round(pj[2],3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X = ",round(atx,3),"\n")
    if (!is.null(atx)) cat(" Predicted P(Y=1)=",round(pj[2],3),"\n")
    print(round(efx,4))
    }


MLefx.loglog<-function(object,atx=NULL,aty=NULL){
    print("Sorry, I don't think I know how to do that yet!")
    }



MLefx.cloglog<-function(object,atx=NULL,aty=NULL){
    print("Sorry, I don't think I know how to do that yet!")
    }



MLefx.cauchy<-function(object,atx=NULL,aty=NULL){
    print("Sorry, I don't think I know how to do that yet!")
    }


MLefx.poisson<-function(object,atx=NULL,aty=NULL){
    nx<-object$nx
    b<-object$par
    
    if (!is.null(aty)){
        efx<-mean(aty)*b
        }
    else
    if (!is.null(atx)){
        atx<-c(1,atx)
        xb<-as.vector(atx%*%b)
        efx<-exp(xb)*b
        }
        else stop("Error in MLefx: must provide either atx or aty values")
    efx<-cbind(b,efx)
    dimnames(efx)<-list(object$varnames,c("Coef","  dE(y)/dx"))    
    cat("\n Marginal Effects for Poisson\n")
    if (!is.null(aty)) cat("\n Evaluated at mean(aty) =",round(mean(aty),3),"\n")
    if (!is.null(atx)) cat("\n Evaluated at X = ",round(atx,3),"\n")
    if (!is.null(atx)) cat("\n Estimated lambda = ",round(exp(xb),3),"\n")
    print(round(efx,4))
    }
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
# Logit Regression via ML in R
# Version 3. 
# copyright Charles H. Franklin, Tuesday, July 4, 2006 at 11:55

  MLlogit<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
 
    llf0<-(log(mean(y))*sum(y) + log(1-mean(y))*sum(1-y))
    neglnl<-function(theta,X,y){
            p<-as.vector(1/(1+exp(-X%*%theta)))
            lnl<-(y*log(p) + (1-y)*log(1-p))
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
        Lambdax<-as.vector(1/(1+exp(-X%*%theta)))
        part1<-as.vector(y-Lambdax)   
        -apply(part1*X,2,sum)
   }
    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    se0<-sqrt(sum((y-X%*%b0)^2)/(length(y)-nx))
    theta<-b0/se0*1.6
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),nx=nx,llf0=llf0))
    result$fitted.values<-1/(1+exp(-X%*%result$par))
    result$xb<-X%*%result$par
    class(result)<-"MLlogit"        
    return(result)
    }
    
## Print method for MLlogit

    print.MLlogit <- function(object){
        coef<-object$par
        print(coef)
        if (object$convergence==0) cat('\n MLlogit converged\n')
        if (!object$convergence==0) cat('\n *** MLlogit failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLlogit

    summary.MLlogit<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Logit\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*(object$llf0--object$value)
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
# Log-log Regression via ML in R
# Version 3. 
# copyright Charles H. Franklin, Tuesday, July 4, 2006 at 12:08

  MLloglog<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 

    llf0<-(log(mean(y))*sum(y) + log(1-mean(y))*sum(1-y))
    neglnl<-function(theta,X,y){
            p<-as.vector(exp(-exp(-X%*%theta)))
            lnl<-y*log(p) + (1-y)*log(1-p)
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
        ememxb<-exp(-exp(-X%*%theta))       
        part1<-as.vector((y-(1-y)*ememxb/(1-ememxb))*exp(-X%*%theta))   
        -apply(part1*X,2,sum)
   }
    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    se0<-sqrt(sum((y-X%*%b0)^2)/(length(y)-nx))
    theta<-b0/se0
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),nx=nx,llf0=llf0))
    result$fitted.values<-as.vector(exp(-exp(-X%*%result$par)))
    result$xb<-X%*%result$par
    class(result)<-"MLloglog"        
    return(result)
    }
    
## Print method for MLloglog

    print.MLloglog <- function(object){
        coef<-object$par
        print(coef)
        if (object$convergence==0) cat('\n MLloglog converged\n')
        if (!object$convergence==0) cat('\n *** MLloglog failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLloglog

    summary.MLloglog<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Log-log model\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*(object$llf0--object$value)
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
## MLlrtest --- Simple function to do LR Tests for MLobjects
# by chf Sunday, July 18, 2004 at 12:26

    MLlrtest<-function(object1,object2){
        lnl1<-object1$value
        lnl2<-object2$value
        df<-abs(length(object1$par)-length(object2$par))
        lr<-2*abs(lnl1-lnl2)
        plr<-1-pchisq(lr,df)
        
        cat("\nML Likelihood Ratio Test")
        cat("\nChi-square:", round(lr,3), "df = ",df)
        cat("\nProb x>chisq:", round(plr,5),"\n")
        }
#MLmnl --- Multinomial Logit via ML in R
#Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:16


  MLmnl<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 


    J<-max(y)
    
    neglnl<-function(theta,X,y){
            options(warn=-1)
            nx<-ncol(X)
            J<-max(y)
            b<-matrix(theta, nx, J-1)
            xb<-X%*%b   #N x J
            denom<-1+rowSums(exp(xb))
            p<-rep(0,length(y))
            for (j in 1:J){
                if (j==1) {
                    p[y==1]<-1/denom[y==1]
                    }
                else 
                    p[y==j]<-exp(xb[y==j,j-1])/denom[y==j]
            }
            lnl<-log(p)
            -sum(lnl)
        }

    theta<-numeric(0)
    for (j in 2:J){
        tx<-X[y==1|y==j,]
        ty<-as.numeric(y[y==1|y==j]==j)
        b0<-solve(crossprod(tx))%*%crossprod(tx,ty)
        theta<-c(theta,b0)
        }
    
    llf0<-sum(table(y)*log(table(y)/length(y)))

    result<-c(optim(theta, neglnl, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),
                nx=nx,J=J,llf0=llf0))
    class(result)<-"MLmnl"        
    return(result)
    }
    
## Print method for MLmnl

    print.MLmnl <- function(object){
        J<-object$J
        coef<-object$par
        names(coef)<-rep(object$varnames,J-1)
        print(coef)
        if (object$convergence==0) cat('\n MLmnl converged\n')
        if (!object$convergence==0) cat('\n *** MLmnl failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLmnl

    summary.MLmnl<-function(object, covar=FALSE){
        J<-object$J
        nx<-object$nx
        coef<-object$par
        names(coef)<-rep(object$varnames,J-1)
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Multinomial Logit\n")
        cat("\n  Estimated Parameters\n")
        for (j in 2:J){
            cat("\n For outcome Y = ",j,"\n")
            print(coef.table[((j-2)*nx+1):((j-1)*nx),])
        }
        
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
#MLnegbin --- Negative Binomial Event Counts via ML in R
# Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:18
# See Cameron & trivedi p. 71

  MLnegbin<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 

    neglnl<-function(theta,X,y){
            options(warn=-1)
            b<-theta[1:(length(theta)-1)]
            a<-theta[length(theta)]
            xb<-as.vector(X%*%b)
            exb<-exp(xb)
            suma=rep(0,length(y))
            for (i in 1:length(y)){
                for (j in 0:(y[i]-1)) suma[i]<-suma[i]+log(j+(1/a))
                }            
            lnl<-suma-lgamma(y+1)-(y+(1/a))*log(1+a*exb)+y*log(a)+y*xb
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
            b<-theta[1:(length(theta)-1)]
            a<-theta[length(theta)]
            xb<-as.vector(X%*%b)
            exb<-exp(xb)
            sumb=rep(0,length(y))
            for (i in 1:length(y)){
                for (j in 0:(y[i]-1)) sumb[i]<-sumb[i]+1/(j+(1/a))
                }
            part1<-((y-exb)/(1+a*exb))*X
            part2<-(1/(a^2))*(log(1+a*exb)-sumb)+(y-exb)/(a*(1+a*exb))  
            gr<-cbind(part1,part2)
            -apply(gr,2,sum)
            }
    
# Fit poisson for both LR test of alpha and to get starts

    poisfit<-MLpoisson(y,X)
    ll
    
    theta<-c(solve(crossprod(X))%*%crossprod(X,y),1/sd(y))
    llf0<-sum(y*log(sum(y)/length(y))-sum(y)/length(y))
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=c(colnames(X),"alpha"),
                nx=nx,llf0=llf0))
    class(result)<-"MLnegbin"        
    return(result)
    }
    
## Print method for MLnegbin

    print.MLnegbin <- function(object){
        coef<-as.vector(object$par)
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n MLnegbin converged\n')
        if (!object$convergence==0) cat('\n *** MLnegbin failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLnegbin

    summary.MLnegbin<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Negative Binomial Count Model\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
#MLologit --- Ordered Logit via ML in R
#Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:20


  MLologit<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
    
    J<-max(y)
    neglnl<-function(theta,X,y){
            options(warn=-1)
            b<-theta[1:nx]
            tau<-c(0,theta[(nx+1):length(theta)]) #tau1=0
            J<-max(y)
            xb<-X%*%b
            p<-rep(0,length(y))
            for (j in 1:J){
                if (j==1) {
                    p[y==1]<-plogis(-xb[y==1])
                    }
                else if (j==J){
                    p[y==J]<-1-plogis(tau[J-1]-xb[y==J])
                    }
                else p[y==j]<-plogis(tau[j]-xb[y==j])-plogis(tau[j-1]-xb[y==j])
            }
            lnl<-log(p)
            -sum(lnl)
        }

    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    theta<-c(b0,1:(J-2))
    llf0<-sum(table(y)*log(table(y)/length(y)))

    result<-c(optim(theta, neglnl, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=c(colnames(X),paste("tau",seq(1:(J-2)),sep=" ")),
                nx=nx,J=J,llf0=llf0))
    class(result)<-"MLologit"        
    return(result)
    }
    
## Print method for MLologit

    print.MLologit <- function(object){
        J<-object$J
        coef<-object$par
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n MLologit converged\n')
        if (!object$convergence==0) cat('\n *** MLologit failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLologit

    summary.MLologit<-function(object, covar=FALSE){
        J<-object$J
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Ordered Logit\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
#MLoprobit --- Ordered Probit via ML in R
#Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:22

  MLoprobit<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 

    
    J<-max(y)
    neglnl<-function(theta,X,y){
            options(warn=-1)
            b<-theta[1:nx]
            tau<-c(0,theta[(nx+1):length(theta)]) #tau1=0
            J<-max(y)
            xb<-X%*%b
            p<-rep(0,length(y))
            for (j in 1:J){
                if (j==1) {
                    p[y==1]<-pnorm(-xb[y==1])
                    }
                else if (j==J){
                    p[y==J]<-1-pnorm(tau[J-1]-xb[y==J])
                    }
                else p[y==j]<-pnorm(tau[j]-xb[y==j])-pnorm(tau[j-1]-xb[y==j])
            }
            lnl<-log(p)
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
            nx<-ncol(X)
            b<-theta[1:nx]
            tau<-c(0,theta[(nx+1):length(theta)]) #tau1=0
            ntau<-length(theta)-nx
            J<-max(y)
            xb<-X%*%b
            gr<-matrix(0,length(y),length(theta))
            for (j in 1:J){
                if (j==1) {
                    gr[y==1]<-dnorm(-xb[y==1])*c((-X[y==1]),rep(0,ntau))/pnorm(-xb[y==1])
                    }
                else if (j==J){
                    gr[y==J]<-dnorm(tau[J-1]-xb[y==J])*c(-X[y==J],rep(0,ntau-1),1)/(1-pnorm(tau[J-1]-xb[y==J]))
                    }
                else {
                    dtj<-rep(0,ntau)
                    dtjm1<-rep(0,ntau)
                    dtj[j-1]<-1
                    if (j>2) dtjm1[j-2]<-1
                    gr[y==j]<-(dnorm(tau[j]-xb[y==j])*c(-X[y==j],dtj) - 
                                dnorm(tau[j-1]-xb[y==j])*c(-X[y==j],dtjm1))/
                                (pnorm(tau[j]-xb[y==j])-pnorm(tau[j-1]-xb[y==j]))
                    }
                }
#                print(-apply(gr,2,sum))
                -apply(gr,2,sum)
            }

    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    theta<-c(b0,1:(J-2))
    llf0<-sum(table(y)*log(table(y)/length(y)))
    result<-c(optim(theta, neglnl, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=c(colnames(X),paste("tau",seq(1:(J-2)),sep=" ")),
                nx=nx,J=J,llf0=llf0))
    class(result)<-"MLoprobit"        
    return(result)
    }
    
## Print method for MLoprobit

    print.MLoprobit <- function(object){
        J<-object$J
        coef<-object$par
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n MLoprobit converged\n')
        if (!object$convergence==0) cat('\n *** MLoprobit failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLoprobit

    summary.MLoprobit<-function(object, covar=FALSE){
        J<-object$J
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Ordered Probit\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
#MLpoisson --- Poisson Event Counts via ML in R
#Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:24

  MLpoisson<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 

    neglnl<-function(theta,X,y){
            options(warn=-1)
            xb<-X%*%theta
            lnl<-y*xb-exp(xb)
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
            xb<-as.vector(X%*%theta)
            gr<-y*X - exp(xb)*X
            -apply(gr,2,sum)
            }
    
    theta<-solve(crossprod(X))%*%crossprod(X,y)
    llf0<-sum(y*log(sum(y)/length(y))-sum(y)/length(y))
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X)),
                nx=nx,llf0=llf0)
    class(result)<-"MLpoisson"        
    return(result)
    }
    
## Print method for MLpoisson

    print.MLpoisson <- function(object){
        coef<-as.vector(object$par)
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n MLpoisson converged\n')
        if (!object$convergence==0) cat('\n *** MLpoisson failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLpoisson

    summary.MLpoisson<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Poisson Count Model\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
# Probit Regression via ML in R
# Version 3. 
# copyright Charles H. Franklin, Tuesday, July 4, 2006 at 14:05
MLprobit<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
    
    llf0<-(log(mean(y))*sum(y) + log(1-mean(y))*sum(1-y))
    
    neglnl<-function(theta,X,y){
            p<-as.vector(pnorm(X%*%theta))
            lnl<-(y*log(p) + (1-y)*log(1-p))
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
        phix<-as.vector(dnorm(X%*%theta))
        Phix<-as.vector(pnorm(X%*%theta))
        part1<-as.vector(y*phix/Phix-(1-y)*phix/(1-Phix))   
        -apply(part1*X,2,sum)
   }
    
    b0<-solve(crossprod(X))%*%crossprod(X,y)
    se0<-sqrt(sum((y-X%*%b0)^2)/(length(y)-nx))
    theta<-b0/se0    
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X),nx=nx,llf0=llf0))
    result$fitted.values<-pnorm(X%*%result$par)
    result$xb<-X%*%result$par
    result$call<-call
    class(result)<-"MLprobit"        
    return(result)
    }
    
## Print method for MLprobit

    print.MLprobit <- function(object){
        coef<-object$par
        print(coef)
        if (object$convergence==0) cat('\n MLprobit converged\n')
        if (!object$convergence==0) cat('\n *** MLprobit failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLProbit

    summary.MLprobit<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Probit\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*(object$llf0--object$value)
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }
##MLwald--- Function to do Wald test on ML objects
## by chf Sunday, July 18, 2004 at 07:17
#Usage:
#    R<-matrix(c(0,1,0,0,0,
#                0,0,1,0,0,
#                0,0,0,1,0,
#                0,0,0,0,1),4,5,byrow=T)
#    q<-c(0,0,0,0)
#    wald(MLobject,R,q)

    MLwald<-function(object,R,q) {
        if (!is.matrix(R)) stop("Restrictions must be a matrix")
        b<-object$par
        vc<-solve(object$hessian)
        w<-t(R%*%b-q)%*%solve(R%*%vc%*%t(R))%*%(R%*%b-q)
        pw<-1-pchisq(w,length(q))
        cat("\nWald Test\n")
        cat("\nR = \n")
        print(R)
        cat("\nq = ",q)
        cat("\nChi-square:", round(w,3),"  df = ",length(q))
        cat("\nProb x>chisq:", round(pw,5),"\n")
       }
       
#MLzip --- Zero Inflated Poisson Event Counts via ML in R
#Version 3
# copyright Charles H. Franklin Tuesday, July 4, 2006 at 13:26

  MLzip<-function (formula, data, subset, weights, na.action, start = NULL, offset, 
    model = TRUE, method = "BFGS", x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "start", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    y <- model.response(mf, "any")
    
    if (length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if (!is.null(nm)) 
            names(y) <- nm
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    start <- model.extract(mf, "start")
    nx<-ncol(X)

#Likelihood code etc follows 
    
    neglnl<-function(theta,X,y){
            options(warn=-1)
            xb<-X%*%theta
            lnl<-y*xb-exp(xb)
            -sum(lnl)
        }
        
    grad<-function(theta,X,y){
            xb<-as.vector(X%*%theta)
            gr<-y*X - exp(xb)*X
            -apply(gr,2,sum)
            }
    
    theta<-solve(crossprod(X))%*%crossprod(X,y)
    llf0<-sum(y*log(sum(y)/length(y))-sum(y)/length(y))
    result<-c(optim(theta, neglnl, gr=grad, hessian=T, 
                method=method, control=(maxit=1000), X=X, y=y),
                list(varnames=colnames(X)),
                nx=nx,llf0=llf0)
    class(result)<-"MLpoisson"        
    return(result)
    }
    
## Print method for MLpoisson

    print.MLpoisson <- function(object){
        coef<-as.vector(object$par)
        names(coef)<-object$varnames
        print(coef)
        if (object$convergence==0) cat('\n MLpoisson converged\n')
        if (!object$convergence==0) cat('\n *** MLpoisson failed to converge ***\n')       
        invisible(object)
        }


##Summary function for MLpoisson

    summary.MLpoisson<-function(object, covar=FALSE){
        coef<-object$par
        names(coef)<-object$varnames
        nx<-object$nx
        maxl<-object$value
        vc<-solve(object$hessian)
        colnames(vc)<-names(coef)
        rownames(vc)<-names(coef)
        se<-sqrt(diag(vc))
        zscore<-coef/se
        pz<-2*pnorm(-abs(coef/se))
        niter<-object$counts
        converge<-object$convergence
        dn<-c("Estimate","Std. Error")
        coef.table<-cbind(coef,se,zscore,pz)
        dimnames(coef.table)<-list(names(coef),c(dn,"z-value",
            "Pr(>|z|)"))
        cat("\n Maximum Likelihood Poisson Count Model\n")
        cat("\n  Estimated Parameters\n")
        print(coef.table)
        cat("\nLog-Likelihood: ",-object$value,"\n")
            
        lrtest<--2*((object$llf0--object$value))
            plr<-1-pchisq(lrtest,length(coef)-1)
            cat("\nLikelihood Ratio Test")
            cat("\nChi-square: ", lrtest)
            cat("\nProb x>chisq: ", plr,"\n")
            
        cat("\nNumber of calls to likelihood and gradient: ",niter)
        if (converge==0) {cat("\nEstimation Converged\n")}
            else {cat("\nEstimation DID NOT Converge. Condition code: ",converge,"\n")}
            
        if (covar) {
            cat("\nVariance-Covariance Matrix for Parameters\n")
            print(vc) }
        }