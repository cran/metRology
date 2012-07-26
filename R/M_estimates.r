MM.estimate<-function(x, ...) {
        UseMethod("MM.estimate")
}

MM.estimate.default<-function(x, u, c=4.685, ...) {
        require(MASS)
        rl<-rlm(x~1, weights=1/u^2, c=c, method="MM")
        srl<-summary(rl)
        rv<-.construct.loc.est(x=coef(srl)[1], u=coef(srl)[2], xi=x, ui=u, u.eff=u*rl$s, 
                w=rl$w, method="MM", method.details=rl)
        return(rv)
}


huber.estimate<-function(x, ...) {
        UseMethod("huber.estimate")
}

huber.estimate.default<-function(x, u, k= 1.345, ...) {
        require(MASS)
        rl<-rlm(x~1, weights=1/u^2, k=k, method="M")
        srl<-summary(rl)
        rv<-.construct.loc.est(x=coef(srl)[1], u=coef(srl)[2], xi=x, ui=u, u.eff=u*rl$s, 
                w=rl$w, method="Huber", method.details=rl)
        return(rv)
}

