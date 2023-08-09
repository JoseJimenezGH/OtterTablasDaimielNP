#==============================================================================#
#                                                                              #
#                               FUNCIONES SCR                                  #
#                                Jose Jimenez                                  #
#                             07/05/2020 8:20:49                               #
#                                                                              #
#==============================================================================#

# BONDAD DEL AJUSTE
# ====================
SCRgof<-function (out, nx = 6, ny = 6, traplocs = NULL, buffer = 2, Xl = NULL,
    Xu = NULL, Yl = NULL, Yu = NULL)
{
    S <- out$s
    Sxout <- S[, , 1]
    Syout <- S[, , 2]
    z <- out$w
    Xl <- min(traplocs[, 1]) - buffer
    Xu <- max(traplocs[, 1]) + buffer
    Yl <- min(traplocs[, 2]) - buffer
    Yu <- max(traplocs[, 2]) + buffer
    niter <- nrow(z)
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    area <- (Xu - Xl) * (Yu - Yl)
    Sxout2 <- cut(Sxout[z == 1], breaks = xg)
    Syout2 <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout2, Syout2)/niter
    op <- par(mar = c(4, 4, 4, 6))
    on.exit(par(op))
    image(xg, yg, Dn, col = terrain.colors(10))
    image.scale(Dn, col = terrain.colors(10))
    title("Tamano de la poblacion local \n(individuos por unidad de superficie)")
    stat2 <- statsim2 <- stat <- statsim <- rep(NA, niter)
    for (i in 1:niter) {
        inside <- (Sxout[i, ] < Xu) & (Sxout[i, ] > Xl) & (Syout[i,
            ] < Yu) & (Syout[i, ] > Yl)
        Dn <- table(cut(Sxout[i, ][z[i, ] == 1 & inside], breaks = xg),
            cut(Syout[i, ][z[i, ] == 1 & inside], breaks = yg))
        Dnv <- Dn[1:length(Dn)]
        E <- mean(Dnv)
        stat[i] <- (var(Dnv)/mean(Dnv))
        stat2[i] <- sum((sqrt(Dnv) - sqrt(E))^2)
        Sxsim <- runif(sum(z[i, ][inside]), Xl, Xu)
        Sysim <- runif(sum(z[i, ][inside]), Yl, Yu)
        Dnsim <- table(cut(Sxsim, breaks = xg), cut(Sysim, breaks = yg))
        Dnsimv <- Dnsim[1:length(Dnsim)]
        statsim[i] <- (var(Dnsimv)/mean(Dnsimv))
        statsim2[i] <- sum((sqrt(Dnsimv) - sqrt(E))^2)
    }
    out <- cbind(data = stat2, newdata = statsim2)
    cat("Indice de dispersion observado: ", mean(stat), fill = TRUE)
    cat("Indice de dispersion simulado: ", mean(statsim), fill = TRUE)
    cat("P-valor del Indice of dispersion: ", mean(statsim > stat),
        fill = TRUE)
    cat("P-valor(2) freeman-tukey: ", mean(statsim2 > stat2), fill = TRUE)
    invisible(out)
}

# AJUSTE GAUSSIANO PARA HRA
#=============================
pGauss1<-function (parms, Dmat) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    p <- plogis(parms[1]) * exp(-(1/(2 * parms[2] * parms[2])) * 
        Dmat * Dmat)
    p
}


# DISTANCIA ENTRE PUNTOS DE MUESTREO Y CENTROS DE ACTIVIDAD
#=============================================================
e2dist1 <- function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


# SIMULACION DE DATOS
#=====================
sim.data <- function(N=30, K=5, sigma=.6, lam0=1.5, side=10,
                     xlims=c(0, 15), ylims=c(0, 15))
{
    coords <- seq(3, 12, length=side)
    X <- cbind(x=rep(coords, each=side), y=rep(coords, times=side))
    J <- nrow(X) # Numero de localizaciones de trampas
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy) # Centros de actividad
    D <- e2dist1(S, X) # Matriz de distancias
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma)) # Ratio de encuentros esperado
    z <- array(NA, c(N, J, K))
    n <- matrix(NA, nrow=J, ncol=K)
    for(k in 1:K) {
        z[,,k] <- rpois(N*J, lam) # encuentros latentes
        n[,k] <- colSums(z[,,k])  # conteos observados
    }
    list(n=n, N=N, sigma=sigma, lam0=lam0, X=X, S=S, z=z, xlims=xlims,
         ylims=ylims)
}



# MUESTREO INCONDICIONAL UTILIZANDO MCMC
#========================================
# Empleamos este metodo cuando no disponemos de ningun animal marcado
# Hacemos dos muestreos de Gibbs para ver como se superponen las cadenas
# Esta version no actualiza z
spNmix <- function(n, X, M, niters, xlims, ylims, tune=c(0.1, 0.1, 2),
                   monitorS=FALSE)
{
    K <- ncol(n)
    # valores iniciales
    S <- cbind(runif(M, xlims[1], xlims[2]),
               runif(M, ylims[1], ylims[2]))
    D <- e2dist1(S, X)
    sigma <-runif(1, .4, 1)
    lam0 <- runif(1, .1, 1)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    w <- rbinom(M, 1, .7)
    psi <- runif(1, .2, .8)
    lamv.curr <- colSums(lam*w)
    # para el caso de que el primer sigma sea rechazado
    ll <- sum(dpois(n, lamv.curr, log=TRUE))

    # matriz para almacenar los muestreos
    out <- matrix(NA, nrow=niters, ncol=4)
    colnames(out) <- c("sigma","lam0","psi","N")

    Sout <- NULL
    if(monitorS)
        Sout <- array(NA, c(M, 2, niters))

    cat("\nValores de inicio =", c(sigma, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 ==0) {
            cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("actual =", out[iter-1,], "\n")
            cat("  Ratios de aceptacion\n")
            cat("    S =", Sups/M, "\n")
            cat("    w =", wUps/M, "\n")
        }

        # actualiza sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            lamv.cand <- colSums(lam.cand*w)
            lamv.curr <- colSums(lam*w)
            ll<- sum(dpois(n, lamv.curr, log=TRUE))
            llcand<- sum(dpois(n, lamv.cand, log=TRUE) )
            if(runif(1) < exp( llcand  - ll ) ){
                ll <- llcand
                lamv.curr <- lamv.cand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }

        # actualiza lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            lamv.cand <- colSums(lam.cand*w)
            llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
            if(runif(1) < exp( llcand - ll ) ) {
                ll <- llcand
                lamv.curr<-lamv.cand
                lam0<-lam0.cand
                lam<-lam.cand
            }
        }


        # actualiza w
        wUps <- 0
        for(i in 1:M) {
            wcand <- w
            wcand[i] <- if(w[i]==0) 1 else 0
            lamv.cand <- colSums(lam*wcand)
            llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp((llcand+prior.cand) - (ll+prior))) {
                w <- wcand
                lamv.curr <- lamv.cand
                ll <- llcand
                wUps <- wUps+1
            }
        }

        # actualiza psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # actualiza S
        Sups <- 0
        for(i in 1:M) {
            Scand <-c(rnorm(1, S[i,1], tune[3]), rnorm(1, S[i,2], tune[3]))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
                     Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
            lam.cand <- lam
            lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
            lamv.cand <- colSums(lam.cand*w)
            llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
            if(runif(1)< exp(llcand - ll)) {
                ll <- llcand
                lamv.curr <- lamv.cand
                S[i,] <- Scand
                lam <- lam.cand
                D[i,] <- dtmp
                Sups <- Sups+1
            }
        }
        out[iter,] <- c(sigma,lam0,psi,sum(w))
        if(monitorS)
            Sout[1:sum(w),,iter] <- S[w==1,]
    }
    last <- list(S=S, lam=lam, w=w)
    list(out=out, last=last, Sout=Sout)
}




# MUESTREO CONDICIONAL UTILIZANDO MCMC
#======================================
# Esta versiun actualiza z. Condicional. Es util cuando uno o mas
# individuos esta disponible para ser identificado individualmente, en cuyo
# caso la frecuencia de encuentros zjk es observable para estos individuos
# Le he anadido un monitorizado de S.

funcZ <- function(y, X, Zknown, M, niters, xlims, ylims,
  tune=c(0.1, 0.01, 2), monitorS=FALSE) {

    R <- nrow(y)
    T <- ncol(y)
    S <- cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1], ylims[2]))
    D <- e2dist1(S, X)
    theta <- runif(1, 3, 7)
    lam0 <- runif(1, .1, 0.5)
    lam <- lam0*exp(-(D*D)/(2*theta*theta))
    w <- rbinom(M,1,.5)
    psi <- runif(1,.2,.8)

    Z <- array(NA, c(M,R,T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if(!missing(Zknown) && !is.null(Zknown)) {
        nMarked <- nrow(Zknown)
        marked[1:nMarked] <- TRUE
        Z[1:nMarked,,] <- Zknown
    }
    w[marked] <- 1
    Zdata <- !is.na(Z)
    for(r in 1:R) {
        for(t in 1:T) {
            if(y[r,t]==0) {
                Z[,r,t] <- 0
                next
            }
            unmarked <- !Zdata[,r,t]
            nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
            probs <- lam[,r]*w
            probs <- probs[unmarked]
            Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs)
        }
    }

    out <- matrix(NA,nrow=niters,ncol=5)
    colnames(out) <- c("sigma", "lam0", "psi", "N", "DS")

    Sout <- NULL
    if(monitorS)
        Sout <- array(NA, c(M, 2, niters))

    cat("\nValores de inicio =", c(theta, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   actual =", out[iter-1,], "\n")
        }

        # Actualiza theta
        theta.cand <- rnorm(1, theta, tune[1])
        ll <- sum(dpois(Z, lam*w, log=TRUE)) # This needs to be here
        if(theta.cand>0){
            lam.cand <- lam0*exp(-(D*D)/(2*theta.cand*theta.cand))
            # w is recycled over lam, R times
            # lam*w is recycled over Z, T times
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp( llcand  - ll ) ){
                ll<-llcand
                lam<-lam.cand
                theta<-theta.cand
            }
        }

        # Actualiza lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*theta*theta))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp(llcand - ll)) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

        ### Actualiza "w" aqui
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i] | marked[i])
                next
            wcand <- w

            if(w[i]==0) {
                wcand[i] <- 1
                ll <- 0
                llcand <- sum(dpois(Z[i,,], lam[i,]*wcand[i], log=TRUE))
            } else {
                wcand[i] <- 0
                ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- 0
            }

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand[i]
                wUps <- wUps+1
            }
        }

        # Actualiza Z
        for(r in 1:R) {
            probs <- lam[,r]*w # These get normalized by rmultinom
            for(t in 1:T) {
                if(y[r,t]==0) {
                    Z[,r,t] <- 0
                    next
                }
                unmarked <- !Zdata[,r,t]
                nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown == 0)
                    next
                Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs[unmarked])
            }
        }

        # Actualiza psi
#        psi<-rbeta(1,1+sum(w),1+M-sum(w))
        psi <- rbeta(1,1+sum(w[!marked]), 1+sum(!marked)-sum(w[!marked]))

        # actualiza Ds
        DS<- sum(w)/(100*(xlims[2]-xlims[1])*(ylims[2]-ylims[1]))

        # Actualiza S
        Sups <- 0
        for(i in 1:M) {
            Scand <- c(rnorm(1, S[i,1], tune[3]),
                       rnorm(1, S[i,2], tune[3]))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
                     Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(inbox) {
                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
                lam.cand <- lam
                lam.cand[i,] <- lam0*exp(-(dtmp*dtmp)/(2*theta*theta) )

                if(w[i]==0) {
                    ll <- llcand <- 0
                } else {
                    ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                    llcand <- sum(dpois(Z[i,,], lam.cand[i,]*w[i],
                                        log=TRUE))
                }

                if(runif(1) < exp(llcand - ll)) {
                    S[i,] <- Scand
                    lam <- lam.cand
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
        }

        if(iter %% 100 == 0) {
            cat("   Ratios de aceptacion\n")
            cat("     w =", wUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(theta,lam0,psi,sum(w),DS)
        if(monitorS)
            Sout[1:sum(w),,iter] <- S[w==1,]
    }

    list(out=out, Sout=Sout)
}


e2dist<-function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


image.scale<-function (z, col, x, y = NULL, size = NULL, digits = 2,
  labels = c("breaks", "ranges"))
{
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2])
    my <- mean(usr[3:4])
    dx <- diff(usr[1:2])
    dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05 * dx/2
    else if (is.list(x)) {
        if (length(x$x) == 2)
            size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    }
    else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
            size <- 0.618 * dy/n
            y <- my + 0.618 * dy/2
        }
        else size <- (y - my) * 2/n
    if (length(size) == 1)
        size <- rep(size, 2)
    if (is.null(y))
        y <- my + n * size[2]/2
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
        col = rev(col), xpd = TRUE)
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format = "f", digits = digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n + 1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj = ifelse(size[1] >
        0, 0, 1), xpd = TRUE)
}


SCRdensity<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col = terrain.colors(ncolors))
    image.scale(Dn, col = terrain.colors(ncolors))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


# Con este comando sacamos un distribución real del espacio entre el eje de 
# las X y el de las Y. Para que sea regular el pixel, lo ajusto con asratio
# plot(X, asp=1)
# asratio <- (Yu-Yl)/(Xu-Xl)
# nx<-75
# ny<-75* asratio

SCRdensityJJ1<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=gray.colors(ncolors, start=1, end=0),asp=1)
    image.scale(Dn, col=gray.colors(ncolors, start=1, end=0))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


SCRdensityJJ2 <- function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=gray.colors(ncolors, start=1, end=0),asp=1)
    image.scale(Dn, col=gray.colors(ncolors, start=1, end=0))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


SCRdensityJJ3 <- function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=terrain.colors(10),asp=1)
    image.scale(Dn, col=terrain.colors(10))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


SimSCR0<-function (N = 100, K = 20, p0=0.1, sigma=0.5, discard0 = TRUE,
  array3d = FALSE, rnd = 2013){
    set.seed(rnd)
    traplocs <- cbind(sort(rep(1:10, 10)), rep(1:10, 10))
    Dmat <- e2dist(traplocs, traplocs)
    ntraps <- nrow(traplocs)
    delta <- 2.5*sigma
    Xl <- min(traplocs[, 1] - delta)
    Xu <- max(traplocs[, 1] + delta)
    Yl <- min(traplocs[, 2] - delta)
    Yu <- max(traplocs[, 2] + delta)
    sx <- runif(N, Xl, Xu)
    sy <- runif(N, Yl, Yu)
    S <- cbind(sx, sy)
    D <- e2dist(S, traplocs)
    alpha0 <- logit(p0)
    sigma <- 0.5
    alpha1 <- 1/(2 * sigma * sigma)
    probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
    Y <- matrix(NA, nrow = N, ncol = ntraps)
    for (i in 1:nrow(Y)) {
        Y[i, ] <- rbinom(ntraps, K, probcap[i, ])
    }
    if (discard0) {
        totalcaps <- apply(Y, 1, sum)
        Y <- Y[totalcaps > 0, ]
    }
    dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), sep = ""))
    if (array3d) {
        Y <- array(NA, dim = c(N, ntraps, K))
        for (i in 1:nrow(Y)) {
            for (j in 1:ntraps) {
                Y[i, j, 1:K] <- rbinom(K, 1, probcap[i, j])
            }
        }
        if (discard0) {
            Y2d <- apply(Y, c(1, 2), sum)
            ncaps <- apply(Y2d, 1, sum)
            Y <- Y[ncaps > 0, , ]
        }
    }

    list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl,
        Yu), N = N, p0=p0, alpha1 = alpha1, sigma = sigma,
        K = K, S=S)
}


simSCR0Poiss<-function (discard0 = TRUE, N=N, K=25,  p0=0.1,  sigma=0.5,
  buffer=buffer, array3d = FALSE, rnd = 2013)
{
    set.seed(rnd)
    traplocs <- cbind(sort(rep(1:15, 15)), rep(1:15, 15))
    Dmat <- e2dist(traplocs, traplocs)
    ntraps <- nrow(traplocs)
    plot(traplocs)
    buffer <- 2.5*sigma
    Xl <- min(traplocs[, 1] - buffer)
    Xu <- max(traplocs[, 1] + buffer)
    Yl <- min(traplocs[, 2] - buffer)
    Yu <- max(traplocs[, 2] + buffer)
    sx <- runif(N, Xl, Xu)
    sy <- runif(N, Yl, Yu)
    S <- cbind(sx, sy)
    plot(traplocs, pch=3, col='blue', xlim=c(Xl, Xu), ylim=c(Yl, Yu))
    points(S, pch=16, col='red')
    D <- e2dist(S, traplocs)
    alpha1 <- 1/(2 * sigma * sigma)
    alpha0<-log(p0)
    muy <- exp(alpha0) * exp(-alpha1 * D * D)
    Y <- matrix(NA, nrow = N, ncol = ntraps)
    for (i in 1:nrow(Y)) {
      Y[i, ] <- rpois(ntraps, K * muy[i, ])
    }
    if (discard0) {
      totalcaps <- apply(Y, 1, sum)
      Y <- Y[totalcaps > 0, ]
    }
    dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), sep = ""))
    if (array3d) {
      Y <- array(NA, dim = c(N, K, ntraps))
        for (i in 1:nrow(Y)) {
          for (j in 1:ntraps) {
            Y[i, 1:K, j] <- rpois(K, muy[i, j])
          }
        }
      if (discard0) {
        Y2d <- apply(Y, c(1, 3), sum)
        ncaps <- apply(Y2d, 1, sum)
        Y <- Y[ncaps > 0, , ]
      }
    }
    Y<-aperm(Y,c(1,3,2))
    list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl,Yu),
      N = N, p0 = p0, sigma = sigma,  S=S, K = K, buffer=buffer)
}


sim.pID.Data <- function (N = N, K = K, sigma = sigma, lam0 = lam0,
    knownID = knownID, X = X, xlims = xlims, ylims = ylims,
    obsmod = c("pois", "bern"),nmarked = c("known", "unknown"), rat = 1,
    tel = 0, nlocs = 0)
{
    if (tel > knownID)
        stop("tel cannot be bigger than knownID")
    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)
    npts <- dim(X)[1]
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    D <- e2dist(S, X)
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    Y <- array(NA, c(N, npts, K))
    for (i in 1:N) {
        for (j in 1:npts) {
            if (identical(obsmod, "bern")) {
                Y[i, j, ] <- rbinom(K, 1, lam[i, j])
            }
            else if (identical(obsmod, "pois")) {
                Y[i, j, ] <- rpois(K, lam[i, j])
            }
        }
    }
    n <- apply(Y, c(2, 3), sum)
    Yknown <- Y[1:knownID, , ]
    if (identical(nmarked, "unknown")) {
        iobs <- which(apply(Yknown > 0, 1, any))
        Yobs <- Yknown[iobs, , ]
    }
    else if (identical(nmarked, "known")) {
        Yobs <- Yknown
    }
    YknownR <- Yobs
    counter <- array(0, c(dim(Yobs)[1], dim(X)[1], K))
    for (i in 1:dim(Yobs)[1]) {
        for (j in 1:dim(X)[1]) {
            for (k in 1:K) {
                if (identical(obsmod, "bern")) {
                  if (YknownR[i, j, k] == 1) {
                    IDed <- rbinom(1, 1, rat)
                    if (IDed == 0) {
                      YknownR[i, j, k] <- 0
                      counter[i, j, k] <- 1
                    }
                  }
                }
                else if (identical(obsmod, "Ypois")) {
                  if (Yobs[i, j, k] > 0) {
                    IDed <- sum(rbinom(Yobs[i, j, k], 1, rat))
                    YknownR[i, j, k] <- IDed
                    if (IDed != Yobs[i, j, k]) {
                      counter[i, j, k] <- Yobs[i, j, k] - IDed
                    }
                  }
                }
            }
        }
    }
    n <- n - apply(counter, 2:3, sum)
    if (tel > 0) {
        itel <- sort(sample(1:knownID, tel, replace = F))
        locs <- list()
        for (i in 1:tel) {
            lx <- rnorm(nlocs, S[itel[i], 1], sigma)
            ly <- rnorm(nlocs, S[itel[i], 2], sigma)
            locs[[i]] <- cbind(lx, ly)
        }
    }
    else {
        locs <- NULL
        itel <- NULL
    }
    list(n = n, Y = Y, Yknown = Yknown, Yobs = Yobs, YknownR = YknownR,
        counter = sum(counter), locs = locs, telID = itel, S=S)
}


hra<-function (func, parms, plot = TRUE, xlim, ylim, ng = 100, target.area = NULL, 
    tol = 0.001) 
{
    s <- c((xlim[2] - xlim[1])/2, (ylim[2] - ylim[1])/2)
    x1 <- rep(seq(xlim[1], xlim[2], , ng), ng)
    x2 <- sort(rep(seq(ylim[1], ylim[2], , ng), ng))
    delta <- min(diff(x1[1:10]))
    x1 <- rep(seq(xlim[1] - delta/2, xlim[2] + delta/2, delta), 
        ng)
    x2 <- sort(rep(seq(ylim[1] - delta/2, ylim[2] + delta/2, 
        delta), ng))
    X <- cbind(x1, x2)
    D <- sqrt((s[1] - x1)^2 + (s[2] - x2)^2)
    p <- func(parms, D)
    if (plot) {
        spatial.plot(X, p)
    }
    psi <- p/sum(p)
    if (is.null(target.area)) {
        x0 <- 0.2
        repeat {
            in.hr <- D <= x0
            total <- sum(psi[in.hr])
            if (total >= 0.95) {
                print(x0)
                break
            }
            x0 <- x0 * (1 + tol)
        }
        radius <- x0
        cat("radio para alcanzar el 95% del area: ", radius, fill = TRUE)
        area <- pi * radius^2
        cat("Area de campeo: ", area, fill = TRUE)
        return(area)
    }
    if (!is.null(target.area)) {
        if (is.null(target.area)) {
            cat("need target.area", fill = TRUE)
            goose.egg <- NULL
            return(goose.egg)
        }
        obj <- function(beta2) {
            p <- func(c(parms[1], beta2), D)
            psi <- p/sum(p)
            x0 <- 0.1
            repeat {
                in.hr <- D <= x0
                total <- sum(psi[in.hr])
                if (total >= 0.95) {
                  break
                }
                x0 <- x0 * (1 + tol)
            }
            hr.area <- pi * x0 * x0
            ss <- (hr.area - target.area)^2
            ss
        }
        tmp <- optimize(obj, interval = c(0.01, 5))
        beta2 <- tmp$minimum
        cat("Valor del parm[2] para alcanzar el 95% del area de campeo ", 
            target.area, " : ", beta2, fill = TRUE)
        return(beta2)
    }
}

spatial.plot<- function (x, y, add = FALSE, cx = 1, col = "gray") 
{
    nc <- as.numeric(cut(y, 10))
    if (!add) 
        plot(x, pch = " ", asp = 1)
    if (col == "gray") {
        cc <- seq(3, 17, , 10)/20
        cc <- gray(cc)
    }
    else cc <- terrain.colors(10)
    points(x, pch = 20, col = cc[nc], cex = cx)
    image.scale(y, col = cc)
}

#Hay 5 tipos de spiderplot:
# spiderplot: clasico
# spiderplotJJ: posibilidad de buffer alrededor de las trampas
# spiderplotJJ2: diferentes colores por animal, diferentes grosores en los segmentos y buffer.
# spiderplotJJ3: spiderplot a un plot existente, con diferentes colores, grosores y buffer
# spiderplotJJ4: spiderplot clasico a un plot existente
# spiderplotJJ5: spiderplot con puntos colores variables y conector gris a un plot existente

spiderplot<-function (y, traplocs)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy
    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ<-function (y, traplocs, buffer)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy
    
    Xl<-min(traplocs[,1])-buffer
    Xu<-max(traplocs[,1])+buffer
    Yl<-min(traplocs[,2])-buffer
    Yu<-max(traplocs[,2])+buffer
    xlim<-c(Xl,Xu)
    ylim<-c(Yl,Yu)
    
    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5, xlim=xlim, ylim=ylim)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ2<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
            las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
            cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ3<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ4<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = 1, lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = 'violet')
    }
    par(op)
}

spiderplotJJ5<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = "grey", lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)                            # Tamaño  # fondo  # grosor
        points(mu.x[[s]], mu.y[[s]], pch = 21, cex = 1.75, bg = clr, lwd=2.25)
    }
    par(op)
}


SCR23darray<-function (edf, tdf)
{
    nind <- max(edf[, 2])
    ntraps <- nrow(tdf)
    nperiods <- ncol(tdf) - 3
    per.id <- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    ind.id <- edf[, 2]
    trap.id <- edf[, 4]
    if (length(per.id) != length(min(per.id):max(per.id))) {
        x <- 1:nperiods
        names(x) <- as.character(per.id)
        per.id <- x[as.character(edf[, 3])]
    }
    else {
        per.id <- edf[, 3]
    }
    y <- array(0, c(nind, ntraps, nperiods))
    tmp <- cbind(ind.id, trap.id, per.id)
    y[tmp] <- 1
    y
}

SCRsmy<-function (y3d)
{
    nind <- dim(y3d)[1]
    totcaps <- nperiods <- sprecaps <- rep(NA, nind)
    for (i in 1:nind) {
        x <- y3d[i, , ]
        ntraps <- sum(apply(x, 2, sum) > 0)
        ncaps <- sum(x)
        nperiods[i] <- sum(apply(x, 1, sum) > 0)
        sprecaps[i] <- ifelse(ntraps > 1, 1, 0) * ncaps
        totcaps[i] <- sum(x)
    }
    cat("Total de capturas: ", sum(totcaps), fill = TRUE)
    cat("Recapturas espaciales: ", sum(sprecaps), fill = TRUE)
    cat("Eventos de captura ordinarios: ", sum(nperiods), fill = TRUE)
    cat("Capturas perdidas en el modelo no espacial: ", sum(totcaps) -
        sum(nperiods), fill = TRUE)
}

scrPID<-function (n, X, y, M, obsmod = c("pois", "bern"), niters, npics,
    xlims, ylims, a, b, inits, delta)
{
    obsmod <- match.arg(obsmod)
    J <- nrow(n)
    K <- ncol(n)
    S <- inits$S
    D <- e2dist(S, X)
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    psi <- inits$psi
    z <- rbinom(M, 1, psi)
    Y <- array(NA, c(M, J, K))
    nMarked <- nrow(y)
    marked <- rep(FALSE, M)
    marked[1:nMarked] <- TRUE
    Y[1:nMarked, , ] <- y
    z[marked] <- 1
    Ydata <- !is.na(Y)
    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            unmarked <- !Ydata[, j, k]
            nUnknown <- n[j, k]
            probs <- lam[, j] * z
            probs <- probs[unmarked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois"))
                Y[unmarked, j, k] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[unmarked, j, k] <- 0
                guys <- sample(which(unmarked), nUnknown, prob = probs)
                Y[guys, j, k] <- 1
            }
        }
    }
    cr <- rep(1, M)
    if (missing(npics)) {
        crat <- 1
    }
    else {
        crat <- npics[1]/npics[2]
    }
    cr[marked] <- crat
    out <- matrix(NA, nrow = niters, ncol = 5)
    colnames(out) <- c("sigma", "lam0", "c", "psi", "N")
    cat("\nValores de inicio =", c(sigma, lam0, crat, psi, sum(z)),
        "\n\n")
    for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"),
                "\n")
            cat("   actual =", out[iter - 1, ], "\n")
        }
        if (identical(obsmod, "pois")) {
            ll <- sum(dpois(Y, lam * cr * z, log = TRUE))
        }
        else if (identical(obsmod, "bern")) {
            ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE))
        }
        if (!missing(npics)) {
            crat <- rbeta(1, 1 + npics[1], 1 + npics[2] - npics[1])
            cr[marked] <- crat
        }
        sigma.cand <- rnorm(1, sigma, delta[1])
        if (sigma.cand > 0) {
            if (!missing(a) && !missing(b)) {
                prior <- dgamma(sigma, a, b, log = TRUE)
                prior.cand <- dgamma(sigma.cand, a, b, log = TRUE)
            }
            else {
                prior <- prior.cand <- 0
            }
            lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand *
                sigma.cand))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam * cr * z, log = TRUE))
                llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE))
                llcand <- sum(dbinom(Y, 1, lam.cand * cr * z,
                  log = TRUE))
            }
            if (runif(1) < exp((llcand + prior.cand) - (ll +
                prior))) {
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }
        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern"))
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma *
                sigma))
            if (identical(obsmod, "pois"))
                llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE))
            else if (identical(obsmod, "bern"))
                llcand <- sum(dbinom(Y, 1, lam.cand * cr * z,
                  log = TRUE))
            if (runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }
        zUps <- 0
        seen <- apply(Y > 0, 1, any)
        for (i in 1:M) {
            if (seen[i] | marked[i])
                next
            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand,
                  log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i],
                  log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] *
                  zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
            if (runif(1) < exp((llcand + prior.cand) - (ll +
                prior))) {
                z[i] <- zcand
                zUps <- zUps + 1
            }
        }
        for (j in 1:J) {
            zip <- lam[, j] * z
            for (k in 1:K) {
                if (n[j, k] == 0) {
                  Y[!marked, j, k] <- 0
                  next
                }
                unmarked <- !Ydata[, j, k]
                nUnknown <- n[j, k]
                if (nUnknown == 0)
                  next
                probs <- zip[unmarked]
                probs <- probs/sum(probs)
                if (identical(obsmod, "pois"))
                  Y[unmarked, j, k] <- rmultinom(1, nUnknown,
                    probs)
                else if (identical(obsmod, "bern")) {
                  Y[unmarked, j, k] <- 0
                  guy <- sample(which(unmarked), nUnknown, prob = probs)
                  Y[guy, j, k] <- 1
                }
            }
        }
        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + (M - sum(marked)) -
            sum(z[!marked]))
        Sups <- 0
        for (i in 1:M) {
            Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1,
                S[i, 2], delta[3]))
            inbox <- Scand[1] >= xlims[1] & Scand[1] <= xlims[2] &
                Scand[2] >= ylims[1] & Scand[2] <= ylims[2]
            if (!inbox)
                next
            dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] -
                X[, 2])^2)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2 * sigma *
                sigma))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam.cand * cr[i] *
                  z[i], log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam.cand *
                  cr[i] * z[i], log = TRUE))
            }
            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i, ] <- Scand
                lam[i, ] <- lam.cand
                D[i, ] <- dtmp
                Sups <- Sups + 1
            }
        }
        if (iter%%100 == 0) {
            cat("   Ratios de aceptacion\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, sum(z))
    }
    return(out)
}

scrPID.um<-function (n, X, y, M, mmax, obsmod = c("pois", "bern"), niters,
    npics, xlims, ylims, inits, delta)
{
    obsmod <- match.arg(obsmod)
    J <- nrow(n)
    K <- ncol(n)
    S <- inits$S
    D <- e2dist(S, X)
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    nObs <- nrow(y)
    Y <- array(0, c(M + mmax, J, K))
    Y[1:nObs, , ] <- y
    marked <- rep(FALSE, M + mmax)
    marked[1:mmax] <- TRUE
    psi <- inits$psi
    psim <- inits$psim
    z <- rbinom(M + mmax, 1, psi)
    z[1:nObs] <- 1
    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            nUnknown <- n[j, k]
            probs <- lam[!marked, j] * z[!marked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois"))
                Y[!marked, j, k] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[!marked, j, k] <- 0
                guys <- sample((mmax + 1):(M + mmax), nUnknown,
                  prob = probs)
                Y[guys, j, k] <- 1
            }
        }
    }
    cr <- rep(1, M + mmax)
    if (missing(npics)) {
        crat <- 1
    }
    else {
        crat <- npics[1]/npics[2]
    }
    cr[marked] <- crat
    out <- matrix(NA, nrow = niters, ncol = 7)
    colnames(out) <- c("sigma", "lam0", "c", "psi", "psim", "m",
        "N")
    cat("\nValores de inicio =", c(sigma, lam0, crat, psi, psim,
        sum(z[marked]), sum(z)), "\n\n")
    for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"),
                "\n")
            cat("   actual =", out[iter - 1, ], "\n")
        }
        if (identical(obsmod, "pois")) {
            ll <- sum(dpois(Y, lam * cr * z, log = TRUE))
        }
        else if (identical(obsmod, "bern")) {
            ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE))
        }
        if (!missing(npics)) {
            crat <- rbeta(1, 1 + npics[1], 1 + npics[2] - npics[1])
            cr[marked] <- crat
        }
        sigma.cand <- rnorm(1, sigma, delta[1])
        if (sigma.cand > 0) {
            lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand *
                sigma.cand))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam * cr * z, log = TRUE))
                llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE))
                llcand <- sum(dbinom(Y, 1, lam.cand * cr * z,
                  log = TRUE))
            }
            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }
        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern"))
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma *
                sigma))
            if (identical(obsmod, "pois"))
                llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE))
            else if (identical(obsmod, "bern"))
                llcand <- sum(dbinom(Y, 1, lam.cand * cr * z,
                  log = TRUE))
            if (runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }
        zUpsm <- zUps <- 0
        for (i in (nObs + 1):mmax) {
            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                llz <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcandz <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] *
                  zcand, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                llz <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcandz <- sum(dbinom(Y[i, , ], 1, lam[i, ] *
                  cr[i] * zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psim, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psim, log = TRUE)
            if (runif(1) < exp((llcandz + prior.cand) - (llz +
                prior))) {
                z[i] <- zcand
                zUpsm <- zUpsm + 1
            }
        }
        seen <- apply(Y > 0, 1, any)
        for (i in (mmax + 1):(M + mmax)) {
            if (seen[i])
                next
            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand,
                  log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i],
                  log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] *
                  zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
            rat <- (llcand + prior.cand) - (ll + prior)
            if (runif(1) < exp(rat)) {
                z[i] <- zcand
                zUps <- zUps + 1
            }
        }
        for (j in 1:J) {
            zip <- lam[!marked, j] * z[!marked]
            for (k in 1:K) {
                if (n[j, k] == 0) {
                  Y[!marked, j, k] <- 0
                  next
                }
                nUnknown <- n[j, k]
                probs <- zip/sum(zip)
                if (identical(obsmod, "pois"))
                  Y[!marked, j, k] <- rmultinom(1, nUnknown,
                    probs)
                else if (identical(obsmod, "bern")) {
                  Y[!marked, j, k] <- 0
                  guy <- sample((mmax + 1):(M + mmax), nUnknown,
                    prob = probs)
                  Y[guy, j, k] <- 1
                }
            }
        }
        psim <- rbeta(1, 1 + sum(z[marked]), 1 + mmax - sum(z[marked]))
        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + M - sum(z[!marked]))
        Sups <- 0
        for (i in 1:(M + mmax)) {
            Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1,
                S[i, 2], delta[3]))
            inbox <- Scand[1] >= xlims[1] & Scand[1] <= xlims[2] &
                Scand[2] >= ylims[1] & Scand[2] <= ylims[2]
            if (!inbox)
                next
            dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] -
                X[, 2])^2)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2 * sigma *
                sigma))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam.cand * cr[i] *
                  z[i], log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] *
                  z[i], log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam.cand *
                  cr[i] * z[i], log = TRUE))
            }
            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i, ] <- Scand
                lam[i, ] <- lam.cand
                D[i, ] <- dtmp
                Sups <- Sups + 1
            }
        }
        if (iter%%100 == 0) {
            cat("   Ratios de aceptacion\n")
            cat("     z =", zUps/M, "\n")
            cat("     zm =", zUpsm/mmax, "\n")
            cat("     S =", Sups/(M + mmax), "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, psim, sum(z[marked]),
            sum(z))
    }
    return(out)
}

scrPID.tel<-function (n, X, y, M, locs, telID, obsmod = c("pois", "bern"),
    npics, niters, xlims, ylims, a, b, inits, delta)
{
    library(mvtnorm)
    obsmod <- match.arg(obsmod)
    R <- nrow(n)
    T <- ncol(n)
    S <- inits$S
    Sin <- t(sapply(locs, colMeans))
    S[telID, ] <- Sin
    D <- e2dist(S, X)
    ntot <- length(locs)
    tel <- rep(FALSE, M)
    tel[telID] <- TRUE
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    lam[lam == 0] <- 1e-300
    psi <- inits$psi
    z <- rbinom(M, 1, psi)
    Y <- array(NA, c(M, R, T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if (!missing(y)) {
        nMarked <- nrow(y)
        marked[1:nMarked] <- TRUE
        Y[1:nMarked, , ] <- y
    }
    z[marked] <- 1
    Ydata <- !is.na(Y)
    for (r in 1:R) {
        for (t in 1:T) {
            if (n[r, t] == 0) {
                Y[!marked, r, t] <- 0
                next
            }
            unmarked <- !Ydata[, r, t]
            nUnknown <- n[r, t]
            if (nUnknown < 0)
                browser()
            probs <- lam[, r] * z
            probs <- probs[unmarked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois"))
                Y[unmarked, r, t] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[unmarked, r, t] <- 0
                guys <- sample(which(unmarked), nUnknown, prob = probs)
                Y[guys, r, t] <- 1
            }
        }
    }
    cr <- rep(1, M)
    if (missing(npics)) {
        crat <- 1
    }
    else {
        crat <- npics[1]/npics[2]
    }
    cr[marked] <- crat
    out <- matrix(NA, nrow = niters, ncol = 5)
    colnames(out) <- c("sigma", "lam0", "c", "psi", "N")
    cat("\nValores de inicio =", c(sigma, lam0, crat, psi, sum(z)),
        "\n\n")
    for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"),
                "\n")
            cat("   actual =", out[iter - 1, ], "\n")
        }
        if (identical(obsmod, "pois")) {
            ll <- sum(dpois(Y, lam * cr * z, log = TRUE))
        }
        else if (identical(obsmod, "bern")) {
            ll <- sum(dbinom(Y, 1, lam * cr * z, log = TRUE))
        }
        sigma.cand <- rnorm(1, sigma, delta[1])
        if (sigma.cand > 0) {
            lls <- lls.cand <- rep(NA, ntot)
            for (x in 1:ntot) {
                lls[x] <- sum(dmvnorm(x = locs[[x]], mean = c(S[telID[x],
                  1], S[telID[x], 2]), sigma = cbind(c(sigma^2,
                  0), c(0, sigma^2)), log = T))
                lls.cand[x] <- sum(dmvnorm(x = locs[[x]], mean = c(S[telID[x],
                  1], S[telID[x], 2]), sigma = cbind(c(sigma.cand^2,
                  0), c(0, sigma.cand^2)), log = T))
            }
            if (runif(1) < exp(sum(lls.cand) - sum(lls))) {
                sigma <- sigma.cand
                lam <- lam0 * exp(-(D * D)/(2 * sigma.cand *
                  sigma.cand))
            }
        }
        if (!missing(npics)) {
            crat <- rbeta(1, 1 + npics[1], 1 + npics[2] - npics[1])
            cr[marked] <- crat
        }
        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern"))
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma *
                sigma))
            if (identical(obsmod, "pois"))
                llcand <- sum(dpois(Y, lam.cand * cr * z, log = TRUE))
            else if (identical(obsmod, "bern"))
                llcand <- sum(dbinom(Y, 1, lam.cand * cr * z,
                  log = TRUE))
            if (runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }
        zUps <- 0
        seen <- apply(Y > 0, 1, any)
        for (i in 1:M) {
            if (seen[i] | marked[i])
                next
            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand,
                  log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i],
                  log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] *
                  zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
            if (runif(1) < exp((llcand + prior.cand) - (ll +
                prior))) {
                z[i] <- zcand
                zUps <- zUps + 1
            }
        }
        for (r in 1:R) {
            zip <- lam[, r] * z
            for (t in 1:T) {
                if (n[r, t] == 0) {
                  Y[!marked, r, t] <- 0
                  next
                }
                unmarked <- !Ydata[, r, t]
                nUnknown <- n[r, t]
                if (nUnknown == 0)
                  next
                probs <- zip[unmarked]
                probs <- probs/sum(probs)
                if (identical(obsmod, "pois"))
                  Y[unmarked, r, t] <- rmultinom(1, nUnknown,
                    probs)
                else if (identical(obsmod, "bern")) {
                  Y[unmarked, r, t] <- 0
                  guy <- sample(which(unmarked), nUnknown, prob = probs)
                  Y[guy, r, t] <- 1
                }
            }
        }
        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + (M - sum(marked)) -
            sum(z[!marked]))
        Sups <- Skups <- 0
        for (i in 1:M) {
            if (tel[i]) {
                Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1,
                  S[i, 2], delta[3]))
            }
            else {
                Scand <- c(rnorm(1, S[i, 1], delta[4]), rnorm(1,
                  S[i, 2], delta[4]))
            }
            inbox <- Scand[1] >= xlims[1] & Scand[1] <= xlims[2] &
                Scand[2] >= ylims[1] & Scand[2] <= ylims[2]
            if (!inbox)
                next
            dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] -
                X[, 2])^2)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2 * sigma *
                sigma))
            if (tel[i]) {
                ll <- sum(dmvnorm(x = locs[[sum(tel[1:i])]],
                  mean = c(S[i, 1], S[i, 2]), sigma = cbind(c(sigma^2,
                    0), c(0, sigma^2)), log = T))
                llcand <- sum(dmvnorm(x = locs[[sum(tel[1:i])]],
                  mean = c(Scand[1], Scand[2]), sigma = cbind(c(sigma^2,
                    0), c(0, sigma^2)), log = T))
                if (runif(1) < exp(llcand - ll)) {
                  ll <- llcand
                  S[i, ] <- Scand
                  lam[i, ] <- lam.cand
                  D[i, ] <- dtmp
                  Skups <- Skups + 1
                }
            }
            else {
                if (identical(obsmod, "pois")) {
                  ll <- sum(dpois(Y[i, , ], lam[i, ] * cr[i] *
                    z[i], log = TRUE))
                  llcand <- sum(dpois(Y[i, , ], lam.cand * cr[i] *
                    z[i], log = TRUE))
                }
                else if (identical(obsmod, "bern")) {
                  ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * cr[i] *
                    z[i], log = TRUE))
                  llcand <- sum(dbinom(Y[i, , ], 1, lam.cand *
                    cr[i] * z[i], log = TRUE))
                }
                if (runif(1) < exp(llcand - ll)) {
                  ll <- llcand
                  S[i, ] <- Scand
                  lam[i, ] <- lam.cand
                  D[i, ] <- dtmp
                  Sups <- Sups + 1
                }
            }
        }
        if (iter%%100 == 0) {
            cat("   Ratios de aceptacion\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/(M - length(locs)), "\n")
            cat("     Sk =", Skups/length(locs), "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, sum(z))
    }
    return(out)
}



# Simula un SMR con identificacion parcial y telemetria
sim.pID.data<-function (N = N, K = K, sigma = sigma, lam0 = lam0, knownID = knownID,
    X = X, xlims = xlims, ylims = ylims, obsmod = c("pois", "bern"),
    nmarked = c("known", "unknown"), rat = 1, tel = 0, nlocs = 0)
{
    if (tel > knownID)
        stop("Los animales con telemetria deben ser menores a knownID")
    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)
    npts <- dim(X)[1]
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    D <- e2dist(S, X)
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    Y <- array(NA, c(N, npts, K))
    for (i in 1:N) {
        for (j in 1:npts) {
            if (identical(obsmod, "bern")) {
                Y[i, j, ] <- rbinom(K, 1, lam[i, j])
            }
            else if (identical(obsmod, "pois")) {
                Y[i, j, ] <- rpois(K, lam[i, j])
            }
        }
    }
    n <- apply(Y, c(2, 3), sum)
    Yknown <- Y[1:knownID, , ]
    if (identical(nmarked, "unknown")) {
        iobs <- which(apply(Yknown > 0, 1, any))
        Yobs <- Yknown[iobs, , ]
    }
    else if (identical(nmarked, "known")) {
        Yobs <- Yknown
    }
    YknownR <- Yobs
    counter <- array(0, c(dim(Yobs)[1], dim(X)[1], K))
    for (i in 1:dim(Yobs)[1]) {
        for (j in 1:dim(X)[1]) {
            for (k in 1:K) {
                if (identical(obsmod, "bern")) {
                  if (YknownR[i, j, k] == 1) {
                    IDed <- rbinom(1, 1, rat)
                    if (IDed == 0) {
                      YknownR[i, j, k] <- 0
                      counter[i, j, k] <- 1
                    }
                  }
                }
                else if (identical(obsmod, "Ypois")) {
                  if (Yobs[i, j, k] > 0) {
                    IDed <- sum(rbinom(Yobs[i, j, k], 1, rat))
                    YknownR[i, j, k] <- IDed
                    if (IDed != Yobs[i, j, k]) {
                      counter[i, j, k] <- Yobs[i, j, k] - IDed
                    }
                  }
                }
            }
        }
    }
    n <- n - apply(counter, 2:3, sum)
    if (tel > 0) {
        itel <- sort(sample(1:knownID, tel, replace = F))
        locs <- list()
        for (i in 1:tel) {
            lx <- rnorm(nlocs, S[itel[i], 1], sigma)
            ly <- rnorm(nlocs, S[itel[i], 2], sigma)
            locs[[i]] <- cbind(lx, ly)
        }
    }
    else {
        locs <- NULL
        itel <- NULL
    }
    list(n = n, Y = Y, Yknown = Yknown, Yobs = Yobs, YknownR = YknownR,
        counter = sum(counter), locs = locs, telID = itel, S=S)
}



SCR0bayesJJ<-function (dataobj, M = 200, engine = "jags", ni = 2000, nb = 1000)
{
    if (sum(engine == c("jags", "winbugs")) == 0) {
        return("use jags o winbugs!")
    }
    y <- data$Y
    if (length(dim(y)) != 2)
    return("Los datos deben ser una matriz 2D, nind x ntraps")
        traplocs <- data$traplocs
    nind <- nrow(y)
    X <- data$traplocs
    K <- data$K
    J <- nrow(X)
    xlim <- data$xlim
    ylim <- data$ylim
    area <- (max(xlim) - min(xlim)) * (max(ylim) - min(ylim))
    y <- rbind(y, matrix(0, nrow = M - nind, ncol = ncol(y)))
    z <- c(rep(1, nind), rep(0, M - nind))
    cat("\nmodel {\nalpha0~dnorm(0,.1)\nlogit(p0)<- alpha0\nalpha1~dnorm(0,.1)\npsi~dunif(0,1)\n\nfor(i in 1:M){\n z[i] ~ dbern(psi)\n s[i,1]~dunif(xlim[1],xlim[2])\n s[i,2]~dunif(ylim[1],ylim[2])\n for(j in 1:J){\n    d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)\n    y[i,j] ~ dbin(p[i,j],K)\n    p[i,j]<- z[i]*p0*exp(- alpha1*d[i,j]*d[i,j])\n }\n}\nN<-sum(z[])\nD<- N/area\n}\n",
        file = "SCR0b.txt")
    sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1],
        ylim[2]))
    for (i in 1:nind) {
        if (sum(y[i, ]) == 0)
            next
        sst[i, 1] <- mean(X[y[i, ] > 0, 1])
        sst[i, 2] <- mean(X[y[i, ] > 0, 2])
    }
    data <- list(y = y, X = X, K = K, M = M, J = J, xlim = xlim,
        ylim = ylim, area = area)
    inits <- function() {
        list(alpha0 = rnorm(1, -4, 0.4), alpha1 = runif(1, 1,
            2), s = sst, z = z)
    }
    parameters <- c("alpha0", "alpha1", "N", "D")
    nthin <- 1
    nc <- 3
    if (engine == "winbugs") {
        library("R2WinBUGS")
        out <- bugs(data, inits, parameters, "SCR0b.txt", n.thin = nthin,
            n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE,
            working.dir = getwd())
    }
    if (engine == "jags") {
        library("rjags")
        jm <- jags.model("SCR0b.txt", data = data, inits = inits,
            n.chains = nc, n.adapt = nb)
        out <- coda.samples(jm, parameters, n.iter = ni - nb,
            thin = nthin)
    }
    return(out)
}

scrUN<-function (n, X, M, updateY = TRUE, niters, xlims, ylims, inits = NULL,
    priors = NULL, tune = c(0.2, 0.1, 2))
{
    obsmod <- "pois"
    J <- nrow(n)
    K <- ncol(n)
    if (!is.null(inits)) {
        for (i in 1:length(inits)) assign(names(inits)[i], inits[[i]])
        pn <- c("sigma", "lam0", "s", "z")
        if (updateY)
            pn <- c(pn, "y")
        noi <- !(pn %in% ls())
        if (any(noi)) {
            cat("\nAtencion: Generando valores iniciales para:",
                pn[noi], "\n")
        }
    }
    if (!("s" %in% ls()))
        s <- cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1],
            ylims[2]))
    if (!("sigma" %in% ls()))
        sigma <- runif(1, min(abs(diff(X[, 1]))), max(abs(diff(X[,
            1]))))
    if (!("lam0" %in% ls()))
        lam0 <- runif(1, 0.1, 0.5)
    if (!("psi" %in% ls()))
        psi <- runif(1, 0.9, 0.99)
    dist <- e2dist(s, X)
    lam <- lam0 * exp(-(dist * dist)/(2 * sigma * sigma))
    if (!("y" %in% ls())) {
        y <- array(0L, c(M, J, K))
        for (j in 1:J) {
            for (k in 1:K) {
                y[sample(1:M, n[j, k], replace = TRUE, lam[,
                  j]/sum(lam[, j])), j, k] <- 1
            }
        }
    }
    if (!("z" %in% ls()))
        z <- ifelse(rowSums(y) > 0, 1L, 0L)
    if (!all(dim(s) == c(M, 2)))
        stop("Las dimensiones de 2 deben ser ", c(M, 2))
    if (!all(dim(y) == c(M, J, K)) & updateY)
        stop("Las dimensiones de y deben ser ", c(M, J, K))
    if (length(z) != M)
        stop("la longitud de z debe ser ", M)
    prior.names <- names(priors)
    if (any(!(prior.names %in% c("sigma", "lam0", "psi"))))
        stop("los nombres (probabilidad a priori) debe ser uno de: 'sigma', 'lam0', 'psi'")
    if ("sigma" %in% prior.names) {
        prior.sigma <- TRUE
        sigma.prior <- function(sig) {
            prior.args <- priors[["sigma"]][[2]]
            prior.args$x <- sig
            prior.args$log <- TRUE
            do.call(priors[["sigma"]][[1]], prior.args)
        }
    }
    else prior.sigma <- FALSE
    if ("lam0" %in% prior.names) {
        prior.lam0 <- TRUE
        lam0.prior <- function(lam0) {
            prior.args <- priors[["lam0"]][[2]]
            prior.args$x <- lam0
            prior.args$log <- TRUE
            do.call(priors[["lam0"]][[1]], prior.args)
        }
    }
    else prior.lam0 <- FALSE
    if ("psi" %in% prior.names)
        prior.psi <- TRUE
    else prior.psi <- FALSE
    if (prior.psi) {
        psi.a <- priors[["psi"]][[2]][["shape1"]]
        psi.b <- priors[["psi"]][[2]][["shape2"]]
    }
    else {
        psi.a <- 1
        psi.b <- 1
    }
    out <- matrix(NA, nrow = niters, ncol = 4)
    nought <- ifelse(obsmod == "pois", "lam0", "p0")
    colnames(out) <- c("sigma", nought, "psi", "N")
    cat("\nValores de inicio =", c(sigma, lam0, psi, sum(z)), "\n\n")
    if (updateY) {
        for (iter in 1:niters) {
            if (iter%%100 == 0) {
                cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"),
                  "\n")
                cat("   actual =", out[iter - 1, ], "\n")
            }
            ll <- sum(dpois(y, lam * z, log = TRUE))
            sigma.cand <- rnorm(1, sigma, tune[1])
            if (sigma.cand > 0) {
                if (prior.sigma) {
                  prior <- sigma.prior(sigma)
                  prior.cand <- sigma.prior(sigma.cand)
                }
                else {
                  prior <- prior.cand <- 0
                }
                lam.cand <- lam0 * exp(-(dist * dist)/(2 * sigma.cand *
                  sigma.cand))
                llcand <- sum(dpois(y, lam.cand * z, log = TRUE))
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  ll <- llcand
                  lam <- lam.cand
                  sigma <- sigma.cand
                }
            }
            lam0.cand <- rnorm(1, lam0, tune[2])
            if (lam0.cand >= 0) {
                lam.cand <- lam0.cand * exp(-(dist * dist)/(2 *
                  sigma * sigma))
                llcand <- sum(dpois(y, lam.cand * z, log = TRUE))
                if (prior.lam0) {
                  prior <- lam0.prior(lam0)
                  prior.cand <- lam0.prior(lam0.cand)
                }
                else {
                  prior <- prior.cand <- 0
                }
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  ll <- llcand
                  lam0 <- lam0.cand
                  lam <- lam.cand
                }
            }
            zUps <- 0
            seen <- rowSums(y) > 0
            for (i in 1:M) {
                if (seen[i])
                  next
                zcand <- z
                zcand[i] <- ifelse(z[i] == 0, 1, 0)
                ll <- sum(dpois(y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(y[i, , ], lam[i, ] * zcand[i],
                  log = TRUE))
                prior <- dbinom(z[i], 1, psi, log = TRUE)
                prior.cand <- dbinom(zcand[i], 1, psi, log = TRUE)
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  z[i] <- zcand[i]
                  zUps <- zUps + 1
                }
            }
            for (j in 1:J) {
                probs <- lam[, j] * z
                for (k in 1:K) {
                  if (n[j, k] == 0) {
                    y[, j, k] <- 0
                    next
                  }
                  probs <- probs/sum(probs)
                  y[, j, k] <- rmultinom(1, n[j, k], probs)
                }
            }
            psi <- rbeta(1, psi.a + sum(z), psi.b + M - sum(z))
            sups <- 0
            for (i in 1:M) {
                scand <- c(rnorm(1, s[i, 1], tune[3]), rnorm(1,
                  s[i, 2], tune[3]))
                inbox <- (scand[1] >= xlims[1] & scand[1] <=
                  xlims[2]) & (scand[2] >= ylims[1] & scand[2] <=
                  ylims[2])
                if (!inbox)
                  next
                dtmp <- sqrt((scand[1] - X[, 1])^2 + (scand[2] -
                  X[, 2])^2)
                lam.cand <- lam
                lam.cand[i, ] <- lam0 * exp(-(dtmp * dtmp)/(2 *
                  sigma * sigma))
                if (z[i] == 0) {
                  ll <- llcand <- 0
                }
                else if (z[i] == 1) {
                  ll <- sum(dpois(y[i, , ], lam[i, ] * z[i],
                    log = TRUE))
                  llcand <- sum(dpois(y[i, , ], lam.cand[i, ] *
                    z[i], log = TRUE))
                }
                if (runif(1) < exp(llcand - ll)) {
                  ll <- llcand
                  s[i, ] <- scand
                  lam[i, ] <- lam.cand[i, ]
                  dist[i, ] <- dtmp
                  sups <- sups + z[i]
                }
            }
            if (iter%%100 == 0) {
                cat("   Ratios de aceptacion\n")
                cat("     z =", zUps/M, "\n")
                cat("     s =", sups/sum(z), "\n")
            }
            out[iter, ] <- c(sigma, lam0, psi, sum(z))
        }
    }
    else if (!updateY) {
        for (iter in 1:niters) {
            if (iter%%100 == 0) {
                cat("iteracion", iter, format(Sys.time(), "%H:%M:%S"),
                  "\n")
                cat("   actual =", out[iter - 1, ], "\n")
            }
            ll <- sum(dpois(n, colSums(lam * z), log = TRUE))
            sigma.cand <- rnorm(1, sigma, tune[1])
            if (sigma.cand > 0) {
                if (prior.sigma) {
                  prior <- sigma.prior(sigma)
                  prior.cand <- sigma.prior(sigma.cand)
                }
                else {
                  prior <- prior.cand <- 0
                }
                lam.cand <- lam0 * exp(-(dist * dist)/(2 * sigma.cand *
                  sigma.cand))
                llcand <- sum(dpois(n, colSums(lam.cand * z),
                  log = TRUE))
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  ll <- llcand
                  lam <- lam.cand
                  sigma <- sigma.cand
                }
            }
            lam0.cand <- rnorm(1, lam0, tune[2])
            if (lam0.cand >= 0) {
                lam.cand <- lam0.cand * exp(-(dist * dist)/(2 *
                  sigma * sigma))
                llcand <- sum(dpois(n, colSums(lam.cand * z),
                  log = TRUE))
                if (prior.lam0) {
                  prior <- lam0.prior(lam0)
                  prior.cand <- lam0.prior(lam0.cand)
                }
                else {
                  prior <- prior.cand <- 0
                }
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  ll <- llcand
                  lam0 <- lam0.cand
                  lam <- lam.cand
                }
            }
            zUps <- 0
            for (i in 1:M) {
                zcand <- z
                zcand[i] <- ifelse(z[i] == 0, 1, 0)
                llcand <- sum(dpois(n, colSums(lam * zcand),
                  log = TRUE))
                prior <- dbinom(z[i], 1, psi, log = TRUE)
                prior.cand <- dbinom(zcand[i], 1, psi, log = TRUE)
                if (runif(1) < exp((llcand + prior.cand) - (ll +
                  prior))) {
                  z[i] <- zcand[i]
                  zUps <- zUps + 1
                  ll <- llcand
                }
            }
            psi <- rbeta(1, psi.a + sum(z), psi.b + M - sum(z))
            sups <- 0
            for (i in 1:M) {
                scand <- c(rnorm(1, s[i, 1], tune[3]), rnorm(1,
                  s[i, 2], tune[3]))
                inbox <- (scand[1] >= xlims[1] & scand[1] <=
                  xlims[2]) & (scand[2] >= ylims[1] & scand[2] <=
                  ylims[2])
                if (!inbox)
                  next
                dtmp <- sqrt((scand[1] - X[, 1])^2 + (scand[2] -
                  X[, 2])^2)
                lam.cand <- lam
                lam.cand[i, ] <- lam0 * exp(-(dtmp * dtmp)/(2 *
                  sigma * sigma))
                if (z[i] == 1) {
                  ll <- sum(dpois(n, colSums(lam * z), log = TRUE))
                  llcand <- sum(dpois(n, colSums(lam.cand * z),
                    log = TRUE))
                }
                else if (z[i] == 0) {
                  ll <- llcand <- 0
                }
                if (runif(1) < exp(llcand - ll)) {
                  ll <- llcand
                  s[i, ] <- scand
                  lam[i, ] <- lam.cand[i, ]
                  dist[i, ] <- dtmp
                  sups <- sups + z[i]
                }
            }
            if (iter%%100 == 0) {
                cat("   Ratios de aceptacion\n")
                cat("     z =", zUps/M, "\n")
                cat("     s =", sups/sum(z), "\n")
            }
            out[iter, ] <- c(sigma, lam0, psi, sum(z))
        }
    }
    last <- list(sigma = sigma, lam0 = lam0, psi = psi, z = z,
        s = s)
    if (updateY)
        last$y <- y
    ret <- list(sims = out, last = last)
    return(ret)
}

make.grid<-function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
nx = 40, ny = NULL, buffer = 0)
  {
    if (is.null(ny))
    ny <- nx
    if (!is.na(ll)) {
    minx <- min(ll[, 1])
    maxx <- max(ll[, 1])
    miny <- min(ll[, 2])
    maxy <- max(ll[, 2])
    bx <- (maxx - minx) * buffer
    by <- (maxy - miny) * buffer
    minx <- minx - bx
    maxx <- maxx + bx
    miny <- miny - by
    maxy <- maxy + by
  }
  x <- sort(rep(seq(minx, maxx, , nx), ny))
  y <- rep(seq(maxy, miny, , ny), nx)
  cbind(x, y)
}


rot<-function (m) 
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

make.scrFrame<-function (caphist, traps, indCovs = NULL, trapCovs = NULL, sigCovs = NULL,
    trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")
{
    if (any(is.null(caphist), is.null(traps)))
        stop("caphist and trap must be provided")
    if (!is.list(caphist))
        stop("caphist must be a list")
    n.sessions <- length(caphist)
    caphist.dimensions <- sapply(caphist, dim)
    if (nrow(caphist.dimensions) == 2)
        caphist.dimensions <- rbind(caphist.dimensions, 1)
    for (i in 1:n.sessions) {
        caphist[[i]] <- array(caphist[[i]], dim = caphist.dimensions[,
            i])
        all.zero <- apply(apply(caphist[[i]], c(1, 3), sum),
            1, sum)
        if (any(all.zero == 0)) {
            cat("At least one individual has an all-zero encounter history",
                fill = TRUE)
            cat("Make sure this is ok...", fill = TRUE)
        }
    }
    if (!is.null(indCovs)) {
        if (!is.list(indCovs))
            stop("indCovs must be a list")
        if (any(!sapply(indCovs, is.data.frame)))
            stop("indCovs must be a list of dataframes")
        if (length(indCovs) != length(caphist))
            stop("number of sessions in indCovs does not match caphist")
        check.dim <- sapply(indCovs, nrow)
        if (any(check.dim != caphist.dimensions[1, ]))
            stop("number of individuals in indCovs does not match caphist")
        if (!("rmv" %in% indCovs[[1]])) {
            for (i in 1:length(indCovs)) {
                indCovs[[i]]$removed <- dim(caphist[[i]])[3]
            }
        }
    }
    else {
        indCovs <- list()
        for (i in 1:length(caphist)) {
            indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],
                dim(caphist[[i]])[1]))
        }
    }
    if (!is.list(traps))
        stop("traps must be a list")
    if (length(traps) != length(caphist))
        stop("number of sessions in traps does not match caphist")
    check.dim <- sapply(traps, nrow)
    if (!all(check.dim == caphist.dimensions[2, ]))
        stop("number of traps does not match caphist")
    if (!is.null(trapCovs)) {
        if (!is.list(trapCovs))
            stop("trapCovs must be a list")
        if (any(!sapply(trapCovs, is.list)))
            stop("trapCovs must be a list of lists")
        if (any(!unlist(sapply(trapCovs, function(x) sapply(x,
            is.data.frame)))))
            stop("trapCovs must be a list of dataframes")
        if (length(trapCovs) != length(caphist))
            stop("number of sessions in trapCovs does not match caphist")
        check.dim <- lapply(trapCovs, function(x) sapply(x, nrow))
        for (i in 1:length(check.dim)) {
            if (!all(check.dim[[i]] == caphist.dimensions[2,
                i]))
                stop("number of traps does not match caphist")
        }
    }
    if (!is.null(sigCovs)) {
        if (nrow(sigCovs) != length(caphist))
            stop("number of rows in sigCovs does not match number of sessions")
        if (!"session" %in% colnames(sigCovs)) {
            sigCovs$session <- factor(1:n.sessions)
        }
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    else {
        sigCovs <- data.frame(session = factor(1:n.sessions))
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    if (!is.null(trapOperation)) {
        if (!is.list(trapOperation))
            stop("trapOperation must be a list")
        if (length(trapOperation) != length(caphist))
            stop("number of sessions in trapOperation does not match caphist")
        check.dim <- sapply(trapOperation, nrow)
        if (!all(check.dim == caphist.dimensions[2, ]))
            stop("number of traps does not match caphist")
    }
    max.dist <- NULL
    for (i in 1:length(caphist)) {
        for (j in 1:nrow(caphist[[i]])) {
            if (dim(caphist[[i]])[3] > 1) {
                where <- apply(caphist[[i]][j, , ], 1, sum) >
                  0
            }
            else {
                where <- caphist[[i]][j, , ] > 0
            }
            if (sum(where) > 1)
                max.dist <- c(max.dist, max(0, dist(traps[[i]][where,
                  c("X", "Y")]), na.rm = T))
        }
    }
    mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
    mdm <- max(max.dist, na.rm = T)
    if (!is.null(telemetry)) {
        if (!is.list(telemetry$fixfreq))
            stop("telemetry$fixfreq must be a list")
        fixfreq.dimensions <- sapply(telemetry$fixfreq, dim)
        if (nrow(fixfreq.dimensions) == 2)
            fixfreq.dimensions <- rbind(fixfreq.dimensions, 1)
        if (!is.null(telemetry$indCovs)) {
            if (!is.list(telemetry$indCovs))
                stop("telemetry$indCovs must be a list")
            if (any(!sapply(telemetry$indCovs, is.data.frame)))
                stop("telemetry$indCovs must be a list of dataframes")
            if (length(telemetry$indCovs) != length(telemetry$fixfreq))
                stop("number of sessions in telemetry$indCovs does not match telemetry$fixfreq")
            check.dim <- sapply(telemetry$indCovs, nrow)
            if (any(check.dim != fixfreq.dimensions[1, ]))
                stop("number of individuals in telemetry$indCovs does not match telemetry$fixfreq")
            if (any(!names(indCovs[[1]]) %in% c(names(telemetry$indCovs[[1]]),
                "removed")))
                stop("indCovs do not match between capture and telemetry data")
        }
        if (!is.null(telemetry$cap.tel)) {
            if (!is.list(telemetry$cap.tel))
                stop("telemetry$indCovs must be a list")
            warning("make sure captured individuals w/ collars sorted first!")
        }
    }
    if (!is.null(rsfDF)) {
        library(FNN)
        rsfCovs <- names(rsfDF[[1]][, -c(1:2), drop = F])
        if (is.null(trapCovs)) {
            trapCovs <- list()
            length(trapCovs) <- n.sessions
            for (s in 1:n.sessions) {
                trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                  c("X", "Y")], traps[[s]][, c("X",
                  "Y")], 1)$nn.index)
                trapCovs[[s]] <- list()
                length(trapCovs[[s]]) <- caphist.dimensions[3,
                  s]
                for (k in 1:caphist.dimensions[3, s]) {
                  trapCovs[[s]][[k]] <- data.frame(rsfDF[[s]][trap.grid,
                    rsfCovs])
                  names(trapCovs[[s]][[k]]) <- rsfCovs
                }
            }
        }
        else {
            for (s in 1:n.sessions) {
                if (any(!rsfCovs %in% trapCovs[[s]][[1]])) {
                  miss.rsfCovs <- rsfCovs[which(!rsfCovs %in%
                    trapCovs[[s]][[1]])]
                  trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                    c("X", "Y")], traps[[s]][, c("X",
                    "Y")], 1)$nn.index)
                  for (k in 1:caphist.dimensions[3, s]) {
                    newtrapCovs <- data.frame(rsfDF[[s]][trap.grid,
                      miss.rsfCovs])
                    names(newtrapCovs) <- miss.rsfCovs
                    trapCovs[[s]][[k]] <- data.frame(trapCovs[[s]][[k]],
                      newtrapCovs)
                  }
                }
            }
        }
    }
    scrFrame <- list(caphist = caphist, traps = traps, indCovs = indCovs,
        trapCovs = trapCovs, sigCovs = sigCovs, trapOperation = trapOperation,
        occasions = caphist.dimensions[3, ], type = type, mmdm = mmdm,
        mdm = mdm, telemetry = telemetry)
    class(scrFrame) <- "scrFrame"
    return(scrFrame)
}



e2dist <- function (x, y) {  # Function from scrbook package to calculate the distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

e2dist.2 <- function (x, y) {  # Function from scrbook package to calculate the squared distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- (x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

##################################################################################################################################
### 1.  Function to generate generalized SMR data with linear, gridded, or random trap (within the state-space), trap designs ####
##################################################################################################################################

fnc.create.SMR.data <- function (N = N, K.trap = K.trap, K.camera = K.camera,
                                 sigma = sigma, g0trap = g0trap, lam0 = lam0,
                                 n.trap = n.trap, n.camera = n.camera,
                                 trap.design = trap.design,
                                 xlims = xlims, ylims = ylims,
                                 xlim.R = xlim.R, ylim.R = ylim.R,
                                 obsmod = c("pois", "bern"),
                                 n.collar = 10000, nlocs = 100,
                                 plot.map = TRUE, rnd = 2018) {
  # N = True number of animals in the state-space
  # K.trap = number of marking (trapping) occasions
  # K = number of resighting (eg camera) occasions
  # sigma = scale parameter for home range size
  # g0trap = probability of detection at the home range center for marking
  # lam0 = encounter rate at home range center for resighting
  # n.trap = total number of traps for marking.  Used to create a square
  #     grid of traps
  # n.camera = total number of cameras for resighting.  Used to create a square
  #     grid of detectors
  # trap.design = one of 'linear', 'grid', 'random'.
  # xlims = min and max coordinates of easting for the state-space
  # ylims = min and max coordinates of northing for the state-space
  # xlim.R = min and max coordinates of easting for the study area
  # ylim.R = min and max coordinates of northing for the study area.
  # obsmod = type of observation model for the resighting process.  We used
  #     'pois' for the manuscript.
  # n.collar = maximum number of animals to fit with telemetry tags or collars.
  #     Telemetry data is created for the lesser of n.collar or n.marked
  # nlocs = number of telemetry locations for each marked animal.  We used 100
  #     for each marked animal in the manuscript

  # VALUES RETURNED
  # y.trap.all = all capture observations including known animals that were
  #     undetected.  Array of individuals x trap x trap.occasion
  # y.trap.marked = SCR observations of captured animals.  Array of
  #     individuals x trap x trap.occasion
  # y.camera.all = all resight observations including animals not detected
  #     and/or not trapped.  Array of individuals x detector x occasion
  # y.camera.marked = resight encounter history of marked animals.  Array of
  #     individuals x detector x occasion
  # n.camera.unmarked = number of detections of unmarked animals at location
  #     j (row) and occasion k (column)
  # i.marked = row numbers of marked animals in y.trap.all and y.camera.all
  # i.unmarked = row numbers of unmarked animals in y.camera.all
  # telemetry.array = Telemetry locations of marked individuals.  Array of
  #     individual x location number x coordinates(x,y).
  # X.trap = locations of traps for marking.  Matrix of trap x coordinates
  # X.camera = locations of detectors for resighting.  Matrix of detector x
  #     coordinates

  set.seed(rnd)

  if ( ! trap.design %in% c('grid', 'linear', 'random')){
    stop('trap.design should be one of:  grid, linear, or random')
  }

  obsmod <- match.arg(obsmod)
  n.trap.row <- sqrt(n.trap)
  n.camera.row <- sqrt(n.camera)
  if (trap.design == 'grid'){
    coor0 <- seq(xlim.R[1], xlim.R[2], length = sqrt(n.trap))
    # Nets for a square grid of traps
    X.trap <- cbind(rep(coor0, each=length(coor0)),rep(coor0,
      times=length(coor0)))
  }

  if (trap.design == 'linear'){
    coor0 <- seq(xlim.R[1], xlim.R[2], length = n.trap)
    X.trap <- cbind(coor0, 0.5)
  }

  if (trap.design == 'random') {
    n.rand <- n.trap
    X.trap <- cbind( runif(n.rand, xlims[1], xlims[2]),
                     runif(n.rand, xlims[1], xlims[2]))
  }
  J.trap <- nrow(X.trap)  # nN

  ## Camera coordinates
  coor0 <- seq(xlim.R[1], xlim.R[2], length = n.camera.row)
  # Nets for a square grid of traps
  X.camera <- cbind(rep(coor0, each=length(coor0)),rep(coor0,
    times=length(coor0)))
  t.jitter <- (coor0[2] - coor0[1])/3
  J.camera <- nrow(X.camera)  # nN

    # Activity Centers
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    #### MARK DATA
    D.trap <- e2dist(S, X.trap)
    ptrap <- g0trap * exp(-(D.trap * D.trap)/(2 * sigma * sigma))
    y.trap.all <-array(0, c(N, J.trap, K.trap))
    for (i in 1:N){
      for (j in 1:J.trap){
        y.trap.all[i, j, ] <- rbinom(K.trap, 1, ptrap[i, j])
      }
    }
    n.trap.ind <- apply(y.trap.all, 1, sum)# number of captures per individual
    marked <- ifelse(n.trap.ind > 0, 1, 0) # is each animal marked (0 or 1)
    i.marked = (1:N)[marked == 1]          # ID for marked individuals
    n.marked <- sum(marked)                # number captured and marked
    y.trap <- y.trap.all[marked == 1, , ]  # capture-recapture data for marked ind.

    #### RESIGHT DATA
    # Distance between each home range center (row) and camera (column)
    D <- e2dist(S, X.camera)
    # Encounter rates
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    # Array for resighting data
    y.camera.all <- array(NA, c(N, J.camera, K.camera))
    for (i in 1:N) {
      for (j in 1:J.camera) {
        if (identical(obsmod, "bern")) {
          y.camera.all[i, j, ] <- rbinom(K.camera, 1, lam[i, j])
        }
        else if (identical(obsmod, "pois")) {
          y.camera.all[i, j, ] <- rpois(K.camera, lam[i, j])
        }
      }
    }
    y.camera.unmarked <- y.camera.all * (1 - marked)
    i.unmarked <- (1:N)[rowSums(y.camera.unmarked) > 0]
    # Sum of detections by Camera (row) and Occasion (column)
    n.camera.unmarked <- apply(y.camera.unmarked, c(2, 3), sum)
    # Resight data of marked animals
    y.camera.marked <- y.camera.all[marked == 1, , ]
    # Number of detections per individual
    n.ind <- apply(y.camera.all, c(1), sum)
    # was each individual detected by camera. Used in plot below.
    det.camera <- ifelse(n.ind > 0, 1, 0)

    # Telemetry data
    n.collar <- min(c(n.marked, n.collar))
    telemetry.array <- array(NA, dim=c(n.collar, nlocs, 2))
    if (nlocs > 0 & n.collar > 0) {
      for (i in 1:n.collar) {
        telemetry.array[i, , 1] <- rnorm(nlocs, S[i.marked[i], 1], sigma)
        telemetry.array[i, , 2] <- rnorm(nlocs, S[i.marked[i], 2], sigma)
      }
    }

    # Plot marked and unmarked animals
    if (plot.map == TRUE){
      par(mfrow = c(1,1))
      plot(S, col = 'red', pch = 19, xlim = xlims, ylim = ylims, asp = 1,
           xlab = 'X', ylab = 'Y',
           main = paste('Spatial Mark-Resight Simulated Data \n g0trap =', g0trap,
                        ' lam0camera =', lam0, ' sigma =', sigma, ' N =', N,
                        ' \n Trap Design:', trap.design,  sep = ' '))
      points(S, col = 'red', pch = 19, cex=1.75)
      points(X.trap, col='black', pch = 17, cex = 1.2)
      points(X.camera, col = 'black', pch = '+', cex = 0.8)
      points(S[det.camera == 1, 1], S[det.camera == 1, 2], pch=19, col='green',
        cex = 1.15)
      points(S[marked == 1, 1], S[marked == 1, 2], pch=19, col='purple', cex = 0.85)
      legend('topright', cex = 0.8,
        legend = c( 'Traps for Marking', 'Cameras', 'Home Range Centre',
        '  Marked', '  Detected Camera'),
        pch = c(17, 3, 19, 19, 19), pt.cex = c(1.2, 0.8, 1.75, 0.85, 1.15),
        col = c('black','black', 'red','purple', 'green' ),
        bg = 'gray95')
    } # End of Plot
    # Collect Results
    list(y.trap.all = y.trap.all, y.trap = y.trap, i.marked = i.marked,
         y.camera.all = y.camera.all, y.camera.marked = y.camera.marked,
         i.unmarked = i.unmarked,
         n.camera.unmarked = n.camera.unmarked,
         telemetry.array = telemetry.array, X.trap = X.trap, X.camera = X.camera)

}   # END OF FUNCTION  fnc.create.SMR.data


library(coda)
library(rstan)
nimSummary = function(d = NULL, trace=FALSE, exclude.params = NULL, digits=2){
  if(is.null(exclude.params)==FALSE){
    require(stringr)
    tmp1 = ifelse(is.na(as.numeric(str_extract(attributes(d[[1]])$dimnames[[2]],"[1-9]+"))),attributes(d[[1]])$dimnames[[2]],substr(attributes(d[[1]])$dimnames[[2]],1,as.numeric(str_locate(attributes(d[[1]])$dimnames[[2]], "\\[")[, 1])-1))
    d.remove = lapply(d, function(x) which(tmp1 %in% exclude.params))
    d2 = lapply(d, function(x) x[,-d.remove[[1]]])
  }else
    if(is.null(exclude.params)){
      d2 = d
      d.remove = list()
      d.remove = 0
    }
  if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
    d3 = as.data.frame(do.call(c, d2))
    #d3 = d3[,-which(apply(d3, 2, function(x) any(x=="Inf")))]
    Means = mean(d3[,1], na.rm=TRUE)
    Median= median(d3[,1], na.rm=TRUE)
    SDs = sd(d3[,1], na.rm=TRUE)
    q2.5 = quantile(d3[,1], 0.025, na.rm=TRUE)
    q50 = quantile(d3[,1], 0.50, na.rm=TRUE)
    q97.5 = quantile(d3[,1], 0.975, na.rm=TRUE)
    over.zero = round(mean(d3[,1]>0),2)
    n.eff = effectiveSize(mcmc.list(lapply(d2, as.mcmc)))
    Rhat = round(gelman.diag(mcmc.list(lapply(d2, as.mcmc)), multivariate = FALSE)[[1]][,1],3)
  }else
    if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
      d3=do.call(rbind,d2)
      #d3 = d3[,-which(apply(d3, 2, function(x) any(x=="Inf")))]
      Means = apply(d3, 2,function(x) mean(x,na.rm=TRUE))
      Median= apply(d3, 2,function(x) median(x,na.rm=TRUE))
      SDs = apply(d3, 2,function(x) sd(x,na.rm=TRUE))
      q2.5 = apply(d3, 2,function(x) quantile(x, 0.025,na.rm=TRUE))
      q50 = apply(d3, 2,function(x) quantile(x, 0.50,na.rm=TRUE))
      q97.5 = apply(d3, 2,function(x) quantile(x, 0.975,na.rm=TRUE))
      over.zero = round(apply(d3, 2, function(x) mean(x>0,na.rm=TRUE)),2)
      n.eff = effectiveSize(mcmc.list(lapply(d2, as.mcmc)))
      Rhat = round(gelman.diag(mcmc.list(lapply(d2, as.mcmc)), multivariate = FALSE)[[1]][,1],3)
      niter = attributes(d[[1]])$dim[1]
      param = attributes(d[[1]])$dim[2]
      arr<-array(0,c(niter,param,3))
      arr[,,1] <- as.matrix(outNim[[1]])
      arr[,,2] <- as.matrix(outNim[[2]])
      arr[,,3] <- as.matrix(outNim[[3]])
      #convert array to iterations x chains x parameters
      arr <- aperm(arr, c(1,3,2))
      #Calculate split-chain Rhat and other parameters
      #no warmup/burn-in iterations saved by nimble, so set warmup to 0
      mo=rstan::monitor(arr, warmup=0, print=FALSE)
    }
  if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
    par(mfrow=c(ncol(d3),2))
    for(i in 1:ncol(d3)){
      plot(1:dim(d2[[1]])[1],d2[[1]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, lapply(d2,function(x) apply(x, 2, range)))[,i]))
      for(j in 2:length(d2)){
        lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,i],main="",xlab=colnames(d3)[i])
    }
  }else
    if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
      par(mfrow=c(1,2))
      plot(1:length(d2[[1]]),d2[[1]],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(d3[,1]))
      for(j in 2:length(d2)){
        lines(1:length(d2[[j]]),d2[[j]],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,1],main="",xlab=attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]])
    }
  tmp.frame = data.frame(Mean=Means,SD=SDs,q2.5=q2.5,q50=q50,q97.5=q97.5,f0=over.zero,n.eff=n.eff,Bulk_ESS=mo$Bulk_ESS,Tail_ESS=mo$Tail_ESS, Rhat=mo$Rhat)
  if(nrow(tmp.frame)==1){
    row.names(tmp.frame) = attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]]
  }
  return(paste(cat("For each parameter, Bulk_ESS and Tail_ESS are crude measures of effective \nsample size for bulk and tail quantities respectively (an ESS > 100 per chain \nis considered good), and Rhat is the potential scale reduction factor on rank \nnormalized split chains (at convergence, Rhat <= 1.05).\n\n"), return(round(tmp.frame, digits=digits))))
}

spiderplot.Over<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{

   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = 1, lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = 'red')
    }
    par(op)
}

