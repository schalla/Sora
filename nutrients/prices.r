require(lattice)
require(plyr)
require(reshape)
require(rjags)
require(rv)
require(latticeExtra)

setwd("~/Desktop/sora-old/nutrients/")
d <- read.csv("prices.csv")
daily <- read.csv("daily.csv")
d$X <- NULL; daily$X <- NULL

runModel <- function(d, reps = 1000) {
  data <- as.list(d)
  data$N <- dim(d)[1]
  data$L <- length(levels(d$name))
  data$T <- length(levels(d$type))

  m <- jags.model("prices.bug", data)
  update(m, 2000)
  s <- coda.samples(m,
                    c("A.global",
                      "m.global",
                      "m", "A",
                      "A.name",
                      "m.name",
                      "beta",
                      "m.sd",
                      "sigma.A",
                      "sigma.m"),
                    reps)
  return(s)
}

modelReps <- function(d, s,
                      amt.reps = c(10, 100, 500, 1000, 2000, 5000, 10000)) {

  # Do some class coersion
  if (class(s) == "mcmc.list") { s <- s[[1]] }
  if (class(s) == "mcmc") {s <- as.data.frame(s) }

  s$sample <- seq(dim(s)[1])
  
  # Create a variable matrix which can be used to build the reps
  vars <- ddply(s, .(sample), function(s1) {
    # Deconstruct the variables
    V <- splitbyname(s1)
    name.vars <- data.frame(A.name = as.double(V$A.name),
                            m.name = as.double(V$m.name),
                            name = levels(d$name))
    type.vars <- data.frame(A.type = as.double(V$A),
                            m.type = as.double(V$m),
                            type = levels(d$type))
    betas <- as.double(V$beta)
    globals <- list(A = as.double(V$A.global), m = as.double(V$m.global))

    # Rebuild factor frame
    skeleton <- unique(subset(d, select = c("name", "type", "acs", "purity")))
    out <- merge(skeleton, name.vars, by = c("name"))
    out <- merge(out, type.vars, by = c("type"))
    out$beta1 <- betas[1]
    out$beta2 <- betas[2]
    out
  })

  # Generate a reps matrix
  #  for each sample...
  reps <- ddply(vars, .(sample), function (s1) {
    # for each observation
    ddply(s1, .(name), function(obs) {
      transform(obs,
                amt = amt.reps,
                per = exp(A.type + A.name + beta1*acs) * amt.reps^(m.type + m.name + beta2*acs)
                )
    })
  })

  transform(reps, price = per*amt)
}

## makeReps <- function(d, rv, samps = 100, amt.reps = c(1, 10, 100, 1000, 10000)) {
##   link <- unique(subset(d, select = c("name", "type", "acs")))
##   vars <- splitbyname(rv)
##   n.names <- length(levels(link$name))
  
##   p <- dlply(link, .(name), function(d1) {
##     n <- d1$name
##     t <- link[link$name == n,]$type
##     acs <- link[link$name == n,]$acs

##     A.type <- vars$A[t]
##     A.name <- vars$A.name[as.numeric(n)]
##     m.type <- vars$m[t]
##     m.name <- vars$m.name[as.numeric(n)]

##     exp(A.type + A.name + vars$beta[1]*acs) *
##       amt.reps^(1 + m.type + m.name + vars$beta[2]*acs)
##   })
  
##   (out <- ldply(p, function(r) rvsample(r, samps)))
##   out$sample <- gl(samps, n.names)
##   out <- melt(out, id.var = c("name", "sample"))
##   names(out) <- c("name", "sample", "amt", "price")
##   out$amt <- factor(out$amt, labels = amt.reps)
##   out$amt <- as.numeric(as.character(out$amt))
##   merge(link, out, by = "name")
## }

compare <- function(d, sd = NULL) {
  if (is.null(sd)) {sd = runModel(d, 100)}
  reps <- modelReps(d, sd)
  reps$type <- factor(reps$type, labels = levels(d$type))
  a <- xyplot(log(price/amt) ~ log(amt) | type, d, groups = name, type = 'l', ylim = c(-10, 0), lwd = 2, xlim = c(2,10), col = 5)
  b <- xyplot(log(price/amt) ~ log(amt) | type, reps, groups = interaction(sample,name), col = 1, alpha = 1/10, lwd = 1, lty = 1, type = 'l', ylim = c(-10, 0), xlim = c(2,10))
  print(b, split = c(1, 1, 1, 1), more = TRUE)
  print(a, split = c(1, 1, 1, 1), more = FALSE)
}

test <- function(d) {
  return(rvsims(runModel(d)[[1]]))
}

evalPrice <- function(d, daily, rv) {
  link <- colwise(as.numeric)(unique(subset(d, select = c("name", "type"))))
  link <- merge(link, colwise(as.numeric)(daily), by = c("type"))
  
  vars <- splitbyname(rv)

  total = 0
  for (t in unique(link$type)) {
    ns <- subset(link, type == t)$name
    which.n <- ns[which.min(E(vars$m.name)[ns])]
    amt <- subset(link, type == t & name == which.n)$amt
    
    A.type <- vars$A[t]
    A.name <- vars$A.name[which.n]

    m.type <- vars$m[t]
    m.name <- vars$m.name[which.n]

    total = total + exp(A.type + A.name) * amt^(1 + m.type + m.name)
  }
  total
}

catalogPrices <- function(d, daily, rv) {
  link <- colwise(as.numeric)(unique(subset(d, select = c("name", "type"))))
  link <- merge(link, colwise(as.numeric)(daily), by = c("type"))
  
  vars <- splitbyname(rv)

  frame <- data.frame()
  for (t in unique(link$type)) {
    ns <- subset(link, type == t)$name
    which.n <- ns[which.min(E(vars$m.name)[ns])]
    amt <- subset(link, type == t & name == which.n)$amt
    
    A.type <- vars$A[t]
    A.name <- vars$A.name[which.n]

    m.type <- vars$m[t]
    m.name <- vars$m.name[which.n]

    unitprice <- exp(A.type + A.name) * amt^(m.type + m.name)

    col <- list(type = t,
                supplier = "Alfa Aesar",
                name = which.n,
                qty = amt,
                unit = E(unitprice),
                total = E(unitprice*amt))
    frame <- rbind(frame, col)
  }
  frame
}
