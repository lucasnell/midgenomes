

#' code: https://github.com/cathyqqtao/ddBD-tree-prior
#' permalink for ddBD.R script pasted below:
#' https://github.com/cathyqqtao/ddBD-tree-prior/blob/6c38f67dea22f118517adf8a60247f5d272a902a/code/ddBD.R
#'
#' paper: https://doi.org/10.1093/bioinformatics/btab307

#' It doesn't work on R 4.0.*, so I have to run it inside the midgenomes
#' docker container using the following command:
#'
#' ```
#' docker run -it --rm=true --platform linux/amd64 -v <midgenomes DIRECTORY>/_data:/data \
#'     lucasnell/midgenomes:v1.0.10 /bin/bash
#' ```
#'
#' Once inside the container, I run the following bash code:
#'
#' ```
#' . /app/.bashrc
#' conda activate phylo-env
#' R
#' ```



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#               CODE FROM ddBD-tree-prior/code/ddBD.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

############ PLEASE RUN THE FOLLOWING CODE ################
library("ape")
library("stats4")
library("FNN")

LL.BD = function(tt, lambda, mu, rho){
    t1 = 1
    A = exp((mu-lambda)*tt)
    A1 = exp((mu-lambda)*t1)
    Prob.t = (rho*(lambda-mu))/(rho*lambda + (lambda*(1-rho) - mu) * A)
    Prob.t1 = (rho*(lambda-mu))/(rho*lambda + (lambda*(1-rho) - mu) * A1)
    vt1 = 1 - 1/rho*Prob.t1*A1
    p1 = (1/rho)*Prob.t^2*A
    gt = lambda*p1/vt1
    -sum(log(gt))
}


BD.density = function(t, birth.rate, death.rate, rho, root.age=1){
    A = exp((death.rate - birth.rate)*t)
    A1 = exp((death.rate - birth.rate)*root.age)
    Prob.t = (rho*(birth.rate - death.rate))/(rho*birth.rate + (birth.rate*(1-rho) - death.rate) * A)
    Prob.t1 = (rho*(birth.rate - death.rate))/(rho*birth.rate + (birth.rate*(1-rho) - death.rate) * A1)
    vt1 = 1 - 1/rho * Prob.t1*A1
    p1 = 1/rho * Prob.t^2 * A
    gt = birth.rate*p1 / vt1

    return(gt)
}


ddBD = function(tr, outgroup, root.time = 1, measure = c("SSE","KL")){
    b.rate.try = seq(1,10,1)+0.1
    d.rate.try = seq(1,10,1)
    s.fr.try = c(0.001,0.01,0.1,0.5,0.9)
    paras.try = expand.grid(b.rate.try, d.rate.try, s.fr.try)
    paras.try = paras.try[-which(paras.try$Var1<paras.try$Var2), ] # force birth rate >= death rate

    tr = ape::drop.tip(tr, outgroup)
    t = ape::branching.times(tr)/max(ape::branching.times(tr))
    t[t<0] = 0
    t.den = density(t)
    t.den.x= t.den$x
    t.den.y = t.den$y
    t.den.x.2 = t.den.x[t.den.x>=0 & t.den.x<=1]
    t.den.y.2 = t.den.y[t.den.x>=0 & t.den.x<=1]

    err = numeric()
    kl.dist = numeric()

    for (i in 1:nrow(paras.try)){
        bd.density = BD.density(t = t.den.x.2, birth.rate = paras.try[i,1], death.rate = paras.try[i,2], rho = paras.try[i,3], root.age = 1)

        err = c(err, sqrt(sum((t.den.y.2-bd.density)^2)))
        kl.dist = c(kl.dist, mean(FNN::KL.divergence(t.den.y.2, bd.density, k=5)))
    }

    err.sort = sort(err)
    kl.sort = sort(kl.dist)
    attempt = 0
    inf.paras = NULL

    while (is.null(inf.paras) && attempt <= 10){
        attempt = attempt + 1

        if (measure == "KL"){
            paras.start = paras.try[match(kl.sort[attempt], kl.dist), ]
        }else{
            paras.start = paras.try[match(err.sort[attempt], err), ]
        }

        names(paras.start) = c("lambda", "mu", "rho")
        inf.paras = tryCatch(stats4::mle(LL.BD, start = list(lambda = as.numeric(paras.start[1]), mu = as.numeric(paras.start[2]), rho = as.numeric(paras.start[3])),
                                         fixed = list(tt = t), method = "L-BFGS-B", lower = c(0, 0, 0), upper = c(Inf, Inf, 1)), error=function(e){})
    }

    if (attempt >= 10){
        return("Sorry, the best parameter setting cannot be found.")
    }else{
        output = c(inf.paras@coef[1:2]/root.time, inf.paras@coef[3])

        return(output)
    }

}

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#               CODE END
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


rel_tr <- read.tree("/data/phylo/chir_mega_relTimes.nwk")

#' Notes:
#' * Because RelTime-ML in MEGA removes the outgroup from the relative-time
#'   tree, our new outgroups are from family Culicidae.
#' * The root time (in units of 100 Ma) is an approximation of the time
#'   from timetree.org.

ddBD(rel_tr, outgroup = c("Asteph", "Aaegyp", "Cquinq"), root.time = 2.2, measure = "SSE")
#    lambda        mu       rho
# 1.0694326 1.0694338 0.7733908


# # previous version:
# ddBD(rel_tr, outgroup = "Culicoides_sonorensis", root.time = 2.2, measure = "SSE")
# #    lambda        mu       rho
# # 3.6612684 3.6612434 0.2195345



