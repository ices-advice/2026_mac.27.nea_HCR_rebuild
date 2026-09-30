### ------------------------------------------------------------------------ ###
### run fmle for FLSR in parallel ####
### ------------------------------------------------------------------------ ###
fmle_parallel <- function(sr, cl, seed = NULL, ...) {

  ### split into parts
  it_parts <- split(1:dim(sr)[6], cut(1:dim(sr)[6], length(cl)))
  names(it_parts) <- NULL

  ### split sr into junks
  sr_chunks <- lapply(it_parts, function(x) {
    FLCore::iter(sr, x)
  })

  ### set RNG seed
  if (!is.null(seed)) {
    clusterEvalQ(cl = cl, expr = {set.seed(seed)})
  }

  ### fit model to each junk
  fits <- foreach(sr_chunk = sr_chunks, .packages = "FLCore") %dopar% {

    fmle(sr_chunk, ...)

  }

  ### combine into single object
  ### slots rec, ssb, covar, logerrror, model, logl, gr, distribution, inital,
  ### name, desc, range
  ### do not need to be changed
  logLik(sr) <- logLik(fits[[1]])
  sr@logLik[seq(dims(sr)$iter)] <- c(sapply(fits, logLik))

  fitted(sr) <- do.call(FLCore::combine, lapply(fits, "fitted"))

  residuals(sr) <- do.call(FLCore::combine, lapply(fits, "residuals"))

  vcov(sr) <- array(data = unlist(lapply(fits, vcov)),
                    dim = c(dim(fits[[1]]@vcov)[1],
                            dim(fits[[1]]@vcov)[2],
                            dims(sr)$iter),
                    dimnames = list(dimnames(fits[[1]]@vcov)[[1]],
                                    dimnames(fits[[1]]@vcov)[[2]],
                                    iter = ac(dimnames(sr)$iter)))

  params(sr) <- propagate(params(sr), dims(sr)$iter)
  params(sr)[] <- unlist(lapply(fits, function(x) {
    split(x@params@.Data, slice.index(x@params, MARGIN = 2))
  }))

  hessian(sr) <- array(data = unlist(lapply(fits, hessian)),
                       dim = c(dim(fits[[1]]@hessian)[1],
                               dim(fits[[1]]@hessian)[2],
                               dims(sr)$iter),
                       dimnames = list(dimnames(fits[[1]]@hessian)[[1]],
                                       dimnames(fits[[1]]@hessian)[[2]],
                                       iter = ac(dimnames(sr)$iter)))

  details(sr) <- details(fits[[1]])

  return(sr)

}

### ------------------------------------------------------------------------ ###
### function for estimating OM MSY reference values ####
### ------------------------------------------------------------------------ ###
est_MSY <- function(om, oem, control, args, parallel = FALSE,
                    path, scenario = "",
                    vals_ini = seq(0, 1, 0.1),
                    lower = 0, upper = 0.5, tol = 0.001,
                    plot = TRUE, x_label = "F (ages 4-8)",
                    save = TRUE) {

  ### check if some results already exist
  if (isTRUE(save)) {
    if (isTRUE(file.exists(paste0(path, scenario, "_MSY_trace.rds")))) {
      res_trace_ini <- readRDS(paste0(path, scenario, "_MSY_trace.rds"))
    } else {
      res_trace_ini <- list()
    }
  }
  ### create object in new environment for storing results
  trace_env <- new.env()
  assign(x = "res_trace", value = res_trace_ini, envir = trace_env)

  ### define function for running projection and returning stats
  input <- list(om = om, oem = oem, control = control, args = args,
                parallel = parallel, verbose = FALSE, progress = FALSE)
  mp_catch <- function(input, Ftrgt, minimise = FALSE) {#browser()
    res_trace_i <- get("res_trace", envir = trace_env)
    ### if run already exists, do not run again
    if (isTRUE(Ftrgt %in% sapply(res_trace_i, function(x) x$Ftrgt))) {
      res_i <- res_trace_i[[which(Ftrgt == sapply(res_trace_i,
                                                  function(x) x$Ftrgt))[[1]]]]
      catch_i <- res_i$catch
      ssb_i <- res_i$ssb
      tsb_i <- res_i$tsb
      rec_i <- res_i$rec
    } else {
      ### run projection
      input$control$hcr@args$ftrg <- Ftrgt
      res_i <- do.call(mp, input)
      catch_i <- median(tail(catch(res_i@om@stock), 10), na.rm = TRUE)
      ssb_i <- median(tail(ssb(res_i@om@stock), 10), na.rm = TRUE)
      tsb_i <- median(tail(tsb(res_i@om@stock), 10), na.rm = TRUE)
      rec_i <- median(tail(rec(res_i@om@stock), 10), na.rm = TRUE)
    }
    ### print results of current run
    cat(paste0("Ftrgt=", Ftrgt, "; C=", catch_i, "; SSB=", ssb_i, "; R=", rec_i,
               "\n"))
    ### save results in res_trace
    res_add <- list(list(Ftrgt = Ftrgt, catch = catch_i,
                         ssb = ssb_i, tsb = tsb_i, rec = rec_i))
    assign(value = append(get("res_trace", envir = trace_env), res_add),
           x = "res_trace", envir = trace_env)
    if (isTRUE(minimise)) catch_i <- -catch_i
    return(catch_i)
  }

  ### first, check some values
  names(vals_ini) <- vals_ini
  res_ini <- lapply(vals_ini, mp_catch, input = input)

  ### use optimise - 1D golden-section search
  res_optimise_MSY <- optimise(f = mp_catch, input = input,
                               interval = c(lower, upper),
                               lower = lower, upper = upper,
                               maximum = TRUE,
                               tol = tol)

  ### get results and format
  res_trace_list <- unique(get("res_trace", envir = trace_env))
  res_trace <- as.data.frame(do.call(rbind, res_trace_list))
  res_trace <- as.data.frame(apply(res_trace, 2, unlist))
  res_trace <- unique(res_trace)

  if (isTRUE(save)) {
    saveRDS(res_trace_list, file = paste0(path, scenario, "_MSY_trace.rds"))
    write.csv(res_trace, file = paste0(path, scenario, "_MSY_trace.csv"),
              row.names = FALSE)
  }

  if (isTRUE(plot)) {
    ### plot
    p <- res_trace %>%
      dplyr::select(Ftrgt, catch, ssb, rec) %>%
      pivot_longer(cols = c("catch", "ssb", "rec")) %>%
      mutate(value = value/1000,
             name = factor(name, levels = c("catch", "ssb", "rec"),
                           labels = c("Catch (1000t)", "SSB (1000t)",
                                      "Recruitment (millions)"))) %>%
      ggplot(aes(x = Ftrgt, y = value)) +
      geom_point(size = 0.8) +
      stat_smooth(aes(alpha = "Loess smoother"), linewidth = 0.5,
                  se = FALSE, span = 0.3, n = 100, show.legend = TRUE) +
      scale_alpha_manual("", values = 1) +
      facet_wrap(~ name, scales = "free_y", strip.position = "left") +
      labs(x = x_label, y = "") +
      ylim(c(0, NA)) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      theme_bw() +
      theme(legend.position = "inside",
            legend.position.inside = c(0.85, 0.2),
            legend.background = element_blank(),
            legend.key = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.title.y = element_blank())
    ggsave(paste0(path, scenario,  "_MSY_search.png"), plot = p,
           width = 17, height = 6, units = "cm", dpi = 300, type = "cairo")
    ggsave(paste0(path, scenario, "_MSY_search.pdf"), plot = p,
           width = 17, height = 6, units = "cm")
  } else {
    p <- NULL
  }

  return(list(plot = p, result = res_trace))

}

### ------------------------------------------------------------------------ ###
### custom recruit-per-spawner recruitment model - NOT USED ####
### ------------------------------------------------------------------------ ###
### custom recruitment function: R ~ a*SSB
rps <- function(){

  ### log-likelihood for fitting model
  logl <- function(a, rec, ssb){
    loglAR1(log(rec),
            FLQuant(log(a * ssb),
                    dimnames = dimnames(ssb)))
  }

  ### recruitment model
  model <- rec ~ a * ssb

  ### initial parameter values
  initial <- structure(
    function(rec, ssb){
      return(
        FLPar(a = mean(c(rec)/c(ssb), na.rm = TRUE))
      )
    },
    lower = rep(0),
    upper = rep(Inf)
  )

  return(list(logl = logl, model = model, initial = initial))
}

### ------------------------------------------------------------------------ ###
### check lag-1 auto-correlation of time series ####
### ------------------------------------------------------------------------ ###
check_AR1 <- function(res, plot = FALSE, round = 4) {
  res_acf <- acf(res, plot = FALSE,
                 na.action = na.exclude)
  if (isTRUE(plot)) plot(res_acf)
  res_rho <- res_acf$acf[2]
  ci <- qnorm((1 + 0.95)/2)/sqrt(res_acf$n.used)
  paste0("lag-1 auto-correlation: ", round(res_rho, round), ";  ",
         ifelse(abs(res_rho) >= ci,
                "statistically significant",
                "NOT statistically significant"))
}

### ------------------------------------------------------------------------ ###
### generate residuals by fitting kernel density smoother ####
### ------------------------------------------------------------------------ ###
dens_resids <- function(res,
                        template,
                        res_ar1 = FALSE, ### include AR1?
                        res_ar1_log = FALSE,
                        by_iter = TRUE,
                        density.bw = "nrd0", ### default
                        iter_seed = TRUE,
                        seed,
                        cut_2SD = FALSE ### cut off beyond +-2SD
) {


  ### generate residuals by iteration
  if (isTRUE(by_iter)) {

    out <- lapply(seq(dim(template)[6]), function(iter_i) {
      ### get residuals for current iteration
      if (is.vector(res)) {
        res_i <- res
      } else {
        res_i <- c(FLCore::iter(res, iter_i))
      }
      res_i <- res_i[!is.na(res_i)]
      ### number of new residuals needed (number of year)
      n_res <- dim(template)[2]
      ### fit kernel density smoother to residuals
      dens <- density(res_i, bw = density.bw)
      ### set random number seet for reproducibility
      if (isTRUE(iter_seed)) set.seed(iter_i)
      ### sample residuals
      mu <- sample(x = res_i, size = n_res, replace = TRUE)
      ### "smooth" resampled residuals
      ### => equivalent to sampling from density distribution
      res_new <- rnorm(n = n_res, mean = mu, sd = dens$bw)

      if (isTRUE(cut_2SD)) {
        res_i_mean <- mean(res_i, na.rm = TRUE)
        res_i_sd <- sd(res_i, na.rm = TRUE)
        res_new[res_new < (res_i_mean - 2*res_i_sd)] <-
          res_i_mean - 2*res_i_sd
        res_new[res_new > (res_i_mean + 2*res_i_sd)] <-
          res_i_mean + 2*res_i_sd
      }
      ### include autocorrelation
      # if (isTRUE(res_ar1)) {
      #   ### use iteration-specific auto-correlation
      #   sr_acf_i <- acf(res_i, lag.max = 1, plot = FALSE,
      #                   na.action = na.exclude)
      #   sr_rho_i <- sr_acf_i$acf[2]
      #   res_ac <- rep(0, n_res)
      #   if (isTRUE(res_ar1_log)) {
      #     res_ac[1] <- exp(sr_rho_i * log(tail(res_i, 1)) +
      #       sqrt(1 - sr_rho_i^2) * log(res_new[1]))
      #     for (r in 2:length(res_ac)) {
      #       res_ac[r] <- exp(sr_rho_i * log(res_ac[r - 1]) +
      #         sqrt(1 - sr_rho_i^2) * log(res_new[r]))
      #     }
      #   } else {
      #     res_ac[1] <- sr_rho_i * log(tail(res_i, 1)) +
      #       sqrt(1 - sr_rho_i^2) * log(res_new[1])
      #     for (r in 2:length(res_ac)) {
      #       res_ac[r] <- sr_rho_i * res_ac[r - 1] +
      #         sqrt(1 - sr_rho_i^2) * res_new[r]
      #     }
      #   }
      #   res_new <- res_ac
      # }

      return(res_new)
    })
    template[] <- unlist(out)

    ### do not consider iterations (faster)
  } else {
    set.seed(seed)
    ### fit kernel density smoother to residuals
    dens <- density(res, bw = density.bw)
    ### sample residuals
    mu <- sample(x = res, size = prod(dim(template)),
                 replace = TRUE)
    ### "smooth" resampled residuals
    ### => equivalent to sampling from density distribution
    template[] <- rnorm(n = prod(dim(template)),
                        mean = mu, sd = dens$bw)
  }

  return(template)
}

### ------------------------------------------------------------------------ ###
### randNums ####
### ------------------------------------------------------------------------ ###
# function that generates n random integers with fixed sum s within boundaries lb and hb

randNums <- function(n,minVal,maxVal,s){
  hit <- FALSE
  while(hit==FALSE){
    total <- 0
    count <- 1
    nums <- rep(NA, n)
    
    while(total < s && count < n){
      r <- 0
      while(r<1){
      r <- round(runif(1, min=minVal-0.51, max=maxVal+0.49))
      }
      total <- total + r
      
      nums[count] <- r
      count <- count + 1
    }
    #nums[count:n] <- 0
    if(total == s){
      hit <- TRUE
    }
  }
  return(nums)
}

### ------------------------------------------------------------------------ ###
### randBlocks ####
### ------------------------------------------------------------------------ ###
randBlocks <- function(fecYears,projPeriod,nits){
  
  fecYears <- c(fecYears,max(fecYears)+1)
  
  nyrs <- length(projPeriod)
  
  saveBlcks_nYears     <- matrix(NA,nrow=nits,ncol=nyrs)
  saveBlcks_startYear  <- matrix(NA,nrow=nits,ncol=nyrs)
  saveBlcks_rev        <- matrix(NA,nrow=nits,ncol=nyrs)
  
  for(idxIter in 1:nits){
    saveBlcks_nYears[idxIter,] <- randNums(nyrs,1,length(fecYears)-1,nyrs)
  }
  
  # randomize the starting years
  for(idxIter in 1:nits){
    nYearsTemp <- saveBlcks_nYears[idxIter,
                                   !is.na(saveBlcks_nYears[idxIter,])] # select only non NaN values
    # loop on non NaN values
    for(idxYear in 1:length(nYearsTemp)){
      
      saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
                                                          min=min(an(fecYears)), 
                                                          max=(max(an(fecYears))-nYearsTemp[idxYear]))) # random number for the starting year
      
      #if(round(runif(1,0,1)) == 0){
      #  saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
      #                                                      min=min(an(fecYears)), 
      #                                                      max=(max(an(fecYears))-nYearsTemp[idxYear]))) # random number for the starting year
      #  saveBlcks_rev[idxIter,idxYear] <- 0
      #}else{
      #  saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
      #                                                      min=(min(an(fecYears))+nYearsTemp[idxYear]), 
      #                                                      max=max(an(fecYears)))) # random number for the starting year
      #  saveBlcks_rev[idxIter,idxYear] <- 1
      #}
    }
  }
  
  # fill in 2 arrays with start and end year for each block
  saveBlcksYrs <- array(NA,dim=c(nits,nyrs,2),dimnames=list(nits=1:nits,blocks=1:nyrs,strtstp=c("start","stop")))
  
  # simply loop on start years and number of years with a randomizer for a possible reverse order 
  for(idxIter in 1:nits){
    for(idxYear in 1:nyrs){
      if(!is.na(saveBlcks_startYear[idxIter,idxYear])){
        #if(saveBlcks_rev[idxIter,idxYear] == 1){ # reverse
        #  tempVal <- c( saveBlcks_startYear[idxIter,idxYear],
        #                saveBlcks_startYear[idxIter,idxYear] - saveBlcks_nYears[idxIter,idxYear])
        #  tempVal[2] <- tempVal[2]+1
        #}else{ # normal order
        #  tempVal <- c(saveBlcks_startYear[idxIter,idxYear],
        #               saveBlcks_startYear[idxIter,idxYear] + saveBlcks_nYears[idxIter,idxYear])
        #  tempVal[2] <- tempVal[2]-1
        #}
        
        tempVal <- c(saveBlcks_startYear[idxIter,idxYear],
                     saveBlcks_startYear[idxIter,idxYear] + saveBlcks_nYears[idxIter,idxYear])
        tempVal[2] <- tempVal[2]-1
        
        saveBlcksYrs[idxIter,idxYear,"start"] <- tempVal[1]
        saveBlcksYrs[idxIter,idxYear,"stop"]  <- tempVal[2]
      }
    }
  }
  
  yrStrngsIter <- apply(saveBlcksYrs,1,function(x){mapply(seq,from=na.omit(x[,1]),to=na.omit(x[,2]))})
  
  # if each sampled block was the same length, it creates a matrix instead of a list, gives an error below.
  # so added fix to change any matrices into lists first
  for (i in 1:nits) {
    if (!is.list(yrStrngsIter[[i]]))  {
      tmp <- list()
      for (cC in 1:ncol(yrStrngsIter[[i]])) tmp[[cC]]<- yrStrngsIter[[i]][,cC]
      yrStrngsIter[[i]] <- tmp
      rm(tmp)
    } 
  }
  
  yrStrngs    <- lapply(yrStrngsIter,function(x){do.call(c,x)})
  
  for (i in 1:nits) {
    yrStrngs[[i]] <- yrStrngs[[i]][1:length(projPeriod)]
  }
  
  return(yrStrngs)
}


