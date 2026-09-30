### ------------------------------------------------------------------------ ###
### projection - fwd with process error ####
### ------------------------------------------------------------------------ ###
fwd_procres <- function(om, ctrl,
                        maxF = 5, ### maximum allowed Fbar
                        proc_res = NULL, ### process error noise,
                        ...) {

  stk <- stock(om)
  sr <- sr(om)
  sr.residuals <- residuals(sr(om))

  ### project forward with FLasher::fwd
  stk[] <- fwd(object = stk, control = ctrl, sr = sr,
               deviances = sr.residuals,
               maxF = maxF)

  ### add process error noise, if supplied
  if (!is.null(proc_res)) {

    ### projected years
    yrs_new <- seq(from = ctrl@target$year, to = range(stk)[["maxyear"]])

    if (!is(proc_res, "FLQuant"))
      stop("Process error needs to be provided in an FLQuant object")
    yrs_available <- as.numeric(dimnames(proc_res)$year)
    if (!all(yrs_new %in% yrs_available))
      stop("Process error years missing")

    ### check iterations
    ### may differ when parallelised -> subset
    if (isTRUE(dim(proc_res)[6] > dim(stk)[6])) {
      proc_res <- proc_res[,,,,, dimnames(stk)$iter]
    }

    ### implement process error
    stock.n(stk)[, ac(yrs_new)] <- stock.n(stk)[, ac(yrs_new)] *
      proc_res[, ac(yrs_new)]
    ### update stock biomass
    stock(stk) <- computeStock(stk)

  }

  ### return stock
  stock(om) <- stk
  return(list(om = om))

}

### ------------------------------------------------------------------------ ###
### OEM: shortcut observations ####
### ------------------------------------------------------------------------ ###
### return stock numbers etc without error
oem_shortcut <- function(stk, observations, deviances, tracking, args, ...){

  ### use OEM stock - biological data etc may differ from OM stock
  stk0 <- observations$stk

  ### update/insert fishery data without error
  catch.n(stk0) <- catch.n(stk)
  catch(stk0) <- catch(stk)
  discards.n(stk0) <- discards.n(stk)
  discards(stk0) <- discards(stk)
  landings.n(stk0) <- landings.n(stk)
  landings(stk0) <- landings(stk)

  ### shortcut - pass on N and F as observations without error
  stock.n(stk0) <- stock.n(stk)
  stock(stk0) <- stock(stk)
  harvest(stk0) <- harvest(stk)

  ### cut off observed stock to last data year
  stk0 <- window(stk0, end = args$dy)

  return(list(stk = stk0, idx = observations$idx, observations = observations,
              tracking = tracking))

}

### ------------------------------------------------------------------------ ###
### estimator emulator: SAM shortcut ####
### ------------------------------------------------------------------------ ###
SAM_emulator <- function(stk, idx, tracking,
                         args, ### contains ay (assessment year)
                         SSB_dev,
                         forecast = TRUE,
                         fwd_trgt = c("advice", "advice"), ### what to target in forecast
                         fwd_yrs = 2, ### number of years to add
                         fwd_yrs_average = 3, ### years used for averages
                         fwd_yrs_rec_start = 2014, ### recruitment
                         ...){

  ### no assessment
  stk0 <- stk

  ### implement SSB shortcut error - apply stock numbers at age
  ### entire stock scaled up/down
  if (!missing(SSB_dev)) {

    ### get years and iterations
    yrs_i <- dimnames(stk0)$year
    iter_i <- dimnames(stk0)$iter

    ### reshape SSB error into same format as stock numbers
    SSB_dev_i <- stock.n(stk0) %=% rep(c(SSB_dev[, yrs_i]),
                                       each = prod(dim(stk0)[1]))

    stock.n(stk0) <- stock.n(stk0) * SSB_dev_i

  }

  ### perform forecast to get SSB ay+1
  if (isTRUE(forecast)) {

    ### extend stock
    stk_stf <- stf(stk0, nyears = fwd_yrs, wts.nyears = fwd_yrs_average)

    ### recruitment: median to mimic average of SAM resampling
    rec_fc <- yearMedians(window(rec(stk0), start = fwd_yrs_rec_start))

    ### define STF target(s)
    if (isFALSE(identical(fwd_trgt, c("advice", "advice"))))
      stop("Unknown STF target(s)")

    ### 1st year (intermediate year) - catch
    if (identical(args$ay, args$iy)) {
      ### in first year of projection - extract from OM
      trgt1 <- tracking$data[tracking$metric == "C.om" & tracking$year == ac(args$ay)]
    } else {
      ### in future years, use catch advice
      trgt1 <- tracking$data[tracking$metric == "isys" &
                               tracking$year == ac(args$ay - args$management_lag)]
    }

    ### 2nd year - use same target (irrelevant)

    ### STF control object
    ctrl_stf <- fwdControl(list(year = c(args$ay, args$ay + 1),
                                quant = "catch",
                                value = rep(trgt1,  each = 2)))
    ### STF
    stk_stf[] <- fwd(stk_stf, control = ctrl_stf,
                     sr = FLSR(model = "geomean",
                               params = FLPar(a = c(rec_fc),
                                              iter = dims(stk_stf)$iter)),
                     verbose = FALSE)

    #print(args$ay)
    #browser()

    # fwd_warnings <- NULL
    # tc <- withCallingHandlers(
    #   expr = {
    #     stk_stf[] <- fwd(stk_stf, control = ctrl_stf,
    #                      sr = FLSR(model = "geomean",
    #                                params = FLPar(a = c(rec_fc),
    #                                               iter = dims(stk_stf)$iter)),
    #                      verbose = TRUE)
    #   },
    #   warning = function(w) {
    #     fwd_warnings <<- w$message
    #     invokeRestart("muffleWarning")
    #   }
    # )
    # ### rerun
    # if (isTRUE(grepl(x = fwd_warnings,
    #                  pattern = "^Unsolved targets at control rows*"))) {
    #   i <- 1
    #   browser()
    #   while (i <= 10) {
    #     fwd_warnings <- NULL
    #     tc <- withCallingHandlers(
    #       expr = {
    #         stk_stf[] <- fwd(stk_stf, control = ctrl_stf,
    #                          sr = FLSR(model = "geomean",
    #                                    params = FLPar(a = c(rec_fc),
    #                                                   iter = dims(stk_stf)$iter)),
    #                          verbose = TRUE)
    #       },
    #       warning = function(w) {
    #         fwd_warnings <<- w$message
    #         invokeRestart("muffleWarning")
    #       }
    #     )
    #     ### check if
    #     if ((!isTRUE(grepl(x = fwd_warnings,
    #                       pattern = "^Unsolved targets at control rows*")))) {
    #       break
    #     } else {
    #       i <- i + 1
    #     }
    #     print(paste0("rerun"))
    #   }
    #
    # }
    # fwd_warnings <- NULL

    stk0 <- stk_stf

  }

  ### return assessed stock
  return(list(stk = stk0, tracking = tracking))

}

### ------------------------------------------------------------------------ ###
### hockey stick HCR ####
### ------------------------------------------------------------------------ ###
hcr_hockeystick <- function(stk, ind, tracking, args,
                            Ftarget, Btrigger, Bto0 = 0, Bdrop = 0,
                            Fmin = 1e-6, Fmin0 = 1e-6, ...) {

  ### get SSB in advice year
  SSB <- c(ssb(stk)[, ac(args$ay + args$management_lag)])

  ### template for F multiplier
  Fmult <- rep(1, length(SSB))

  ### find SSB below Btrigger and reduce linearly to 0 at Bto0
  pos_below_Btrigger <- which(SSB < Btrigger)
  if (length(pos_below_Btrigger) > 0)
    Fmult[pos_below_Btrigger] <- (SSB[pos_below_Btrigger] - Bto0)/(Btrigger - Bto0)

  ### find SSB below Bdrop and set to 0
  pos_below_Bdrop <- which(SSB < Bdrop)
  if (length(pos_below_Bdrop) > 0)
    Fmult[pos_below_Bdrop] <- 0

  ### calculate resulting F
  Fout <- Ftarget * Fmult

  ### set F minimum
  pos_F_below_Fmin <- which(Fout < Fmin)
  if (length(pos_F_below_Fmin) > 0)
    Fout[Fout < Fmin] <- Fmin

  ### avoid 0 F target - can cause issues with fwd() later
  Fout[Fout < Fmin0] <- Fmin0

  ### create control object
  advice <- FLQuant(Fout,
                    dimnames = list(year = args$ay + args$management_lag,
                                    iter = seq(dim(stk)[6])))
  ctrl <- fwdControl(target = advice, quant = "f")


  return(list(ctrl = ctrl, tracking = tracking))

}

### ------------------------------------------------------------------------ ###
### implementation system for shortcut (including shortcut error on F) ####
### ------------------------------------------------------------------------ ###
isys_shortcut <- function(stk, args, ctrl, tracking,
                          F_dev,
                          fwd_yrs_average = 3, fwd_yrs_rec_start = 2014,
                          ...) {

  ### add shortcut error to target F
  if (!missing(F_dev))
    ctrl@iters[, "value", ] <- ctrl@iters[, "value", ] *
      c(F_dev[, ac(args$ay + args$management_lag)])

  ### recruitment for STF
  ### value already set in estimator - recycle
  rec_fc <- c(rec(stk)[, ac(args$ay + args$management_lag)])

  ### STF
  stk_stf <- fwd(stk, control = ctrl,
                 sr = FLSR(model = "geomean",
                           params = FLPar(a = c(rec_fc),
                                          iter = dims(stk)$iter)))

  ### extract catch
  catch_target <- catch(stk_stf)[, ac(args$ay + args$management_lag)]
  ### format into control object
  ctrl <- fwdControl(target = catch_target, quant = "catch")

  return(list(ctrl = ctrl, tracking = tracking))

}

### ------------------------------------------------------------------------ ###
### IEM for shortcut ####
### ------------------------------------------------------------------------ ###
iem_shortcut <- function(ctrl, args, tracking, iem_dev, ...) {

  if (!missing(iem_dev)) {

    ctrl@iters[, "value", ] <- ctrl@iters[, "value", ] *
      iem_dev[, ac(ctrl@target$year)]

  }

  return(list(ctrl = ctrl, tracking = tracking))

}

