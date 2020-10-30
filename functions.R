# functions for size spectra analysis from https://github.com/andrew-edwards/fitting-size-spectra

# Statistical functions:
# set.params - create list with parameters for running MLE
# negLL.PLB - negative log-likelihood function for PLB model
# pPLB - probability distribution function for bounded power-law distribution
# mle_b - # Use analytical value of MLE b for PL model
# slope.conf.int - calculate confidence intervals of estimated slope

# Plotting functions:
# logTicks - add axes and tick marks to a log-log plot to represent
# trophic.plots - plot fitted slope for different trophic positions

set.params = function(biomass_kg){
  # Creates a list of parameters for input into the negll.PLB and pPLB functions
  biomass <- biomass_kg
  log.biomass <- log(biomass)
  sum.log.biomass <- sum(log.biomass)
  min.biomass <- min(biomass_kg)
  max.biomass <- max(biomass_kg)
  out.list <- (list(biomass, log.biomass, sum.log.biomass, min.biomass, max.biomass))
  names(out.list) <- c("biomass", "log.biomass", "sum.log.biomass", "min.biomass", "max.biomass")
  return(out.list)
}

negLL.PLB = function(b, x, n, xmin, xmax, sumlogx)
  {
  # Calculates the negative log-likelihood of the parameters b, xmin and xmax
  #  given data x for the PLB model. Returns the negative log-likelihood. Will
  #  be called by nlm or similar, but xmin and xmax are just estimated as the
  #  min and max of the data, not numerically using likelihood.
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   x: vector of values of data (e.g. masses of individual fish)
  #   n: length(x), have as an input to avoid repeatedly calculating it 
  #   xmin: minimum value of x, have as an input to avoid repeatedly calculating
  #   xmax: maximum value of x, have as an input to avoid repeatedly calculating
  #   sumlogx: sum(log(x)) as an input, to avoid repeatedly calculating
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
}

pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100)    
  {
  # Computes probability distribution function, P(X <= x),  for a
  #   bounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x     # so have zeros where x < xmin
    y[x > xmax] = 1  # 1 for x > xmax
    if(b != -1)
        {  xmintobplus1 = xmin^(b+1)
           denom = xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] =
               ( x[x >= xmin & x <= xmax]^(b + 1) - xmintobplus1 ) / denom
        } else
        {  logxmin = log(xmin)
           denom = log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)                              
  }

mle_b = function(region, x, log_x, sum_log_x, x_min, x_max){ # function added by PC
    # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
    # as a starting point for nlm for MLE of b for PLB model.

    PL.bMLE = 1/( log(min(x)) - sum_log_x/length(x)) - 1
    
    PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
        xmin=x_min, xmax=x_max, sumlogx=sum_log_x, hessian = TRUE) #, print.level=2 )
    
    PLB.bMLE = PLB.minLL$estimate
    
    PLB.return = list(PLB.bMLE, PLB.minLL)
    
    return(PLB.return)
}

slope.conf.int = function(PLB.bMLE.b, PLB.minLL.b, input){
  # Confidence intervals
  bvec = seq(PLB.bMLE.b - 0.5, PLB.bMLE.b + 0.5, 0.00001) 
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec)){
    PLB.LLvals[i] = negLL.PLB(bvec[i], x=input$biomass, n=length(input$biomass), xmin=input$min.biomass,
        xmax=input$max.biomass, sumlogx=input$sum.log.biomass)   
  }
  critVal = PLB.minLL.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
  return(c(min(bIn95), max(bIn95)))
}

logTicks = function(xLim, yLim = NULL, tclSmall = -0.2, xLabelSmall = NULL,
                    yLabelSmall = NULL, xLabelBig = NULL, mgpVal=c(1.6,0.5,0)){
    # Add axes and tick marks to a log-log plot to represent unlogged values.
    # Args:
    #  xLim: the x limits for the plot (unlogged scale); if NULL then
    #         do not add anything to x-axis
    #  yLim: the y limits for the plot (unlogged scale); if NULL then
    #         do not add anything to y-axis
    #  tclSmall: size of small tick marks
    #  xLabelSmall: which small tick marks on x-axis to label
    #  yLabelSmall: which small tick marks on y-axis to label
    #  xLabelBig: which big tick marks on the x-axis to label
    #   (when automated they can overlap, so may need to specify).
    #  mgpVal: mgp values for axes. See ?par 
    # Returns:
    #  Adds axes and big and small tick marks to the plot. Returns NULL.
    # Example:
    #  Adapt the following:
    #  plot(..., log="xy", xlab=..., ylab=..., xlim=..., ylim=..., axes=FALSE)
    #  xLim = 10^par("usr")[1:2]
    #  yLim = 10^par("usr")[3:4]
    #  logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500)) 
    # 
    ll = 1:9
    log10ll = log10(ll)
    # box()                                                                       # PC removed
    # x axis
    if(!is.null(xLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      xEncompassLog = c(floor(log10(xLim[1])), ceiling(log10(xLim[2])))
      xBig = 10^c(xEncompassLog[1]:xEncompassLog[2])
      # Big unlabelled, always want these:
      # axis(1, at= xBig, labels = rep("", length(xBig)), mgp = mgpVal)           # PC removed
      # Big labelled:
      if(is.null(xLabelBig)) { xLabelBig = xBig }
      axis(1, at= xLabelBig, labels = xLabelBig, mgp = mgpVal)
      # axis(1, at=c(1, 10, 100), labels = c(1, 10, 100), mgp=c(1.7,0.7,0))
      # Small unlabelled:
      # axis(1, xBig %x% ll, labels=rep("", length(xBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(xLabelSmall))
          {
          axis(1, at=xLabelSmall, labels=xLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }
    # Repeat for y axis:
    if(!is.null(yLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      yEncompassLog = c(floor(log10(yLim[1])), ceiling(log10(yLim[2])))
      # yBig = 10^c(yEncompassLog[1]:yEncompassLog[2])                            # PC removed
      yBig = c(0,1,10,100,1000,10000)
      # Big labelled:
      axis(2, at= yBig, labels = yBig, mgp = mgpVal)
      # Small unlabelled:
      # axis(2, yBig %x% ll, labels=rep("", length(yBig %x% ll)), tcl=tclSmall)   # PC removed
      # Small labelled:
      if(!is.null(yLabelSmall))
          {
          axis(2, at=yLabelSmall, labels=yLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }     
}

trophic.plots = function(PLB.return, PLB.bMLE.b, PLB.minLL.b, df.input, mgpVals, troph_id, panel){
  # plot and find 95% confidence intervals for MLE method.
  PLB.minNegLL.b <- PLB.minLL.b$minimum
  x <- df.input$biomass  
  plot(sort(df.input$biomass, decreasing=TRUE), 1:length(df.input$biomass), log="xy",
    xlab=expression(paste("Body sizes, ", italic(x), " (log[kg])")),
    # xlab = "log(kg)",
    ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
    # ylab = "  ",
    xlim = c(df.input$min.biomass, df.input$max.biomass), ylim = c(1, 10000), axes=FALSE)
  logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10)) # Tick marks.
  x.PLB = seq(min(df.input$biomass), max(df.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.b, xmin = min(x.PLB),
    xmax = max(x.PLB))) * length(df.input$biomass)
  lines(x.PLB, y.PLB, col="red", lwd=2)
  text(x = 0.1, y = 15, labels = troph_id, cex = 1.1, pos = 1, col = "black")
  spectra.text <- as.character(round(PLB.bMLE.b, 2))
  text(x=0.1, y=5, labels = bquote(paste(italic("b = "), .(spectra.text))), 
    cex=1.1, pos=1, col="black")
  mtext(panel, side = 3, cex = 1.4, adj = -0.15)
  # Values of b to test to obtain confidence interval. For the real movement data
  # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  # symmetric interval here.
  bvec = seq(PLB.bMLE.b - 0.5, PLB.bMLE.b + 0.5, 0.00001) 
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec)){
    PLB.LLvals[i] = negLL.PLB(bvec[i], x=df.input$biomass, n=length(df.input$biomass), xmin=df.input$min.biomass,
      xmax=df.input$max.biomass, sumlogx=df.input$sum.log.biomass)   
  }
  critVal = PLB.minNegLL.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
  # To add just the curves at the limits of the 95% confidence interval of b:
  for(i in c(1, length(bIn95))){
    lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
      xmax = max(x.PLB))) * length(df.input$biomass), col="red", lty=2)
  }
}

slope_regAndFunc = function(carn.df, herb.df){
  # Paul - Carnivore, Raja Ampat
  pc.carn.ra.df <- carn.df %>%
    filter(observer == "PC") %>%
    filter(region == "raja_ampat")
  pc.carn.ra.input <- set.params(pc.carn.ra.df$biomass_kg)
  PLB.return.pc.carn.ra <- mle_b(region=NA, x=pc.carn.ra.input$biomass, log_x=pc.carn.ra.input$log.biomass, sum_log_x=pc.carn.ra.input$sum.log.biomass,
                   x_min=pc.carn.ra.input$min.biomass, x_max=pc.carn.ra.input$max.biomass)
  PLB.bMLE.pc.carn.ra.b <- PLB.return.pc.carn.ra[[1]] 
  PLB.minLL.pc.carn.ra.b <- PLB.return.pc.carn.ra[[2]]
  pc.carn.ra.b.In95 <- slope.conf.int(PLB.bMLE.pc.carn.ra.b, PLB.minLL.pc.carn.ra.b$minimum, pc.carn.ra.input)
  divers.tp.df <- data.frame(obs_region = "PC_RA",
                             tp = "Carnivore",
                             b = PLB.bMLE.pc.carn.ra.b,
                             b_lwr = pc.carn.ra.b.In95[1],
                             b_upr = pc.carn.ra.b.In95[2])
  
  # Paul - Herbivore, Raja Ampat
  pc.herb.ra.df <- herb.df %>%
    filter(observer == "PC") %>%
    filter(region == "raja_ampat")
  pc.herb.ra.input <- set.params(pc.herb.ra.df$biomass_kg)
  PLB.return.pc.herb.ra <- mle_b(region=NA, x=pc.herb.ra.input$biomass, log_x=pc.herb.ra.input$log.biomass, sum_log_x=pc.herb.ra.input$sum.log.biomass,
                   x_min=pc.herb.ra.input$min.biomass, x_max=pc.herb.ra.input$max.biomass)
  PLB.bMLE.pc.herb.ra.b <- PLB.return.pc.herb.ra[[1]] 
  PLB.minLL.pc.herb.ra.b <- PLB.return.pc.herb.ra[[2]]
  pc.herb.ra.b.In95 <- slope.conf.int(PLB.bMLE.pc.herb.ra.b, PLB.minLL.pc.herb.ra.b$minimum, pc.herb.ra.input)
  tmp <- data.frame(obs_region = "PC_RA",
                    tp = "Herbivore",
                    b = PLB.bMLE.pc.herb.ra.b,
                    b_lwr = pc.herb.ra.b.In95[1],
                    b_upr = pc.herb.ra.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Paul - Carnivore, Wakatobi
  pc.carn.wa.df <- carn.df %>%
    filter(observer == "PC") %>%
    filter(region == "wakatobi")
  pc.carn.wa.input <- set.params(pc.carn.wa.df$biomass_kg)
  PLB.return.pc.carn.wa <- mle_b(region=NA, x=pc.carn.wa.input$biomass, log_x=pc.carn.wa.input$log.biomass, sum_log_x=pc.carn.wa.input$sum.log.biomass,
                   x_min=pc.carn.wa.input$min.biomass, x_max=pc.carn.wa.input$max.biomass)
  PLB.bMLE.pc.carn.wa.b <- PLB.return.pc.carn.wa[[1]] 
  PLB.minLL.pc.carn.wa.b <- PLB.return.pc.carn.wa[[2]]
  pc.carn.wa.b.In95 <- slope.conf.int(PLB.bMLE.pc.carn.wa.b, PLB.minLL.pc.carn.wa.b$minimum, pc.carn.wa.input)
  tmp <- data.frame(obs_region = "PC_WA",
                    tp = "Carnivore",
                    b = PLB.bMLE.pc.carn.wa.b,
                    b_lwr = pc.carn.wa.b.In95[1],
                    b_upr = pc.carn.wa.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Paul - Herbivore, Wakatobi
  pc.herb.wa.df <- herb.df %>%
    filter(observer == "PC") %>%
    filter(region == "wakatobi")
  pc.herb.wa.input <- set.params(pc.herb.wa.df$biomass_kg)
  PLB.return.pc.herb.wa <- mle_b(region=NA, x=pc.herb.wa.input$biomass, log_x=pc.herb.wa.input$log.biomass, sum_log_x=pc.herb.wa.input$sum.log.biomass,
                   x_min=pc.herb.wa.input$min.biomass, x_max=pc.herb.wa.input$max.biomass)
  PLB.bMLE.pc.herb.wa.b <- PLB.return.pc.herb.wa[[1]] 
  PLB.minLL.pc.herb.wa.b <- PLB.return.pc.herb.wa[[2]]
  pc.herb.wa.b.In95 <- slope.conf.int(PLB.bMLE.pc.herb.wa.b, PLB.minLL.pc.herb.wa.b$minimum, pc.herb.wa.input)
  tmp <- data.frame(obs_region = "PC_WA",
                    tp = "Herbivore",
                    b = PLB.bMLE.pc.herb.wa.b,
                    b_lwr = pc.herb.wa.b.In95[1],
                    b_upr = pc.herb.wa.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Paul - Carnivore, Lombok
  pc.carn.lo.df <- carn.df %>%
    filter(observer == "PC") %>%
    filter(region == "lombok")
  pc.carn.lo.input <- set.params(pc.carn.lo.df$biomass_kg)
  PLB.return.pc.carn.lo <- mle_b(region=NA, x=pc.carn.lo.input$biomass, log_x=pc.carn.lo.input$log.biomass, sum_log_x=pc.carn.lo.input$sum.log.biomass,
                   x_min=pc.carn.lo.input$min.biomass, x_max=pc.carn.lo.input$max.biomass)
  PLB.bMLE.pc.carn.lo.b <- PLB.return.pc.carn.lo[[1]] 
  PLB.minLL.pc.carn.lo.b <- PLB.return.pc.carn.lo[[2]]
  pc.carn.lo.b.In95 <- slope.conf.int(PLB.bMLE.pc.carn.lo.b, PLB.minLL.pc.carn.lo.b$minimum, pc.carn.lo.input)
  tmp <- data.frame(obs_region = "PC_LO",
                    tp = "Carnivore",
                    b = PLB.bMLE.pc.carn.lo.b,
                    b_lwr = pc.carn.lo.b.In95[1],
                    b_upr = pc.carn.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Paul - Herbivore, Lombok
  pc.herb.lo.df <- herb.df %>%
    filter(observer == "PC") %>%
    filter(region == "lombok")
  pc.herb.lo.input <- set.params(pc.herb.lo.df$biomass_kg)
  PLB.return.pc.herb.lo <- mle_b(region=NA, x=pc.herb.lo.input$biomass, log_x=pc.herb.lo.input$log.biomass, sum_log_x=pc.herb.lo.input$sum.log.biomass,
                   x_min=pc.herb.lo.input$min.biomass, x_max=pc.herb.lo.input$max.biomass)
  PLB.bMLE.pc.herb.lo.b <- PLB.return.pc.herb.lo[[1]] 
  PLB.minLL.pc.herb.lo.b <- PLB.return.pc.herb.lo[[2]]
  pc.herb.lo.b.In95 <- slope.conf.int(PLB.bMLE.pc.herb.lo.b, PLB.minLL.pc.herb.lo.b$minimum, pc.herb.lo.input)
  tmp <- data.frame(obs_region = "PC_LO",
                    tp = "Herbivore",
                    b = PLB.bMLE.pc.herb.lo.b,
                    b_lwr = pc.herb.lo.b.In95[1],
                    b_upr = pc.herb.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Ubun - Carnivore, Raja Ampat
  fs.carn.ra.df <- carn.df %>%
    filter(observer == "FS") %>%
    filter(region == "raja_ampat")
  fs.carn.ra.input <- set.params(fs.carn.ra.df$biomass_kg)
  PLB.return.fs.carn.ra <- mle_b(region=NA, x=fs.carn.ra.input$biomass, log_x=fs.carn.ra.input$log.biomass, sum_log_x=fs.carn.ra.input$sum.log.biomass,
                   x_min=fs.carn.ra.input$min.biomass, x_max=fs.carn.ra.input$max.biomass)
  PLB.bMLE.fs.carn.ra.b <- PLB.return.fs.carn.ra[[1]] 
  PLB.minLL.fs.carn.ra.b <- PLB.return.fs.carn.ra[[2]]
  fs.carn.ra.b.In95 <- slope.conf.int(PLB.bMLE.fs.carn.ra.b, PLB.minLL.fs.carn.ra.b$minimum, fs.carn.ra.input)
  tmp <- data.frame(obs_region = "FS_RA",
                    tp = "Carnivore",
                    b = PLB.bMLE.fs.carn.ra.b,
                    b_lwr = fs.carn.ra.b.In95[1],
                    b_upr = fs.carn.ra.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Ubun - Herbivore, Raja Ampat
  fs.herb.ra.df <- herb.df %>%
    filter(observer == "FS") %>%
    filter(region == "raja_ampat")
  fs.herb.ra.input <- set.params(fs.herb.ra.df$biomass_kg)
  PLB.return.fs.herb.ra <- mle_b(region=NA, x=fs.herb.ra.input$biomass, log_x=fs.herb.ra.input$log.biomass, sum_log_x=fs.herb.ra.input$sum.log.biomass,
                   x_min=fs.herb.ra.input$min.biomass, x_max=fs.herb.ra.input$max.biomass)
  PLB.bMLE.fs.herb.ra.b <- PLB.return.fs.herb.ra[[1]] 
  PLB.minLL.fs.herb.ra.b <- PLB.return.fs.herb.ra[[2]]
  fs.herb.ra.b.In95 <- slope.conf.int(PLB.bMLE.fs.herb.ra.b, PLB.minLL.fs.herb.ra.b$minimum, fs.herb.ra.input)
  tmp <- data.frame(obs_region = "FS_RA",
                    tp = "Herbivore",
                    b = PLB.bMLE.fs.herb.ra.b,
                    b_lwr = fs.herb.ra.b.In95[1],
                    b_upr = fs.herb.ra.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Ubun - Carnivore, Lombok
  fs.carn.lo.df <- carn.df %>%
    filter(observer == "FS") %>%
    filter(region == "lombok")
  fs.carn.lo.input <- set.params(fs.carn.lo.df$biomass_kg)
  PLB.return.fs.carn.lo <- mle_b(region=NA, x=fs.carn.lo.input$biomass, log_x=fs.carn.lo.input$log.biomass, sum_log_x=fs.carn.lo.input$sum.log.biomass,
                   x_min=fs.carn.lo.input$min.biomass, x_max=fs.carn.lo.input$max.biomass)
  PLB.bMLE.fs.carn.lo.b <- PLB.return.fs.carn.lo[[1]] 
  PLB.minLL.fs.carn.lo.b <- PLB.return.fs.carn.lo[[2]]
  fs.carn.lo.b.In95 <- slope.conf.int(PLB.bMLE.fs.carn.lo.b, PLB.minLL.fs.carn.lo.b$minimum, fs.carn.lo.input)
  tmp <- data.frame(obs_region = "FS_LO",
                    tp = "Carnivore",
                    b = PLB.bMLE.fs.carn.lo.b,
                    b_lwr = fs.carn.lo.b.In95[1],
                    b_upr = fs.carn.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Ubun - Herbivore, Lombok
  fs.herb.lo.df <- herb.df %>%
    filter(observer == "FS") %>%
    filter(region == "lombok")
  fs.herb.lo.input <- set.params(fs.herb.lo.df$biomass_kg)
  PLB.return.fs.herb.lo <- mle_b(region=NA, x=fs.herb.lo.input$biomass, log_x=fs.herb.lo.input$log.biomass, sum_log_x=fs.herb.lo.input$sum.log.biomass,
                   x_min=fs.herb.lo.input$min.biomass, x_max=fs.herb.lo.input$max.biomass)
  PLB.bMLE.fs.herb.lo.b <- PLB.return.fs.herb.lo[[1]] 
  PLB.minLL.fs.herb.lo.b <- PLB.return.fs.herb.lo[[2]]
  fs.herb.lo.b.In95 <- slope.conf.int(PLB.bMLE.fs.herb.lo.b, PLB.minLL.fs.herb.lo.b$minimum, fs.herb.lo.input)
  tmp <- data.frame(obs_region = "FS_LO",
                    tp = "Herbivore",
                    b = PLB.bMLE.fs.herb.lo.b,
                    b_lwr = fs.herb.lo.b.In95[1],
                    b_upr = fs.herb.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Choro - Carnivore, Lombok
  ch.carn.lo.df <- carn.df %>%
    filter(observer == "ch") %>%
    filter(region == "lombok")
  ch.carn.lo.input <- set.params(ch.carn.lo.df$biomass_kg)
  PLB.return.ch.carn.lo <- mle_b(region=NA, x=ch.carn.lo.input$biomass, log_x=ch.carn.lo.input$log.biomass, sum_log_x=ch.carn.lo.input$sum.log.biomass,
                   x_min=ch.carn.lo.input$min.biomass, x_max=ch.carn.lo.input$max.biomass)
  PLB.bMLE.ch.carn.lo.b <- PLB.return.ch.carn.lo[[1]] 
  PLB.minLL.ch.carn.lo.b <- PLB.return.ch.carn.lo[[2]]
  ch.carn.lo.b.In95 <- slope.conf.int(PLB.bMLE.ch.carn.lo.b, PLB.minLL.ch.carn.lo.b$minimum, ch.carn.lo.input)
  tmp <- data.frame(obs_region = "CH_LO",
                    tp = "Carnivore",
                    b = PLB.bMLE.ch.carn.lo.b,
                    b_lwr = ch.carn.lo.b.In95[1],
                    b_upr = ch.carn.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  # Choro - Herbivore, Lombok
  ch.herb.lo.df <- herb.df %>%
    filter(observer == "ch") %>%
    filter(region == "lombok")
  ch.herb.lo.input <- set.params(ch.herb.lo.df$biomass_kg)
  PLB.return.ch.herb.lo <- mle_b(region=NA, x=ch.herb.lo.input$biomass, log_x=ch.herb.lo.input$log.biomass, sum_log_x=ch.herb.lo.input$sum.log.biomass,
                   x_min=ch.herb.lo.input$min.biomass, x_max=ch.herb.lo.input$max.biomass)
  PLB.bMLE.ch.herb.lo.b <- PLB.return.ch.herb.lo[[1]] 
  PLB.minLL.ch.herb.lo.b <- PLB.return.ch.herb.lo[[2]]
  ch.herb.lo.b.In95 <- slope.conf.int(PLB.bMLE.ch.herb.lo.b, PLB.minLL.ch.herb.lo.b$minimum, ch.herb.lo.input)
  tmp <- data.frame(obs_region = "CH_LO",
                    tp = "Herbivore",
                    b = PLB.bMLE.ch.herb.lo.b,
                    b_lwr = ch.herb.lo.b.In95[1],
                    b_upr = ch.herb.lo.b.In95[2])
  divers.tp.df <- rbind(divers.tp.df, tmp)
  
  divers.tp.df$region <- c("Raja Ampat", "Raja Ampat", "Wakatobi", "Wakatobi", "Lombok", "Lombok", "Raja Ampat", "Raja Ampat",
                           "Lombok", "Lombok", "Lombok", "Lombok")
                          
  divers.tp.df <- divers.tp.df %>%
    mutate(obs_region = fct_relevel(obs_region, "PC_RA","FS_RA","PC_WA","PC_LO","FS_LO","CH_LO")) 
  
  p <- ggplot() +
    geom_pointrange(data = divers.tp.df, aes(x = obs_region, y = b, color = tp, shape = region, ymin = b_lwr, ymax = b_upr)) +
    theme_classic() +
    labs(x = "Observer", y = expression(paste("Size spectra (", italic("b"),")"))) +
    scale_x_discrete(labels = c("PC", "FS", "PC", "PC", "FS", "PS")) +
    theme(legend.title = element_blank())
  
  return(p)
}

slope_reg = function(){
  # Paul - Raja Ampat
  pc.ra.df <- ra %>%
    filter(observer == "PC")
  pc.ra.input <- set.params(pc.ra.df$biomass_kg)
  PLB.return.pc.ra <- mle_b(region="raja_ampat", x=pc.ra.input$biomass, log_x=pc.ra.input$log.biomass, sum_log_x=pc.ra.input$sum.log.biomass,
                   x_min=pc.ra.input$min.biomass, x_max=pc.ra.input$max.biomass)
  PLB.bMLE.pc.ra.b <- PLB.return.pc.ra[[1]] 
  PLB.minLL.pc.ra.b <- PLB.return.pc.ra[[2]]
  pc.ra.b.In95 <- slope.conf.int(PLB.bMLE.pc.ra.b, PLB.minLL.pc.ra.b$minimum, pc.ra.input)
  
  # Ubun - Raja Ampat
  fs.ra.df <- ra %>%
    filter(observer == "FS")
  fs.ra.input <- set.params(fs.ra.df$biomass_kg)
  PLB.return.fs.ra <- mle_b(region="raja_ampat", x=fs.ra.input$biomass, log_x=fs.ra.input$log.biomass, sum_log_x=fs.ra.input$sum.log.biomass,
                   x_min=fs.ra.input$min.biomass, x_max=fs.ra.input$max.biomass)
  PLB.bMLE.fs.ra.b <- PLB.return.fs.ra[[1]] 
  PLB.minLL.fs.ra.b <- PLB.return.fs.ra[[2]]
  fs.ra.b.In95 <- slope.conf.int(PLB.bMLE.fs.ra.b, PLB.minLL.fs.ra.b$minimum, fs.ra.input)
  
  # Paul - Lombok
  pc.lo.df <- lo %>%
    filter(observer == "PC")
  pc.lo.input <- set.params(pc.lo.df$biomass_kg)
  PLB.return.pc.lo <- mle_b(region="lombok", x=pc.lo.input$biomass, log_x=pc.lo.input$log.biomass, sum_log_x=pc.lo.input$sum.log.biomass,
                   x_min=pc.lo.input$min.biomass, x_max=pc.lo.input$max.biomass)
  PLB.bMLE.pc.lo.b <- PLB.return.pc.lo[[1]] 
  PLB.minLL.pc.lo.b <- PLB.return.pc.lo[[2]]
  pc.lo.b.In95 <- slope.conf.int(PLB.bMLE.pc.lo.b, PLB.minLL.pc.lo.b$minimum, pc.lo.input)
  
  # Ubun - Lombok
  fs.lo.df <- lo %>%
    filter(observer == "FS")
  fs.lo.input <- set.params(fs.lo.df$biomass_kg)
  PLB.return.fs.lo <- mle_b(region="lombok", x=fs.lo.input$biomass, log_x=fs.lo.input$log.biomass, sum_log_x=fs.lo.input$sum.log.biomass,
                   x_min=fs.lo.input$min.biomass, x_max=fs.lo.input$max.biomass)
  PLB.bMLE.fs.lo.b <- PLB.return.fs.lo[[1]] 
  PLB.minLL.fs.lo.b <- PLB.return.fs.lo[[2]]
  fs.lo.b.In95 <- slope.conf.int(PLB.bMLE.fs.lo.b, PLB.minLL.fs.lo.b$minimum, fs.lo.input)
  
  # Choro - Lombok
  ch.lo.df <- lo %>%
    filter(observer == "ch")
  ch.lo.input <- set.params(ch.lo.df$biomass_kg)
  PLB.return.ch.lo <- mle_b(region="lombok", x=ch.lo.input$biomass, log_x=ch.lo.input$log.biomass, sum_log_x=ch.lo.input$sum.log.biomass,
                   x_min=ch.lo.input$min.biomass, x_max=ch.lo.input$max.biomass)
  PLB.bMLE.ch.lo.b <- PLB.return.ch.lo[[1]] 
  PLB.minLL.ch.lo.b <- PLB.return.ch.lo[[2]]
  ch.lo.b.In95 <- slope.conf.int(PLB.bMLE.ch.lo.b, PLB.minLL.ch.lo.b$minimum, ch.lo.input)
  
  divers.df <- data.frame(obs_region = c("PC_RA","FS_RA","PC_LO","FS_LO","PS_LO"),
                          region  = c("Raja Ampat", "Raja Ampat", "Lombok", "Lombok", "Lombok"),
                          b = c(PLB.bMLE.pc.ra.b, PLB.bMLE.fs.ra.b, PLB.bMLE.pc.lo.b, PLB.bMLE.fs.lo.b, PLB.bMLE.ch.lo.b),
                          b_lwr = c(pc.ra.b.In95[1], fs.ra.b.In95[1], pc.lo.b.In95[1], fs.lo.b.In95[1], ch.lo.b.In95[1]),
                          b_upr = c(pc.ra.b.In95[2], fs.ra.b.In95[2], pc.lo.b.In95[2], fs.lo.b.In95[2], ch.lo.b.In95[2]))
  
  divers.df <- divers.df %>%
    mutate(obs_region = fct_relevel(obs_region, "PC_RA","FS_RA","PC_LO","FS_LO","PS_LO"))
  
  p <- ggplot() +
    geom_pointrange(data = divers.df, aes(x = obs_region, y = b, ymin = b_lwr, ymax = b_upr, color = region)) +
    theme_classic() +
    labs(x = "Observer", y = "Size spectrum slope") +
    scale_x_discrete(labels = c("PC_RA" = "PC", "FS_RA" = "FS", "PC_LO" = "PC", "FS_LO" = "FS", "PS_LO" = "PS")) +
    theme(legend.title = element_blank())
  
  return(p)
}
  