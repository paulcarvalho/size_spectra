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



  
  
  