###################### Project-specific Functions ############################

######## User functions ########
mk_datsummary = function(lsfnms = mk_lsfnms(), 
                         pngfnm = "datsummary.png", pngtitle = "", 
                         diagColor = NA, datUnits = "", bQuantile = TRUE,
                         fromRes = NA, toRes = NA, bGrid = FALSE,  
                         zrange = 3.0, cutoff = 1.0)
{
  datsummary(lsfnms, 5, 8, pngfnm = pngfnm, pngtitle = pngtitle, 
             rnames = c("M1", "M2toM1", "M2FBP", "M2", "M2H"), 
             cnames = c("R1 Chain A", "R1 Chain B", "R1 Chain C", "R1 Chain D", 
                        "R2 Chain A", "R2 Chain B", "R2 Chain C", "R2 Chain D"), 
             bQuantile = bQuantile, diagColor = diagColor, bColorBar = TRUE, 
             datUnits = datUnits, fromRes = fromRes, toRes = toRes, bGrid = bGrid, 
             zmean = 0.0, zrange = zrange, cutoff = cutoff, ospace = 6)
}


mk_datoutlier = function(lsfnms = mk_lsfnms(), fromRes = NA, toRes = NA, bGrid = FALSE, 
                         pngfnm = "datoutlier.png", pngtitle = "", 
                         diagColor = NA, datUnits = "", bQuantile = TRUE,  
                         zrange = 3.0, cutoff = 1.0)
{
  datoutlier(lsfnms, 5, 8, pngfnm = pngfnm, pngtitle = pngtitle, 
             rnames = c("M1", "M2toM1", "M2FBP", "M2", "M2H"), 
             cnames = c("R1 Chain A", "R1 Chain B", "R1 Chain C", "R1 Chain D", 
                        "R2 Chain A", "R2 Chain B", "R2 Chain C", "R2 Chain D"), 
             bQuantile = bQuantile, diagColor = diagColor, bColorBar = TRUE, 
             datUnits = datUnits, fromRes = fromRes, toRes = toRes, bGrid = bGrid,
             zmean = 0.0, zrange = zrange, cutoff = cutoff, ospace = 6)
}


mk_diffsummary = function(lsfnms = mk_lsfnms(), fromRes = NA, toRes = NA, 
                         pngfnm = "diffsummary.png", pngtitle = "", 
                         diagColor = 0.5, datUnits = "", bQuantile = TRUE, 
                         zrange = 3.0, cutoff = 0.5, bGrid = FALSE)
{
  fnms2diffsummary(lsfnms, pngfnm = pngfnm, pngtitle = pngtitle, 
             rcnames = c("M1", "M2toM1", "M2FBP", "M2", "M2H"),  
             bQuantile = bQuantile, diagColor = diagColor, bColorBar = TRUE, 
             datUnits = datUnits, fromRes = fromRes, toRes = toRes, bGrid = bGrid,
             zmean = 0.0, zrange = zrange, cutoff = cutoff, ospace = 6)
}


mk_diffsummary2 = function(rlsfnms = mk_lsfnms(includeFiles = c(1, 3, 4),
                                              addPattern = "rdat"), 
                           dlsfnms = mk_lsfnms(includeFiles = c(1, 3, 4),
                                              addPattern = "ddat"), 
                           bPNG = TRUE, avgDiffPNG = "",
                           pngfnm = "diffsummary.png", pngtitle = "", 
                           bColorBar = TRUE, bGrid = FALSE, 
                           fromRes = NA, toRes = NA, ospace = 6, 
                           pngWidth = NA, pngHeight = NA, fWidth = 500, 
                           fHeight = 500, bDiagNA = TRUE, diagColor = 0.5, 
                           bQuantile = TRUE, zmean = 0.0, zrange = 3.0, 
                           cutoff = 0.5, datUnits = "",
                           rcnames = c("M1", "M2FBP", "M2"))
{
  lsdats = list()
  for (i in 1:length(rlsfnms))
  {
    rdat = avg2png(rlsfnms[[i]], datUnits = datUnits, bQuantile = FALSE, 
                   bDiagNA = bDiagNA, diagColor = diagColor, 
                   zmean = 0.0, zrange = zrange, cutoff = cutoff)
    rdat = r2sigma(rdat)
    ddat = avg2png(dlsfnms[[i]], datUnits = datUnits, bQuantile = FALSE, 
                   bDiagNA = bDiagNA, diagColor = diagColor, 
                   zmean = 0.0, zrange = zrange, cutoff = cutoff)
    #return(ddat)
    ddat = invert_ddat_v4(ddat)
    idat = rdat * ddat
    lsdats[[i]] = idat
  }
  
  dmat = diffsummary(lsdats, bPNG = bPNG, pngfnm = pngfnm, pngtitle = pngtitle, 
              rcnames = rcnames, bColorBar = bColorBar, datUnits = datUnits, 
              pngWidth = pngWidth, pngHeight = pngHeight, 
              fromRes = fromRes, toRes = toRes, bGrid = bGrid, 
              ospace = ospace, bQuantile = bQuantile, bDiagNA = bDiagNA, 
              diagColor = diagColor, fWidth = fWidth, fHeight = fHeight, 
              zmean = zmean, zrange = zrange, cutoff = cutoff, bReturn = TRUE)
  
  if (avgDiffPNG != "")
  {
    dat = dat2minabs(dmat[[3, 1]], dmat[[3, 2]])
    
    # Subfunction does the work.
    dat = dat2png(dat, bPNG = bPNG, pngfnm = avgDiffPNG, 
                  pngWidth = fWidth, pngHeight = fHeight, bPlot = TRUE,
                  fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                  bColorBar = bColorBar, 
                  bDiagNA = bDiagNA, diagColor = diagColor, 
                  zmean = zmean, zrange = zrange, cutoff = cutoff, 
                  datUnits = datUnits, bQuantile = bQuantile, 
                  pngtitle = "Mean Difference")
  }
}


######### Sub-functions ########

# Make the list of file names. To use a subset, set includeFiles equal
# to an array with only the numbers corresponding to the wanted files.
# The array must be sorted, or it will leave out files.
mk_lsfnms = function(includeFiles = c(1, 2, 3, 4, 5), addPattern = "")
{
  lsfnms = list()
  if (any(includeFiles == 1))
  {
    lsfnms[[which(includeFiles == 1)]] = dir(pattern = paste("^M1\\.", ".*", addPattern, sep = ""))
  }
  if (any(includeFiles == 2))
  {
    lsfnms[[which(includeFiles == 2)]] = dir(pattern = paste("^M2toM1\\.", ".*", addPattern, sep = ""))
  }
  if (any(includeFiles == 3))
  {
    lsfnms[[which(includeFiles == 3)]] = dir(pattern = paste("^M2FBP\\.", ".*", addPattern, sep = ""))
  }
  if (any(includeFiles == 4))
  {
    lsfnms[[which(includeFiles == 4)]] = dir(pattern = paste("^M2\\.", ".*", addPattern, sep = ""))
  }
  if (any(includeFiles == 5))
  {
    lsfnms[[which(includeFiles == 5)]] = dir(pattern = paste("^M2H\\.", ".*", addPattern, sep = ""))
  }
  return(lsfnms)
}


dat2minabs = function(dat1, dat2)
{
  for (i in 1:length(dat1))
  {
    if (abs(dat2[i]) < abs(dat1[i]))
    {
      dat1[i] = dat2[i]
    }
  }
  return(dat1)
}


invert_ddat_v1 = function(dat)
{
  # Invert by quantilization.
  f   = ecdf(dat)
  dat = 1.0 - f(dat)
  
  return(dat)
}


invert_ddat_v2 = function(dat)
{
  # Invert by exponential.
  x   = mean(dat, na.rm = TRUE)
  dat = 2 ^ -(dat / x)
  
  return(dat)
}


invert_ddat_v3 = function(dat)
{
  # Invert by division.
  x   = mean(dat, na.rm = TRUE)
  dat = x / dat
  
  # Outliers.
  dat[dat > 100.0] = 100.0
  
  return(dat)
}


invert_ddat_v4 = function(dat)
{
  # Invert by division.
  x   = mean(dat, na.rm = TRUE)
  dat = x / dat
  
  # Outliers.
  dat[dat > 100.0] = 100.0
  
  # Normalize using sigma equation.
  dat = r2sigma(dat)
  
  return(dat)
}


invert_ddat_v5 = function(dat)
{
  # Invert by fitting to normal distribution.
  x   = mean(dat, na.rm = TRUE)
  xsd = sd(dat, na.rm = TRUE)
  dat = pnorm(dat, mean = x, sd = xsd, lower.tail = FALSE)
  
  return(dat)
}


r2sigma = function(r, rsigma = 1.0)
{
  s = 2.0 / (1.0 + exp(-r / rsigma)) - 1.0
  
  return(s)
}


############################### User Functions ###############################

# lsdats:    List of dat matrices. Each item has a corresponding row and col.
#            Function calculates all possible paired differences.
# rcnames:   Vector of row and column names. Should be same length as lsdats.
# pngWidth:  Number of pixels. Overrides fWidth.
# pngHeight: Number of pixels. Overrides fHeight.
# bPNG:      Write to png or draw to screen.
# fWidth:    Automatically calculates the width of the png.
# fHeight:   Automatically calculates the height of the png.
# bQuantile: Use normalized quantiles instead of the input data.
# bColorBar: Don't change this for now.
diffsummary = function(lsdats, bPNG = TRUE, pngfnm = "diffsummary.png", 
                       pngtitle = "", rcnames = NA, bColorBar = FALSE, 
                       fromRes = NA, toRes = NA, bGrid = FALSE, 
                       pngWidth = NA, pngHeight = NA, ospace = 0, 
                       bQuantile = FALSE, bDiagNA = TRUE, diagColor = 0.5, 
                       fWidth = 500, fHeight = 500, zmean = 0.0, zrange = 0.0, 
                       cutoff = 0.0, datUnits = "", 
                       bDQuantile = TRUE, bReturn = FALSE)
{
  ndats = length(lsdats)
  if (bPNG)
  {
    if (is.na(pngWidth))
    {
      pngWidth = fWidth * ndats
    }
    if (is.na(pngHeight))
    {
      pngHeight = fHeight * ndats
    }
    png(pngfnm, width = pngWidth, height = pngHeight)
  }
  fmat = matrix(data = 1:(ndats * ndats), byrow = TRUE, 
                nrow = ndats, ncol = ndats)
  layout(fmat)
  # Space for labeling the main figure.
  par(oma = c(ospace, ospace, ospace, 0))
  
  # Save the data to a list matrix if required. Creates matrix.
  if (bReturn)
  {
    dmat = as.list(numeric(ndats * ndats))
    dim(dmat) = c(ndats, ndats)
  }
  
  # Loops through all files and plots them.
  for (i in 1:ndats)
  {
    idat = lsdats[[i]]
    
    if (bDQuantile)
    {
      idat = quantilize(idat)
    }
    
    for (j in 1:ndats)
    {
      jdat = lsdats[[j]]
      
      if (bDQuantile)
      {
        jdat = quantilize(jdat)
      }
      
      ddat = diff2png(idat, jdat, datUnits = datUnits, bQuantile = bQuantile, 
                      bDiagNA = bDiagNA, diagColor = diagColor, 
                      zmean = zmean, zrange = zrange, cutoff = cutoff, 
                      fromRes = fromRes, toRes = toRes, bGrid = bGrid)
      # Save a copy of difference to the output matrix.
      if (bReturn)
      {
        dmat[[i, j]] = ddat
      }
    }
  }
  
  # Calculate positions of the labels.
  rcnpos = NA
  if (!is.na(rcnames[1]))
  {
    rcsize = 1 / ndats
    rcnpos = seq(0.0, 1.0, by = rcsize) - (rcsize / 2.0)
    rcnpos = rcnpos[-1]
  }
  
  # Figure labels.
  mtext(pngtitle, outer = TRUE, cex = 2.0, line = ospace %/% 2)
  mtext(rcnames, outer = TRUE, side = 1, cex = 2.0, at = rcnpos, 
        line = ospace %/% 2) 
  mtext(rev(rcnames), outer = TRUE, side = 2, cex = 2.0, at = rcnpos, 
        line = ospace %/% 2)
  
  # Close image.
  if (bPNG)
  {
    dev.off()
  }
  
  # Return matrix.
  if (bReturn)
  {
    return(dmat)
  }
}


# Tries to run diffsummary on a list of file names. The input should be a 
# list where each element is a vector of file names where the average 
# corresponds to one row and column.
fnms2diffsummary = function(lsfnms, bPNG = TRUE, pngfnm = "diffsummary.png", 
                            pngtitle = "", rcnames = NA, bColorBar = FALSE, 
                            fromRes = NA, toRes = NA, bGrid = FALSE, 
                            pngWidth = NA, pngHeight = NA, ospace = 0, 
                            bQuantile = TRUE, bDiagNA = TRUE, diagColor = 0.5, 
                            fWidth = 500, fHeight = 500, zmean = 0.0, 
                            zrange = 0.0, cutoff = 0.0, datUnits = "")
{
  lsdats = list()
  for (i in 1:length(lsfnms))
  {
    idat = avg2png(lsfnms[[i]], datUnits = datUnits, bQuantile = bQuantile, 
                   bDiagNA = bDiagNA, diagColor = diagColor, 
                   zmean = zmean, zrange = zrange, cutoff = cutoff)
    lsdats[[i]] = idat
  }
  
  diffsummary(lsdats, bPNG = bPNG, pngfnm = pngfnm, pngtitle = pngtitle, 
              rcnames = rcnames, bColorBar = bColorBar, datUnits = datUnits, 
              pngWidth = pngWidth, pngHeight = pngHeight, 
              fromRes = fromRes, toRes = toRes, bGrid = bGrid, 
              ospace = ospace, bQuantile = bQuantile, bDiagNA = bDiagNA, 
              diagColor = diagColor, fWidth = fWidth, fHeight = fHeight, 
              zmean = zmean, zrange = zrange, cutoff = cutoff)
}


# Calculates fnm1 - fnm2.
dfnm2png = function(fnm1, fnm2, bPNG = FALSE, pngfnm = "diffcorr.png", 
                    pngWidth = 1400, pngHeight = 1000, bPlot = TRUE, 
                    fromRes = NA, toRes = NA, bGrid = FALSE, 
                    bColorBar = FALSE, bHist = FALSE, bQHist = FALSE, 
                    bDiagNA = TRUE, diagColor = 0.5, bQuantile = TRUE, 
                    zmean = 0.0, zrange = 0.0, cutoff = 0.0, 
                    datUnits = "", pngtitle = "", bDQuantile = TRUE)
{
  # Load files.
  dat1 = as.matrix(read.table(fnm1))
  dat2 = as.matrix(read.table(fnm2))
  
  # Take difference of quantiles?
  if (bDQuantile)
  {
    dat1 = quantilize(dat1)
    dat2 = quantilize(dat2)
  }
  
  dat  = diff2png(dat1, dat2, bPNG = bPNG, pngfnm = pngfnm, 
                  pngWidth = pngWidth, pngHeight = pngHeight, bPlot = bPlot,
                  fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                  bColorBar = bColorBar, bHist = bHist, bQHist = bQHist, 
                  bDiagNA = bDiagNA, diagColor = diagColor, 
                  zmean = zmean, zrange = zrange, cutoff = cutoff, 
                  datUnits = datUnits, pngtitle = pngtitle, 
                  bQuantile = bQuantile)
  
  return(dat)
}


# Calculates dat1 - dat2.
diff2png = function(dat1, dat2, bPNG = FALSE, pngfnm = "diffcorr.png", 
                    pngWidth = 1400, pngHeight = 1000, bPlot = TRUE, 
                    fromRes = NA, toRes = NA, bGrid = FALSE, 
                    bColorBar = FALSE, bHist = FALSE, bQHist = FALSE, 
                    bDiagNA = TRUE, diagColor = 0.5, bQuantile = FALSE, 
                    zmean = 0.0, zrange = 0.0, cutoff = 0.0, 
                    datUnits = "", pngtitle = "")
{
  # If dat is intended to be the difference of quantiles, the function 
  # quantilize(dat) should be used on the inputs.
  dat  = dat1 - dat2
  
  # Subfunction does the work.
  dat = dat2png(dat, bPNG = bPNG, pngfnm = pngfnm, pngtitle = pngtitle, 
                pngWidth = pngWidth, pngHeight = pngHeight, bPlot = bPlot,
                fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                bColorBar = bColorBar, bHist = bHist, bQHist = bQHist, 
                bDiagNA = bDiagNA, diagColor = diagColor, 
                zmean = zmean, zrange = zrange, cutoff = cutoff, 
                datUnits = datUnits)
  
  return(dat)
}


# lsfnms:    List of vectors of filenames. Each list item will go on its own 
#            row. Each item in the vector will go on its own column.
# nrows:     Should be equal to the number of list items.
# ncols:     Should be equal to the number of elements in the largest vector.
# rnames:    Vector of row names. Should have length nrows.
# cnames:    Vector of column names. Should have length ncols.
# pngWidth:  Number of pixels. Overrides fWidth.
# pngHeight: Number of pixels. Overrides fHeight.
# bPNG:      Write to png or draw to screen.
# fWidth:    Automatically calculates the width of the png.
# fHeight:   Automatically calculates the height of the png.
# bAVG:      Adds an additional figure to the end of each row with the average 
#            of all the figures in the row. Automatically adds one to ncols.
# bQuantile: Use normalized quantiles instead of the input data.
# bColorBar: Don't change this for now.
datsummary = function(lsfnms, nrows, ncols, pngfnm = "datsummary.png", 
                      pngtitle = "", rnames = NA, cnames = NA, 
                      pngWidth = NA, pngHeight = NA, bPNG = TRUE, 
                      fromRes = NA, toRes = NA, bGrid = FALSE, 
                      fWidth = 500, fHeight = 500, bAVG = TRUE, 
                      bQuantile = TRUE, bDiagNA = TRUE, diagColor = NA, 
                      bColorBar = FALSE, datUnits = "", 
                      zmean = 0.0, zrange = 0.0, cutoff = 0.0, ospace = 0)
{
  if (bPNG)
  {
    if (is.na(pngWidth))
    {
      pngWidth = fWidth * (ncols + bAVG)
    }
    if (is.na(pngHeight))
    {
      pngHeight = fHeight * nrows
    }
    png(pngfnm, width = pngWidth, height = pngHeight)
  }
  fmat = matrix(data = 1:(nrows * (ncols + bAVG)), byrow = TRUE, 
                nrow = nrows, ncol = ncols + bAVG)
  layout(fmat)
  # Space for labeling the main figure.
  par(oma = c(ospace, ospace, ospace, 0))
  
  # Loops through all files and plots them.
  for (i in 1:nrows)
  {
    fnmlist = lsfnms[[i]]
    nFiles  = length(fnmlist)
    for (j in 1:ncols)
    {
      if (j > nFiles)
      {
        plot.new()
      }
      else
      {
        fnm  = fnmlist[j]
        idat = fnm2png(fnm, datUnits = datUnits, bQuantile = bQuantile, 
                       bDiagNA = bDiagNA, diagColor = diagColor,
                       fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                       zmean = zmean, zrange = zrange, cutoff = cutoff)
      }
    }
    # Average plot for the row.
    if (bAVG)
    {
      idat = avg2png(fnmlist, datUnits = datUnits, bQuantile = bQuantile, 
                     bDiagNA = bDiagNA, diagColor = diagColor,
                     fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                     zmean = zmean, zrange = zrange, cutoff = cutoff)
    }
  }
  
  # Calculate positions of the labels.
  cnpos = NA
  rnpos = NA
  if (!is.na(cnames[1]))
  {
    csize = 1 / (ncols + bAVG)
    cnpos = seq(0.0, 1.0, by = csize) - (csize / 2.0)
    cnpos = cnpos[-1]
    if (bAVG)
    {
      cnames[ncols + bAVG] = "Average"
    }
  }
  if (!is.na(rnames[1]))
  {
    rsize = 1 / nrows
    rnpos = seq(0.0, 1.0, by = rsize) - (rsize / 2.0)
    rnpos = rnpos[-1]
  }
  
  # Figure labels.
  mtext(pngtitle, outer = TRUE, cex = 2.0, line = ospace %/% 2)
  mtext(cnames, outer = TRUE, side = 1, cex = 2.0, at = cnpos, 
        line = ospace %/% 2) 
  mtext(rev(rnames), outer = TRUE, side = 2, cex = 2.0, at = rnpos, 
        line = ospace %/% 2)
  
  if (bPNG)
  {
    dev.off()
  }
}


# Returns outliers.
datoutlier = function(lsfnms, nrows, ncols, pngfnm = "datoutlier.png", 
                      pngtitle = "", rnames = NA, cnames = NA, 
                      pngWidth = NA, pngHeight = NA, bPNG = TRUE, 
                      fromRes = NA, toRes = NA, bGrid = FALSE, 
                      fWidth = 500, fHeight = 500, bAVG = TRUE, 
                      bQuantile = TRUE, bDiagNA = TRUE, diagColor = NA, 
                      bColorBar = FALSE, datUnits = "", 
                      zmean = 0.0, zrange = 0.0, cutoff = 0.0, ospace = 0, 
                      pval = 0.05)
{
  nvals  = 0
  mdiff  = c()
  lsdiff = list()
  # Loops through all files and does t-tests.
  for (i in 1:nrows)
  {
    fnmlist = lsfnms[[i]]
    nFiles  = length(fnmlist)
    lsdat   = list()
    for (j in 1:ncols)
    {
      if (j <= nFiles)
      {
        fnm   = fnmlist[j]
        jdat  = as.matrix(read.table(fnm))
        msize = nrow(jdat)
        if (bDiagNA)
        {
          diag(jdat) = NA
        }
        
        ######### Probably should use sub-function. ##########
        # Convert to quantiles.
        if (bQuantile)
        {
          f    = ecdf(jdat)
          jdat = qnorm(f(jdat))
          dim(jdat) = c(msize, msize)
          
          # Outliers.
          jdat[jdat >  5.0] =  5.0
          jdat[jdat < -5.0] = -5.0
        }
        ######################################################
        lsdat[[j]] = jdat
      }
    }
    
    idiff = matrix(nrow = nFiles, ncol = nFiles - 1)
    for (j in 1:nFiles)
    {
      jdat = lsdat[[j]]
      n    = 1
      for (k in 1:nFiles)
      {
        if (j == k)
        {
          next
        }
        kdat = lsdat[[k]]
        idiff[j, n] = mean(abs(kdat - jdat), na.rm = TRUE)
        n = n + 1
      }
    }
    lsdiff[[i]] = idiff
    
    nvals = nvals + nFiles
    mdiff = c(mdiff, as.vector(lsdiff[[i]]))
  }
  
  pvals = numeric(nvals)
  n = 1
  for (i in 1:nrows)
  {
    fnmlist = lsfnms[[i]]
    nFiles  = length(fnmlist)
    for (j in 1:nFiles)
    {
      pvals[n] = t.test(as.vector(lsdiff[[i]][j,]), mdiff, alternative = "greater")$p.value
      n = n + 1
    }
  }
  
  if (bPNG)
  {
    if (is.na(pngWidth))
    {
      pngWidth = fWidth * (ncols + bAVG)
    }
    if (is.na(pngHeight))
    {
      pngHeight = fHeight * nrows
    }
    png(pngfnm, width = pngWidth, height = pngHeight)
  
    fmat = matrix(data = 1:(nrows * (ncols + bAVG)), byrow = TRUE, 
                  nrow = nrows, ncol = ncols + bAVG)
    layout(fmat)
    
    # Space for labeling the main figure.
    par(oma = c(ospace, ospace, ospace, 0))
    
    # Loops through all files and plots them.
    n = 1
    for (i in 1:nrows)
    {
      fnmlist = lsfnms[[i]]
      nFiles  = length(fnmlist)
      ifnm    = c()
      for (j in 1:ncols)
      {
        if (j > nFiles)
        {
          plot.new()
        }
        else if (pvals[n] < pval)
        {
          plot.new()
          n = n + 1
        }
        else
        {
          fnm  = fnmlist[j]
          ifnm = c(ifnm, fnm)
          idat = fnm2png(fnm, datUnits = datUnits, bQuantile = bQuantile, 
                         bDiagNA = bDiagNA, diagColor = diagColor,
                         fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                         zmean = zmean, zrange = zrange, cutoff = cutoff)
          n = n + 1
        }
      }
    
      
      # Average plot for the row.
      if (bAVG)
      {
        idat = avg2png(ifnm, datUnits = datUnits, bQuantile = bQuantile, 
                       bDiagNA = bDiagNA, diagColor = diagColor,
                       fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                       zmean = zmean, zrange = zrange, cutoff = cutoff)
      }
    }
    
    # Calculate positions of the labels.
    cnpos = NA
    rnpos = NA
    if (!is.na(cnames[1]))
    {
      csize = 1 / (ncols + bAVG)
      cnpos = seq(0.0, 1.0, by = csize) - (csize / 2.0)
      cnpos = cnpos[-1]
      if (bAVG)
      {
        cnames[ncols + bAVG] = "Average"
      }
    }
    if (!is.na(rnames[1]))
    {
      rsize = 1 / nrows
      rnpos = seq(0.0, 1.0, by = rsize) - (rsize / 2.0)
      rnpos = rnpos[-1]
    }
    
    # Figure labels.
    mtext(pngtitle, outer = TRUE, cex = 2.0, line = ospace %/% 2)
    mtext(cnames, outer = TRUE, side = 1, cex = 2.0, at = cnpos, 
          line = ospace %/% 2) 
    mtext(rev(rnames), outer = TRUE, side = 2, cex = 2.0, at = rnpos, 
          line = ospace %/% 2)
    
    dev.off()
  }  
  return(pvals)
}


# fnmlist: Vector of file names as returned by the dir() function.
avg2png = function(fnmlist, bPNG = FALSE, pngfnm = "avgcor.png", 
                   pngWidth = 1400, pngHeight = 1000, bPlot = FALSE, 
                   fromRes = NA, toRes = NA, bGrid = FALSE, 
                   bColorBar = FALSE, bHist = FALSE, bQHist = FALSE, 
                   bDiagNA = TRUE, diagColor = NA, bQuantile = TRUE, 
                   zmean = 0.0, zrange = 0.0, cutoff = 0.0, 
                   datUnits = "", pngtitle = "")
{
  nFiles = length(fnmlist)
  d      = dim(as.matrix(read.table(fnmlist[1])))
  dat    = matrix(data = 0, nrow = d[1], ncol = d[2])
  
  for (i in 1:nFiles)
  {
    fnm = fnmlist[i]
    # This actually doesn't do anything except load the file unless bQuantile 
    # or bDiagNA are set to TRUE.
    idat = fnm2png(fnm, bPlot = FALSE, bDiagNA = bDiagNA, 
                   bQuantile = bQuantile)
    
    # Add to sum.
    dat  = dat + idat
  }
  
  # Mean of dat.
  dat = dat / nFiles
  
  # Subfunction does the work.
  dat = dat2png(dat, bPNG = bPNG, pngfnm = pngfnm, pngtitle = pngtitle, 
                pngWidth = pngWidth, pngHeight = pngHeight, bPlot = bPlot,
                fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
                bColorBar = bColorBar, bHist = bHist, bQHist = bQHist, 
                bDiagNA = bDiagNA, diagColor = diagColor, 
                zmean = zmean, zrange = zrange, cutoff = cutoff, 
                datUnits = datUnits, bQuantile = bQuantile)
  
  return(dat)
}


fnm2png = function(fnm, bPNG = FALSE, pngfnm = "correl.png", pngtitle = "", 
                   pngWidth = 1400, pngHeight = 1000, bPlot = TRUE, 
                   bColorBar = FALSE, bHist = FALSE, bQHist = FALSE, 
                   bDiagNA = TRUE, diagColor = NA, bQuantile = TRUE, 
                   zmean = 0.0, zrange = 0.0, cutoff = 0.0, 
                   datUnits = "", fromRes = NA, toRes = NA, bGrid = FALSE)
{
  # Load file.
  dat = as.matrix(read.table(fnm))
  
  # Subfunction does the work.
  dat = dat2png(dat, bPNG = bPNG, pngfnm = pngfnm, pngtitle = pngtitle, 
        pngWidth = pngWidth, pngHeight = pngHeight, bPlot = bPlot,
        fromRes = fromRes, toRes = toRes, bGrid = bGrid,  
        bColorBar = bColorBar, bHist = bHist, bQHist = bQHist, 
        bDiagNA = bDiagNA, diagColor = diagColor, bQuantile = bQuantile, 
        zmean = zmean, zrange = zrange, cutoff = cutoff, 
        datUnits = datUnits)
  
  return(dat)
}


############################### Subfunctions ################################# 

###############################
#
# Main subfunction
#
############ Inputs
#
# dat:       Either raw matrix from .dat file or processed matrix.
# bPNG:      Save a png?
# pngfnm:    Only used if bPNG is TRUE.
# pngtitle:  Only used if bPNG is TRUE.
# pngWidth:  Number of pixels. Only used if bPNG is TRUE.
# pngHeight: Number of pixels. Only used if bPNG is TRUE.
# bPlot:     Draw a heat map?
# bColorBar: Draw a color bar? Only if bPlot is TRUE.
# bHist:     Draw histogram of input data.
# bQHist:    Draw histogram of normalized quantile data.
# bDiagNA:   Sets diagonal to NA before processing. This means that the 
#            diagonal will have no effect on the histograms, the color scale, 
#            or the normalization to quantiles.
# diagColor: After processing, overrides the value of the diagonal to this 
#            value. At this point, the matrix has been converted to the color 
#            map, so that value should be between 0 and 1.
#
#            For a difference plot, a value of 0.5 makes sense.
#            For a correlation or MI plot, a value of 1.0 makes sense. 
#            For a distance plot, a value of 0.0 makes sense.
#
# bQuantile: Normalizes the input data into quantiles and maps to the cdf.
# zmean:     Set the mean of the output color. This value in the input data 
#            will be mapped to 0.5 in the color map. Should be in units of 
#            the input data if bQuantile is FALSE or in z-score units if 
#            bQuantile is TRUE.
# zrange:    Maximum range of the color plot, so zmean - zrange is mapped to 
#            0.0 and zmean + zrange is mapped to 1.0. Should be in units of 
#            the input data if bQuantile is FALSE or in z-score units if 
#            bQuantile is TRUE.
# cutoff:    Values below the cutoff are not plotted in the color map or in 
#            the color bar. They are set outside of the color range which 
#            shows up as white instead.
# datUnits:  Text. Name of the units of the data matrix. Used as xlab for 
#            histograms and the color bar.
#
############ Outputs
#
# Two png figures are intended. Two histograms in one figure can be obtained 
# with bPlot = FALSE as well as bQuantile = TRUE, bHist = TRUE, and bQHIST = 
# TRUE. Color map with optional color bar and histograms can be obtained with 
# bPlot = TRUE. The final processed data matrix is returned. Be aware that 
# bPlot = TRUE results in mapping the input data to a range [0, 1] color 
# scale, so the returned data matrix is not useful for further processing.
# Getting the normalized z-score of the data can be achieved by setting 
# bPlot = FALSE and bQuantile = TRUE.
#
############ Function
dat2png = function(dat, bPNG = FALSE, pngfnm = "correl.png", pngtitle = "", 
                   pngWidth = 1400, pngHeight = 1000, bPlot = TRUE, 
                   bColorBar = FALSE, bHist = FALSE, bQHist = FALSE, 
                   bDiagNA = TRUE, diagColor = NA, bQuantile = TRUE, 
                   zmean = 0.0, zrange = 0.0, cutoff = 0.0, 
                   datUnits = "", fromRes = NA, toRes = NA, bGrid = FALSE)
{
  # Axes labels.
  msize = nrow(dat)
  x     = c(13:(msize+12))
  y     = c(13:(msize+12))
  
  if (bPNG)
  {
    png(pngfnm, width = pngWidth, height = pngHeight)
    if (bPlot)
    {
      if (bHist)
      {
        layout(matrix(c(0, 0, 0, 3, 3, 4, 1, 2, 0), nrow = 3, ncol = 3), 
               widths = c(1, 10, 5), heights = c(7, 7, 2))
      }
      else
      {
        layout(matrix(c(0, 0, 1, 2), nrow = 2, ncol = 2), 
               widths = c(1, 10), heights = c(14, 2))
      }
    }
    else
    {
      layout(matrix(c(1, 2), nrow = 1, ncol = 2))
    }
  }
  
  # Set diagonal to NA.
  if (bDiagNA)
  {
    diag(dat) = NA
  }
  
  if (bHist)
  {
    # Record the data range.
    if ((zrange > 0.0) & (!bQuantile))
    {
      drange = c(zmean - zrange, zmean + zrange)
    }
    else
    {
      drange = c(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE))
    }
    drawHist(dat, drange, "Distribution", datUnits)
  }
  else
  {
    if (bPNG)
    {
      plot.new()
    }
  }
  
  # Convert to quantiles.
  if (bQuantile)
  {
    dat = quantilize(dat)
    
    # Display histogram after quantilization.
    if (bQHist)
    {
      # Record the data range.
      if (zrange > 0.0)
      {
        drange = c(zmean - zrange, zmean + zrange)
      }
      else
      {
        drange = c(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE))
      }
      
      # Draw histogram.
      drawHist(dat, drange, "Quantiles", paste("Z Score of ", datUnits))
    }
    else
    {
      if (bPNG)
      {
        plot.new()
      }
    }
  }
  else
  {
    if (bPNG)
    {
      plot.new()
    }
  }

  # Save a copy of dat to return.
  ret = dat
  
  # The rest of the function is unnecessary unless a plot is being generated.
  if (bPlot)
  {
    # Convert to the range [0, 1] for color mapping.
    if (zrange > 0.0)
    {
      # Record the data range.
      drange = c(zmean - zrange, zmean + zrange)
      dat    = ((dat - zmean) / (2 * zrange)) + 0.5
      
      # Set values outside of the range to threshold.
      dat[dat > 1.0] = 1.0
      dat[dat < 0.0] = 0.0
    }
    else
    {
      # Record the data range.
      drange = c(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE))
      dat    = (dat - drange[1]) / (drange[2] - drange[1])
    }
    
    # Override color on the diagonal.
    if (!is.na(diagColor))
    {
      diag(dat) = diagColor
    }
    
    if (cutoff > 0.0)
    {
      # Calculate color cutoffs.
      if (zrange > 0.0)
      {
        clim    = c()
        clim[1] = 0.5 - cutoff / (2 * zrange)
        clim[2] = 0.5 + cutoff / (2 * zrange)
      }
      else
      {
        clim    = c()
        clim[1] = 0.5 - cutoff / (drange[2] - drange[1])
        clim[2] = 0.5 + cutoff / (drange[2] - drange[1])
      }
      
      # Values that do not meet the cutoff are set to white.
      dat[(dat > clim[1]) & (dat < clim[2])] = -1.0
      cat("Number that do not meet cutoff: ", 
          sum(dat == -1.0, na.rm = TRUE), "\n")
    }
    else
    {
      clim = NA
    }
    
    # Draw heat map and color bar.
    image(x, y, dat, zlim = c(0, 1), main = pngtitle, useRaster = TRUE, 
          xlab = "Residue Number", ylab="Residue Number", cex.main = 2.0, 
          cex.axis = 1.5, cex.lab = 1.5, col = rev(rainbow(1000, end = 0.66)))
    if (bGrid)
    {
      tick_length = 1
      tick_width  = 1
    }
    else
    {
      tick_length = -0.02
      tick_width  = 4
    }
    if (!is.na(fromRes[1]))
    {
      axis(4, at = fromRes, labels = NA, lwd = 0, lwd.ticks = tick_width, las = 2,
           tck = tick_length, col.ticks = "blue")
    }
    if (!is.na(toRes[1]))
    {
      axis(3, at =   toRes, labels = NA, lwd = 0, lwd.ticks = tick_width, las = 2,
           tck = tick_length, col.ticks = "red")
    }
    if (bColorBar)
    {
      drawColorBar(drange, clim, datUnits)
    }
  }
  
  if (bPNG)
  {
    dev.off()
  }
  return(ret)
}


quantilize = function(dat)
{
  # Store dimensions. 
  dims = dim(dat)
  
  # Do not quantilize the values on the diagonal which are all the same. 
  diag(dat) = NA
  
  # Do nothing to an empty matrix.
  if (sum(abs(dat), na.rm = TRUE) > 0.0)
  {
    # Convert to quantiles.
    f   = ecdf(dat)
    dat = qnorm(f(dat))
    dim(dat) = dims
    
    # Outliers.
    dat[dat >  5.0] =  5.0
    dat[dat < -5.0] = -5.0
  }
  
  return(dat)
}


drawHist = function(dat, drange = NA, histTitle = "", datUnits = "")
{
  if (is.na(drange[1]))
  {
    drange = c(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE))
  }
  hist(dat, main = histTitle, xlim = c(drange[1], drange[2]), yaxt = 'n', 
       xlab = datUnits, ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
       breaks = "fd")
}


drawColorBar = function(drange, clim = NA, datUnits = "")
{
  xdat = seq(drange[1], drange[2], length.out=1000)
  ydat = seq(0, 1, length.out = 2)
  zdat = matrix(nrow = 1000, ncol = 2)
  zdat[,1] = seq(0, 1, length.out=1000)
  zdat[,2] = seq(0, 1, length.out=1000)
  if (!is.na(clim[1]))
  {
    zdat[(zdat > clim[1]) & (zdat < clim[2])] = -1.0
  }
  image(xdat, ydat, zdat, zlim = c(0,1), xlab = datUnits, ylab = "", 
        yaxt = 'n', cex.lab = 1.5, cex.axis = 1.5, useRaster = TRUE, 
        col = rev(rainbow(1000, end = 0.66)))
}


correl2dat = function(cdat, fnm)
{
  isize = dim(cdat)[1]
  jsize = dim(cdat)[2]
  sink(file = fnm)
  for (i in seq(1,isize))
  {
    for (j in seq(1,jsize))
    {
      cat(sprintf("%10.5f ",cdat[i,j]))
    }
    cat(sprintf("\n"))
  }
  sink()
}
