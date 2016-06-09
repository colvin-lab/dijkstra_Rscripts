# fnm:     Name of dijkstra path data file.
# fromRes: Vector of residue numbers to search from.
# toRes:   Vector of residue numbers to search to.
# nres:    Number of residues in file.
dijkstra_search = function(fnm, fromRes, toRes, nres = 518, 
                           fnmPDF = "", titlePDF = "", 
                           leg = NA, legSpace = 0.4, type = "l",  
                           bSHOW = TRUE, yscale = "Frequency")
{
  # Options.
  bPDF = (fnmPDF != "")
  
  nodes = numeric(nres)
  for (i in fromRes)
  {
    print(i)
    
    # Read nres lines from the file and store them to dij_data.
    dij_data  = read_dij(fnm, i, nres, hlines = 1)
    print("finished reading")
    
    # Sum the node frequency.
    nodes = nodes + sum_dij_nodes(dij_data, nres, toRes, yscale)
    print("finished summing")
    print(nodes)
  }


  if (bPDF)
  {
    pdf(file = fnmPDF)
    par(mar = c(5, 4, 3.5, 1), oma = c(0, 1, 2, 0))
    
    if (!is.na(leg[1]))
    {
      par(xpd = TRUE)
      layout(matrix(c(1, 0), nrow = 1, ncol = 2), 
             widths = c(1, legSpace))
      plot.new()
      legend("topright", title = "Tick Labels", legend = leg, 
             fill = c("blue", "red"), inset = c(0.5, 0.0))
      par(xpd = FALSE)
    }
    doPlot(nodes,  fromRes, toRes, type, leg, yscale)
    mtext(titlePDF, cex = 1.2, outer = TRUE)
    dev.off()
  }

  if (bSHOW)
  {
    par(mar = c(5, 4, 3.5, 1), oma = c(0, 1, 2, 0))
  
    if (!is.na(leg[1]))
    {
      par(xpd = TRUE)
      layout(matrix(c(1, 0), nrow = 1, ncol = 2), 
             widths = c(1, legSpace))
      plot.new()
      legend("topright", title = "Tick Labels", legend = leg, 
             fill = c("blue", "red"), inset = c(0.5, 0.0))
      par(xpd = FALSE)
    }
    doPlot(nodes,  fromRes, toRes, type, leg, yscale)
    mtext(titlePDF, cex = 1.2, outer = TRUE)
  }

  return(nodes)
}


doPlot = function(x, fromRes, toRes, type, leg, ylab)
{
  plot(x, xlab = "Residue", ylab = ylab, type = type)
  axis(3, at = fromRes, labels = NA, lwd = 0, lwd.ticks = 2, las = 2, 
       col.ticks = "blue")
  #abline(v = fromRes, col = "blue", lwd = 1)
  axis(3, at =   toRes, labels = NA, lwd = 0, lwd.ticks = 2, las = 2, 
       col.ticks = "red")
  #abline(v =   toRes, col = "red",  lwd = 1)
}



outside_paths = function(infile, outfile, nChain = 4, 
                         nres = c(518, 518, 518, 518))
{
  sink(file = outfile)
  n = sum(nres)
  for (i in 1:nChain)
  {
    iLine = sum(nres[0:(i - 1)]) * n + sum(nres[0:(i - 1)])
    bRes  = sum(nres[0:(i - 1)])
    eRes  = sum(nres[0:i])
    for (j in 1:nres[i])
    {
      # jLine is one less than the line that the k loop starts on.
      jLine = iLine + (j - 1) * n
      for (k in 1:nres[i])
      {
        # kLine is one less than the line being read.
        kLine = jLine + k
        kText = scan(infile, what = "character", sep = ":", quiet = TRUE, 
                     strip.white = TRUE, skip = kLine, nlines = 1)
        kPath = as.numeric(strsplit(kText[6], " +")[[1]])
        pLeg  = length(kPath)
        if (pLeg > 3)
        {
          for (m in 3:(pLeg - 1))
          {
            mRes = kPath[m]
            if ((mRes < bRes) | (mRes >= eRes))
            {
              sRes = as.numeric(kText[1])
              tRes = as.numeric(kText[2])
              cat(sprintf("%6i %6i %6i \n", sRes, tRes, mRes))
            }
          }
        }
      }
    }
  }
  sink()
}


sum_dij_nodes = function(dij_data, nres, jres, yscale)
{
  A = test_yscale(dij_data, yscale, nres)
  nodes = numeric(nres)
  
  # Residue numbers use 0-based index, so add 1.
  for (j in (jres + 1))
  {
    N = length(dij_data$PATH_ATOMS[[j]])
    if (N < 3)
    {
      next
    }
    # Residue numbers use 0-based index, so add 1.
    P = dij_data$PATH_ATOMS[[j]] + 1
    # Sum for all steps on the path except the beginning and end.
    Pm        = P[2:(N - 1)]
    nodes[Pm] = nodes[Pm] + A[Pm]
    #for (k in 2:(N - 1))
    #{
    #  m = P[k]
    #  nodes[m] = nodes[m] + A[m]
    #}
  }
  return(nodes)
}


test_yscale = function(dij_data, yscale, nres)
{
  if (yscale == "Frequency")
  {
    A = rep(1, nres)
  }
  else if (yscale == "R")
  {
    A = dij_data$MEAN_R
  }
  else if (yscale == "S")
  {
    A = dij_data$SIGMA_R
  }
  else if (yscale == "Dijkstra")
  {
    A = dij_data$DIJ_DIST
  }
  else if (yscale == "R Dijkstra")
  {
    A = dij_data$MEAN_R * dij_data$DIJ_DIST
  }
  else if (yscale == "S Dijkstra")
  {
    A = dij_data$MEAN_S * dij_data$DIJ_DIST
  }
  else
  {
    stop("Argument for yscale is invalid.")
  }
  
  return(A)
}


read_dij = function(fnm, ires, nres = 518, hlines = 1)
{
  dij_data = list()
  ncol     = 12
  skip     = ires * nres + hlines
  
  # Matches the header line of the dij file.
  dij_data$START_ATOM  = numeric(nres)
  dij_data$TARGET_ATOM = numeric(nres)
  dij_data$MI_DIST     = numeric(nres)
  dij_data$DIJ_DIST    = numeric(nres)
  dij_data$NORM_DIST   = numeric(nres)
  dij_data$MEAN_R      = numeric(nres)
  dij_data$SIGMA_R     = numeric(nres)
  dij_data$PATH_ATOMS  = vector(mode = "list", length = nres)
  dij_data$PATH_MI     = vector(mode = "list", length = nres)
  dij_data$PATH_DIJ    = vector(mode = "list", length = nres)
  dij_data$PATH_R      = vector(mode = "list", length = nres)
  dij_data$TOTAL_R     = numeric(nres)
  
  dij_scan = scan(fnm, what = "character", sep = ":", strip.white = TRUE, 
                  quiet = TRUE, skip = skip, nlines = nres)
  
  for (i in 1:nres)
  {
    if (i == ires)
    {
      next
    }
    
    # Converts numerical fields.
    dij_data$START_ATOM[i]  = as.numeric(dij_scan[ncol * (i - 1) +  1])
    dij_data$TARGET_ATOM[i] = as.numeric(dij_scan[ncol * (i - 1) +  2])
    dij_data$MI_DIST[i]     = as.numeric(dij_scan[ncol * (i - 1) +  3])
    dij_data$DIJ_DIST[i]    = as.numeric(dij_scan[ncol * (i - 1) +  4])
    dij_data$NORM_DIST[i]   = as.numeric(dij_scan[ncol * (i - 1) +  5])
    dij_data$MEAN_R[i]      = as.numeric(dij_scan[ncol * (i - 1) +  6])
    dij_data$SIGMA_R[i]     = as.numeric(dij_scan[ncol * (i - 1) +  7])
    dij_data$TOTAL_R[i]     = as.numeric(dij_scan[ncol * (i - 1) + 12])
    
    
    # Converts vector fields.
    path_atoms = unlist(strsplit(dij_scan[ncol * (i - 1) +  8], ' +'))
    path_mi    = unlist(strsplit(dij_scan[ncol * (i - 1) +  9], ' +'))
    path_dij   = unlist(strsplit(dij_scan[ncol * (i - 1) + 10], ' +'))
    path_r     = unlist(strsplit(dij_scan[ncol * (i - 1) + 11], ' +'))
    
    n   = length(path_atoms)
    dij_data$PATH_ATOMS[[i]]    = numeric(n)
    dij_data$PATH_ATOMS[[i]][1] = as.numeric(path_atoms[1])
    dij_data$PATH_MI[[i]]       = numeric(n - 1)
    dij_data$PATH_DIJ[[i]]      = numeric(n - 1)
    dij_data$PATH_R[[i]]        = numeric(n - 1)
    for (j in 1:(n - 1))
    {
      dij_data$PATH_ATOMS[[i]][j + 1] = as.numeric(path_atoms[j + 1])
      dij_data$PATH_MI[[i]][j]        = as.numeric(path_mi[j])
      dij_data$PATH_DIJ[[i]][j]       = as.numeric(path_dij[j])
      dij_data$PATH_R[[i]][j]         = as.numeric(path_r[j])
    }
  }
  return(dij_data)
}