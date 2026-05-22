########################################################################
# Parse eye-tracking data in fixations and saccades
#
# Fixations and Saccades are detected using the SMI low speed event detection 
# algorithm (p. 219f, SMI manual)
#
# timeline = time in ms
# x and y = coordinates in pixels
# blk = blink vector (1=blink)
# sampl = Saccade threshold (for smaller saccades, fixations are combined)
# plotdata = Shall data be plotted
# trial = Trial label for the current plot
eventdetection_lowspeed <- function(timeline, x, y, blk, sampl, plotdata=FALSE, trial="") {
  hz <- 1000/mean(diff(timeline))
  
  mindurms <- 80     # Minimum duration of fixation in ms
  mindur   <- round(mindurms/1000*hz)   # in Samples
  maxdisp  <- 50     # Maximim dispersion in pixels
  
  # Set NAs to 0
  x[is.na(x)] <- 0
  y[is.na(y)] <- 0
  
  # Plot raw data?
  if (plotdata) {
    par(mfrow=c(2,1),mar=c(2,2,1,1))
    plot(timeline,x,type="n",ylim=c(min(c(x,y)),max(c(x,y))))
    # Plot interpolated regions
    abline(v=timeline[blk==1],col="grey80")
    lines(timeline,x,col="red2"); lines(timeline,y,col="blue")
    text(max(timeline),max(c(x,y)),trial,adj=c(1,1))
  }
  
  # Moving window for fixation detection
  trialfix <- numeric()    
  wst <- 1; wen <- (wst+mindur-1)
  while (wen<(length(x)-1)) {
    gazewin <- wst:wen
    d <- (max(x[gazewin])-min(x[gazewin])) + (max(y[gazewin])-min(y[gazewin]))
    if ((d<maxdisp) & (sum(blk[gazewin])==0)) {
      while ((d<maxdisp) & (sum(blk[gazewin])==0) & (wen<(length(x)-1))) {
        wen <- wen+1
        gazewin <- wst:wen
        d <- (max(x[gazewin])-min(x[gazewin])) + (max(y[gazewin])-min(y[gazewin]))
      } 
      trialfix <- rbind(trialfix,c(wst,wen-1,mean(x[wst:(wen-1)]),mean(y[wst:(wen-1)])))
      wst <- wen; wen <- (wst+mindur-1)
    } else {
      wst <- wst+1; wen <- (wst+mindur-1)
    }
  }    
  
  # Continue only if fixations were found
  if (length(trialfix)>0) {
    # Calculate fixation duration (in samples)
    trialfix <- cbind(trialfix,trialfix[,2]-trialfix[,1])
    
    trialfix <- cbind(trialfix,0,0)  # Add columns for saccade length and blinks
    if (nrow(trialfix) > 1) {
      for (i in 1:(nrow(trialfix)-1)) {
        sacwin <- trialfix[i,2]:trialfix[i+1,1]
        # Saccade length in pixels
        sacampl <- sqrt((trialfix[i+1,3]-trialfix[i,3])^2+(trialfix[i+1,4]-trialfix[i,4])^2)
        trialfix[i+1,6] <- sacampl
        if (sum(blk[sacwin])!=0) {
          trialfix[i+1,7] <- 1   # Blink vor gefundener Fixation
        }
      }
    }
    
    # Generate final list of fixations and saccades
    finalfix <- numeric()
    finalsac <- numeric()
    
    # Only one fixation
    if (nrow(trialfix) == 1) {
      finalfix <- rbind(finalfix,c(trialfix[1,1],trialfix[1,2],
                                   trialfix[1,3],trialfix[1,4],
                                   trialfix[1,7]))
    }
    
    # Iterativ die Liste der Fixationen durcharbeiten und Fixationen zusammenlegen, die weniger 
    # als sampl auseinanderliegen (nur wenn kein Blink dazwischen liegt):
    # < sampl : Continuous fixation
    # > sampl : Real saccade
    if (nrow(trialfix) > 1) {
      reprocess <- TRUE
      if (sum(trialfix[trialfix[,7]==0,6]<sampl)==0) { reprocess <- FALSE }
      
      st <- 1
      while ((st<=(nrow(trialfix)-1)) & reprocess) {
        # Blink gefunden -> mit naechster Fixation weiterarbeiten
        if (trialfix[st+1,7]==1) {
          st <- st + 1
        } else {
          # Nur kleine Fixationsaenderung -> Fixationen zusammenlegen
          if (trialfix[st+1,6]<sampl) {
            newfix <- c(trialfix[st,1],trialfix[st+1,2],
                        sum(trialfix[st:(st+1),3]*trialfix[st:(st+1),5])/sum(trialfix[st:(st+1),5]),
                        sum(trialfix[st:(st+1),4]*trialfix[st:(st+1),5])/sum(trialfix[st:(st+1),5]))
            newfix <- c(newfix,newfix[2]-newfix[1],0,0)
            
            trialfix[st+1,] <- newfix
            trialfix <- matrix(trialfix[(1:nrow(trialfix))!=st,],ncol=7)
            
            # Winkel neu berechnen
            if (nrow(trialfix)>1) {
              for (i in 1:(nrow(trialfix)-1)) {
                # Saccade length in degrees
                sacw  <- sqrt((trialfix[i+1,3]-trialfix[i,3])^2+(trialfix[i+1,4]-trialfix[i,4])^2)
                trialfix[i+1,6] <- sacw
                trialfix[i+1,7] <- ifelse(sum(blk[trialfix[i,2]:trialfix[i+1,1]])==0,0,1)
              }
            }
          } else {  # grosse Fixationsaenderung -> mit naechster Fixation weiterarbeiten
            st <- st + 1
          }
        }
        
        if (sum(trialfix[trialfix[,7]==0,6]<sampl)==0) { reprocess <- FALSE }
      }
      finalfix <- matrix(trialfix[,c(1,2,3,4,7)],ncol=5)
    }
    
    # Sakkadenliste erstellen
    if (nrow(finalfix)>1) {
      for (i in 1:(nrow(finalfix)-1)) {
        finalsac <- rbind(finalsac,c(finalfix[i,2],finalfix[i+1,1],
                                     finalfix[i,3],finalfix[i,4],
                                     finalfix[i+1,3],finalfix[i+1,4],
                                     finalfix[i+1,5]))
      }
    }
    
    # Fixationen und Sakkaden merken und zurueckgeben
    # vorher: Zeitstempel umrechnen Samples -> Zeit (ms)
    out <- list()
    if (length(finalfix)>0) { 
      finalfix[,1] <- timeline[finalfix[,1]]
      finalfix[,2] <- timeline[finalfix[,2]]				
      # Round position data to 2 decimals
      finalfix[,3:4] <- round(finalfix[,3:4],2)
      out$fixations <- finalfix[,1:4]
    }
    if (length(finalsac)>0) { 
      finalsac[,1] <- timeline[finalsac[,1]]
      finalsac[,2] <- timeline[finalsac[,2]]
      # Round position data to 2 decimals
      finalsac[,3:6] <- round(finalsac[,3:6],2)
      out$saccades  <- finalsac
    }
    
    # Plot detected fixations and saccades
    if (plotdata) {
      plot(timeline,x,type="l",col="red2",ylim=c(min(c(x,y)),max(c(x,y))))
      lines(timeline,y,col="blue")
      segments(finalfix[,1],finalfix[,3],finalfix[,2],finalfix[,3],col="black")
      segments(finalfix[,1],finalfix[,4],finalfix[,2],finalfix[,4],col="black")
      if (length(finalsac[finalsac[,7]==0])>0) { abline(v=finalsac[finalsac[,7]==0,1]) }
    }
  } else {
    out <- list()
    out$fixations <- NA
    out$saccades <- NA
  }
  
  return(out)
}
