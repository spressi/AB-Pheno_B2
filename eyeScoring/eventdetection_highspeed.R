########################################################################
# Parse eye-tracking data in fixations and saccades
#
# High-speed algorithm similar to SR Research
# Saccades are detected with velocity and acceleration thresholds of 
# 30°/s and 8000°/s^2
#
# timeline = time in ms
# x and y = coordinates in pixels
# blk = blink vector (1=blink)
# pixperdeg = How many pixels amount to 1° of visual angle?
# plotdata = Shall data be plotted
# trial = Trial label for the current plot
eventdetection_highspeed <- function(timeline, x, y, blk, pixperdeg, plotdata=FALSE, trial="") {
  hz <- 1000/mean(diff(timeline))
  
  vth <- 30     # Velocity threshold
  ath <- 8000   # Acceleration threshold
  refrac <- 50  # Refractory period after a saccade in ms (to combine successive saccades)

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
  
  # Calculate velocity and acceleration in deg visual angle
  vel <- sqrt(c(0,(diff(x)/pixperdeg)^2+(diff(y)/pixperdeg)^2))/c(1,diff(timeline)/1000)
  acc <- c(0,diff(vel))/c(1,diff(timeline)/1000)
  
  # Define timepoints of potential saccades
  sac <- as.numeric(vel>vth | acc>ath)

  # Fill blink period between saccades
  sac[blk==1] <- 1
  
  # Analyze saccades only when they are available
  if (sum(sac)>0) {
    # Combine successive saccades and index saccade number
    sactp <- which(sac==1)
    saci  <- rep(0,length(sac))
    aktsaci <- 1
    saci[sactp[1]] <- aktsaci
    if (sum(sac)>1) {
      for (sacnr in 2:length(sactp)) {
        if ((timeline[sactp[sacnr]]-timeline[sactp[sacnr-1]])<refrac) {
          sac[sactp[sacnr-1]:sactp[sacnr]] <- 1
          saci[sactp[sacnr-1]:sactp[sacnr]] <- aktsaci
        } else {
          aktsaci <- aktsaci + 1
          saci[sactp[sacnr]] <- aktsaci
        }
      }
    }
    
    # Generate lists of saccades and fixations
    finalsac <- numeric()
    finalfix <- numeric()
    for (sacnr in 1:max(saci)) {
      sacst <- head(timeline[saci==sacnr],1)
      sacen <- tail(timeline[saci==sacnr],1)
      
      stpos <- which(timeline==sacst)[1]
      enpos <- which(timeline==sacen)[1]
      
      # Ignore first/last saccade when starting/ending point could not be determined
      if ((stpos>1) & (enpos<length(timeline))) {
        # Define saccade
        stx <- ifelse(stpos>1,x[stpos-1],x[stpos])
        sty <- ifelse(stpos>1,y[stpos-1],y[stpos])
        
        enx <- ifelse(enpos<length(timeline),x[enpos+1],x[enpos])
        eny <- ifelse(enpos<length(timeline),y[enpos+1],y[enpos])
    
        sacblk <- as.numeric(sum(blk[saci==sacnr])!=0)
        
        finalsac <- rbind(finalsac,c(sacst,sacen,stx,sty,enx,eny,sacblk))
      }
      
      # Define fixation
      # Is fixation available?
      # First fixaion (when available)
      if ((sacnr==1) & (stpos>1)) {
        fixst <- head(timeline,1)
        fixen <- timeline[stpos-1]
        
        fixx <- mean(x[(timeline>=fixst) & (timeline<=fixen)])
        fixy <- mean(y[(timeline>=fixst) & (timeline<=fixen)])
        
        finalfix <- rbind(finalfix,c(fixst,fixen,fixx,fixy))
      } 
      # Last fixation
      if ((sacnr==max(saci)) & (enpos<length(timeline))) {
        fixst <- timeline[enpos+1]
        fixen <- tail(timeline,1)

        fixx <- mean(x[(timeline>=fixst) & (timeline<=fixen)])
        fixy <- mean(y[(timeline>=fixst) & (timeline<=fixen)])
        
        finalfix <- rbind(finalfix,c(fixst,fixen,fixx,fixy))
      } 
      if (sacnr<max(saci)) {
        nextsacst <- head(timeline[saci==sacnr+1],1)
        nextstpos <- which(timeline==nextsacst)[1]
        
        fixst <- timeline[enpos+1]
        fixen <- timeline[nextstpos-1]
        
        fixx <- mean(x[(timeline>=fixst) & (timeline<=fixen)])
        fixy <- mean(y[(timeline>=fixst) & (timeline<=fixen)])
        
        finalfix <- rbind(finalfix,c(fixst,fixen,fixx,fixy))
      }
    }
  } else {   # No saccade available - is there a fixation?
    if (sum(blk)<length(blk)) {
      fixst <- head(timeline,1)
      fixen <- tail(timeline,1)
      
      fixx <- mean(x[(timeline>=fixst) & (timeline<=fixen)])
      fixy <- mean(y[(timeline>=fixst) & (timeline<=fixen)])
      
      finalsac <- numeric()
      finalfix <- numeric()
      finalfix <- rbind(finalfix,c(fixst,fixen,fixx,fixy))
    }
  }
  
  # Fixationen und Sakkaden merken und zurueckgeben
  # vorher: Zeitstempel umrechnen Samples -> Zeit (ms)
  out <- list()
  if (length(finalfix)>0) { 
    # Round position data to 2 decimals
    finalfix[,3:4] <- round(finalfix[,3:4],2)
    out$fixations <- finalfix
  } else {
    out$fixations <- NA
  }
  if (length(finalsac)>0) { 
    # Round position data to 2 decimals
    finalsac[,3:6] <- round(finalsac[,3:6],2)
    out$saccades <- finalsac
  } else {
    out$saccades <- NA
  }

  # Plot detected fixations and saccades
  if ((length(finalfix)>0) & (plotdata)) {
    plot(timeline,x,type="l",col="red2",ylim=c(min(c(x,y)),max(c(x,y))))
    lines(timeline,y,col="blue")
    segments(finalfix[,1],finalfix[,3],finalfix[,2],finalfix[,3],col="black")
    segments(finalfix[,1],finalfix[,4],finalfix[,2],finalfix[,4],col="black")
    if (length(finalsac)>0) {
      if (length(finalsac[finalsac[,7]==0])>0) { abline(v=finalsac[finalsac[,7]==0,1]) }
    }
  }
  
  return(out)
}
