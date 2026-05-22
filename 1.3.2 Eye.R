library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")
#source("1.3.1 Score Eye.R")

library("ggforce") #for drawing circles
library("png") #for reading png
library("DescTools")

{ # Global Parameters -------------------------------------------------------
  exclusions.eye = c() %>% c(exclusions) #a priori exclusions, e.g. calibration not successful
  
  #baseline validation
  baseline = c(-300, 0) #Baseline in ms relative to stimulus onset; min(baseline) = start; max(baseline) = end
  useAllBaselines = list() #manually allow all baselines for some participants?
  saveBaselinePlots = F
  driftPlots = T #c("vp30", "vp33")
  maxDeviation_rel = 3 #max abs value of z-score
  outlierLimit.eye = .5 #maximum percentage of invalid baselines per subject
  maxSpread = 150 #maximum spread of valid baselines (diameter / edge length)
  usePointDistance = T #use point (vector) distance instead of coordinates independently?
  
  #ROI analysis
  validFixTime.trial = .5 #percentage of valid fixation time within trial in order to analyze trial
  validFixTime.subj = .5 #percentage of trials with sufficient valid fixation time in order to analyze subject
  diagnosticDwell = .5 #percentage of trials per subject that need at least one fixation towards the diagnostic ROI
  showTrialPlots = F
  unifyRois = T
  
  screen.width <- 1920 ; screen.height <- 1080 #screen resolution (TODO: model & size of screen)
  pixsize <- 520 / screen.width  # Pixel size in mm (for scan path length)
  distance = 650 #screen center to subject eyes in mm
  
  indifference.zone = 100 #pixels away from center to count as looking to face
  
  z.max = 2 #Winsorize dependent variables to a z value of 2
  q.max = pnorm(z.max * c(-1, 1)) #needed for DescTools::Winsorize function
  q.max = c(0, 1) #switch off Winsorizing (comment out to switch on)
}

{ # Functions ---------------------------------------------------------------
  outlier_remove <- function(x, z=3) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 3) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)] #can contain NAs
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first minimum if several are present)
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)] #can contain NAs
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first maximum if several are present)
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - z*sd(subx))) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + z*sd(subx))) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  #removes values until (absolute) deviation is satisfied
  outlier_remove_spread = function(x, deviation, insample=NULL) {
    if (is.null(insample)) insample = T %>% rep(length(x))
    insample[is.na(x)] = F
    
    spread = x[insample] %>% range() %>% diff()
    
    while (spread > deviation) {
      m = x[insample] %>% mean()
      x[!insample] = m #set outliers to m to ignore them in next step
      out.next = {x - m} %>% abs() %>% which.max()
      
      insample[out.next] = F #crucial: use x, not x[insample] to match vector length of insample
      spread = x[insample] %>% range() %>% diff()
    }
    return(insample)
  }
    
  
  #TODO after outlier removal is finished, look once (?) more if "outliers" can be INCLUDED (can happen for biased distributions)
  outlier_remove_abs <- function(x, deviation) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 2) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)]
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1]
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)]
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1]
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - deviation)) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + deviation)) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  pointDistance = function(x1, y1, x2, y2) {
    return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
  }
  
  outlier_remove_point = function(x, y, z=3, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 3) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= mean(distances[insample]) + z*sd(distances[insample])) {
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  outlier_remove_point_spread = function(x, y, deviation, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point_spread: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 2) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= deviation/2) { #deviation = diameter but with distances only radius can be checked => deviation/2
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  pointInBounds = function(point.x, point.y, x.bounds, y.bounds) {
    if (is.na(point.x) || is.na(point.y) || length(point.x)==0 || length(point.y)==0) return(FALSE)
    return(point.x >= min(x.bounds) & point.x <= max(x.bounds) & 
             point.y >= min(y.bounds) & point.y <= max(y.bounds))
  }
  
  degToCm = function(angle, distance) {
    return(distance*tan(angle*pi/180))
  }
  
  degToPix = function(angle, distance, resolution, screenSize) {
    return(degToCm(angle, distance) * resolution / screenSize)
  }
  
  cmToPix = function(cm, resolution, screenSize) {
    return(cm * resolution / screenSize)
  }
  
  cmToDeg = function(cm, distance) {
    return(tan(cm/distance)*180/pi)
  }
  
  # Calculate distance to monitor for deflected fixations
  distfix <- function(x,y,distance) {
    distcalc <- sqrt(x^2+y^2+distance^2)
    return(distcalc)
  }
  
  # Calculate angle between 2 vectors (to compute scanpath)
  angle <- function(x,y){
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    theta.deg <- theta/pi*180
    as.numeric(theta.deg)
  }
  
  plotRoi = function(roi, matrixRoi=T) {
    #plot(raster::raster(roi)) #error
    
    # raster = reshape2::melt()
    # raster$value = as.factor(raster$value)
    # with(raster, plot(Var2, Var1, col=ifelse(value==1, "green", "white"), asp=1))
    
    if (matrixRoi) #image function expects x then y but matrix is y then x => switch manually
      roi = roi %>% apply(2, rev) %>% t()
    roi %>% ifelse(1, 0) %>% image(col=c("white", "green"), useRaster=T, asp=1)
  }
  
  #creates new variables that are accessible in the global environment (see "<<-")
  validateBaselines = function(fixs, mess, exclusions, 
                               maxDeviation_rel, maxSpread,
                               saveBaselinePlots=FALSE, prefix="", postfix="") {
    if (saveBaselinePlots) dir.create(paste0(path.rds, "BL plots/"), showWarnings=F)
    
    vpn = fixs$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
    vpn = vpn[!(vpn %in% exclusions)] #minus a priori exclusions
    vpn.n = length(vpn)
    
    trials.n = mess %>% group_by(subject) %>% summarize(trials = max(trial) - min(trial) + 1) %>% .$trials %>% max()
    
    baselines.trial = vector("list", length(vpn.n)) #list of evaluations of baselines per trial for every subject
    baselines.summary = data.frame(subject=character(vpn.n), ntrials=numeric(vpn.n), nValid=numeric(vpn.n), 
                                   invalid=numeric(vpn.n), na=numeric(vpn.n), 
                                   sd_x=numeric(vpn.n), sd_y=numeric(vpn.n),
                                   range_x=numeric(vpn.n), range_y=numeric(vpn.n),
                                   stringsAsFactors=FALSE)
    baselines.summary[,] = NA
    
    print("Validating baselines")
    for (vpIndex in 1:vpn.n) { #TODO make inner of loop enclosed function and call it with apply?
      code = vpn[vpIndex]
      cat(code, " ... ")
      
      # Determine trial number
      fix.vp = fixs %>% filter(subject==code)
      msg.vp = mess %>% filter(subject==code)
      trials.min = min(fix.vp$trial)
      trials.max = max(fix.vp$trial)
      #trials.n = trials.max - trials.min + 1 #must be same for every subject
      
      # baseline
      baseline.vp = data.frame(x=numeric(trials.n), y=numeric(trials.n), condition=character(trials.n), stringsAsFactors=FALSE)
      baseline.vp[,] = NA
      
      # Loop over trials to determine trial-by-trial baselines
      for (trial in fix.vp$trial %>% unique()) {
        # Select trial data
        fix.vp.trial = fix.vp[fix.vp$trial==trial,] #fix.vp %>% filter(trial==trial) #doesn't work
        msg.vp.trial = msg.vp[msg.vp$trial==trial,] #msg.vp %>% filter(trial==trial)
        
        # Determine onset (in ms)
        include <- grep(expoID, msg.vp.trial$event)
        if (length(include) > 0) {
          MsgIndex = grep(expoID, msg.vp.trial$event)
          onset <- msg.vp.trial$time[MsgIndex] %>% as.numeric()
          condition = msg.vp.trial$event[MsgIndex] %>% as.character(); condition = condition %>% substring(nchar(expoID)+1, nchar(condition)-4)
          
          # Subtract onset from startamps
          fix.vp.trial$start  <- fix.vp.trial$start - onset
          fix.vp.trial$end <- fix.vp.trial$end - onset
          
          # Caluculate baseline as weighted average of fixations
          fix.vp.trial.bl <- fix.vp.trial[fix.vp.trial$end>min(baseline) & fix.vp.trial$start<max(baseline),]
          if (nrow(fix.vp.trial.bl) > 0) {
            # Restrict fixation data to baseline
            fix.vp.trial.bl$start[1] <- max(c(min(baseline), first(fix.vp.trial.bl$start)))
            fix.vp.trial.bl$end[nrow(fix.vp.trial.bl)] <- min(c(max(baseline), last(fix.vp.trial.bl$end))) #ifelse(tail(fix.vp.trial.bl$end,1)>blen,blen,tail(fix.vp.trial.bl$end,1)) # = min(blen, tail...) ?
            fix.vp.trial.bl$dur <- fix.vp.trial.bl$end - fix.vp.trial.bl$start
            
            # Calculate baseline coordinates
            xbl <- sum(fix.vp.trial.bl$x*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
            ybl <- sum(fix.vp.trial.bl$y*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
            
            # Store values
            baseline.vp[trial-trials.min+1,] = c(xbl,ybl, condition)
          } else {
            # When no valid fixations are available store NA as baseline for current trial
            baseline.vp[trial-trials.min+1, "condition"] = condition
          }
        }
      }
      
      baseline.vp = baseline.vp %>% mutate(x = as.numeric(x), y = as.numeric(y))
      
      # Determine outlier
      if (usePointDistance) {
        baseline.vp = baseline.vp %>% mutate(
          blok = outlier_remove_point(x, y, maxDeviation_rel) %>%
            outlier_remove_point_spread(x, y, maxSpread, insample=.))
      } else {
        blxok = outlier_remove(baseline.vp$x, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$x, maxSpread, .)
        blyok = outlier_remove(baseline.vp$y, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$y, maxSpread, .)
        baseline.vp$blok = blxok & blyok #Baseline is valid when x and y coordinates are ok (i.e. no outlier)
      }
      
      # Store number of valid baselines per subject
      nValid = sum(baseline.vp$blok)
      invalid = 1 - nValid / trials.n
      nas = sum(is.na(baseline.vp$x)) / trials.n
      x.coords = baseline.vp %>% filter(blok) %>% .$x; mean_x = mean(x.coords); sd_x = sd(x.coords); range_x = x.coords %>% range() %>% diff()
      y.coords = baseline.vp %>% filter(blok) %>% .$y; mean_y = mean(y.coords); sd_y = sd(y.coords); range_y = y.coords %>% range() %>% diff()
      baselines.summary[vpIndex, 2:9] = c(trials.n, nValid, invalid, nas, sd_x, sd_y, range_x, range_y); baselines.summary[vpIndex , 1] = code
      baselines.trial[[code %>% as.character()]] = baseline.vp
      
      #plot subject
      if (saveBaselinePlots==TRUE || code %in% saveBaselinePlots || as.character(code) %in% saveBaselinePlots) {
        filename = paste0(path.rds, "BL plots/", prefix, code, postfix, ".png")
          
        borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
        borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
        borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                          x.coords %>% range() %>% mean() + maxSpread/2)
        borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                          y.coords %>% range() %>% mean() + maxSpread/2)
        
        if (driftPlots==T || code %in% driftPlots || as.character(code) %in% driftPlots) {
          blplot = baseline.vp %>% mutate(trial=1:n()) %>% 
            ggplot(aes(x=x, y=y, color=trial)) + 
            xlim(0, screen.width) + ylim(0, screen.height) + #restrict area to screen
            geom_rect(xmin=0, ymin=0, xmax=screen.width, ymax=screen.height, color="black", fill=NA) +
            geom_hline(yintercept=mean_y, linetype="longdash") + geom_vline(xintercept = mean_x, linetype="longdash") #centroid of valid baselines
          
          #add borders depending on whether point distance was used (circles) or not (rectangles)
          if (usePointDistance) {
            distances = pointDistance(x.coords, y.coords, mean_x, mean_y)
            r_abs = mean(distances) + maxDeviation_rel*sd(distances)
            
            blplot = blplot + #circles
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=r_abs), color="red", fill=NA, inherit.aes=F) + #borders for validity (relative)
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=maxSpread/2), color="orange", fill=NA, inherit.aes=F) #borders for validity (absolute)
          } else {
            borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
            borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
            borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                              x.coords %>% range() %>% mean() + maxSpread/2)
            borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                              y.coords %>% range() %>% mean() + maxSpread/2)
            
            blplot = blplot + #rectangles
              geom_rect(xmin=min(borders.rel.x), ymin=min(borders.rel.y), xmax=max(borders.rel.x), ymax=max(borders.rel.y), color="red", fill=NA) + #borders for validity (relative)
              geom_rect(xmin=min(borders.abs.x), ymin=min(borders.abs.y), xmax=max(borders.abs.x), ymax=max(borders.abs.y), color="orange", fill=NA) #borders for validity (absolute)
          }
          
          #add points
          blplot = blplot + 
            geom_point() + scale_color_continuous(low="blue", high="green") + #all points (color coded by trial)
            geom_point(data=baseline.vp %>% filter(blok==F), mapping=aes(x=x, y=y), color="red") + #invalid baselines red
            geom_point(x=screen.width/2, y=screen.height/2, shape="+", size=5, color="black") + #fixation cross
            #theme(panel.border = element_rect(color = "black", fill=NA, size=5)) +
            coord_fixed() +
            ggtitle(paste0(prefix, code, postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          blplot %>% ggsave(filename=filename, plot=., device="png", dpi=300, units="in", width=1920/300, height = 1080/300)
        } else {
          borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
          borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
          borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                            x.coords %>% range() %>% mean() + maxSpread/2)
          borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                            y.coords %>% range() %>% mean() + maxSpread/2)
          
          #x11() #plot in new windows (max of 63)
          png(filename, 
              width=1920, height=1080)
          plot(baseline.vp$x, baseline.vp$y, pch=16, col=ifelse(baseline.vp$blok==0, "red", "black"), xlab="x (px)", ylab="y (px)", xlim=c(0, screen.width), ylim=c(0, screen.height), asp=1)
          title(paste0(prefix, code, postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          abline(h=c(mean_y, screen.height), v=mean_x); abline(v=borders.rel.x, h=borders.rel.y, col="red"); points(x=mean(c(0, screen.width)), y=screen.height/2, pch=3, col="blue")
          
          dev.off()
        }
      }
      
      #invisible(readline(prompt="Baseline created! Press [enter] to continue"))
    }
    
    baselines.trial <<- baselines.trial
    
    print("Baseline validation finished")
    return(baselines.summary)
  }
  
  roiAnalysis = function(fixs, baselines, showTrialPlots=F, offset=0) {
    trials.vp = fixs %>% .$trial %>% unique() %>% length()
    
    eye.vp = data.frame(trial = numeric(trials.vp),
                        condition = numeric(trials.vp),
                        dwell = numeric(trials.vp),
                        dwell.non = numeric(trials.vp),
                        ms = numeric(trials.vp),
                        ms.non = numeric(trials.vp),
                        diagnosticFirst = numeric(trials.vp), 
                        fixN = numeric(trials.vp),
                        mFixTime = numeric(trials.vp),
                        roiSwitch = numeric(trials.vp),
                        scanPath = numeric(trials.vp))
    for (i in 1:(length(bins)-1)) eye.vp = eye.vp %>% mutate(!!paste0("dwell.bin", i) := numeric(trials.vp))
    for (i in 1:(length(bins)-1)) eye.vp = eye.vp %>% mutate(!!paste0("dwell.non.bin", i) := numeric(trials.vp))
    
    bl.mean.x = mean(baselines$x[baselines$blok]) #calculate mean valid baseline for trials that have invalid baseline
    bl.mean.y = mean(baselines$y[baselines$blok])
    
    #for (i in fixs %>% .$trial %>% unique() %>% sort() %>% {. - offset}) {
    for (i in fixs %>% pull(trial) %>% unique() %>% seq()) {
      
      #print(i)
      
      #baseline correction
      trial = baselines[i, ]
      bl.corr.x = ifelse(is.na(trial$blok)==F && trial$blok, trial$x, bl.mean.x) %>% round(digits=1) #round to accuracy of eye tracker
      bl.corr.y = ifelse(is.na(trial$blok)==F && trial$blok, trial$y, bl.mean.y) %>% round(digits=1)
      
      i.total = i + offset
      
      #fixations
      fixations.trial = fixs %>% filter(trial==i.total)
      if (fixations.trial %>% nrow() == 0) {
        eye.vp[i, 1] = i.total; eye.vp[i, -1] = NA
        next
      }
      
      #messages
      condition = fixations.trial %>% transmute(condition = pic %>% sub(".png", "", ., fixed=T)) %>% .$condition %>% head(1)
      
      #rois
      roi = rois[[condition %>% substring(1, 1) %>% as.integer()]]
      roi.non = rois.nondiag[[condition %>% substring(1, 1) %>% as.integer()]]
      
      fixations.trial = fixations.trial %>% 
        mutate(x.cent = x - bl.corr.x,
               y.cent = y - bl.corr.y,
               dist = sqrt(x.cent^2 + y.cent^2), #distance from center
               x = x - bl.corr.x + screen.width/2, 
               y = y - bl.corr.y + screen.height/2,
               x.roi = {x - screen.width/2 + 1 + dim(roi)[2]/2} %>% round(), #old coordinate minus half screen + 1 => middle of screen & roi = 1. this + half roi => 1 = most left column of roi
               y.roi = {screen.height - y - screen.height/2 + 1 + dim(roi)[1]/2} %>% round(), #same as with x but reverse y-direction first (currently up = bigger but down = bigger needed for matrix indexing)
               #inRoi = ifelse(x.roi > 0 & y.roi > 0 & x.roi <= dim(roi)[2] & y.roi <= dim(roi)[1], roi[y.roi, x.roi], FALSE), #if coordinates within bounds of roi matrix, look up in matrix, else FALSE #doesn't work because whole vector is evaluated and throws error even though result FALSE shall be used
               start = start - onset, #start of exposition = time_0
               end = end - onset, #start of exposition = time_0
        )
      
      #roi analysis
      fixations.trial = fixations.trial %>% mutate(
        inRoi = x.roi > 0 & y.roi > 0 & x.roi <= dim(roi)[2] & y.roi <= dim(roi)[1],
        inRoi.non = x.roi > 0 & y.roi > 0 & x.roi <= dim(roi.non)[2] & y.roi <= dim(roi.non)[1]) #fixation on stimulus? (enables indexing without out of bounds error)
      #are fixations on stimuli also within roi?
      fixations.trial$inRoi[fixations.trial$inRoi] = roi[cbind(fixations.trial$y.roi[fixations.trial$inRoi], fixations.trial$x.roi[fixations.trial$inRoi])]
      fixations.trial$inRoi.non[fixations.trial$inRoi.non] = roi.non[cbind(fixations.trial$y.roi[fixations.trial$inRoi.non], fixations.trial$x.roi[fixations.trial$inRoi.non])]
      fixations.trial = fixations.trial %>% mutate(roi = case_when(
        inRoi ~ "diagnostic", inRoi.non ~ "non-diagnostic", TRUE ~ "no ROI") %>% as.factor())
      
      fixations.trial.analysis = fixations.trial %>% filter(end > 0, start < ratingStart) %>% 
        filter(start > 0) #filter out fixations starting before stimulus onset (even if they extend into a stimulus ROI)
      
      if (fixations.trial.analysis %>% nrow() == 0) {
        eye.vp[i, 1] = i.total; eye.vp[i, -1] = NA
        next
      }
      
      fixations.trial.analysis = fixations.trial.analysis %>% 
        mutate(start = ifelse(start < 0, 0, start),
               end = ifelse(end > ratingStart, ratingStart, end),
               dur = end - start,
               x.st.mm = lag(x.cent) * pixsize, y.st.mm = lag(y.cent) * pixsize, #previous fixation = saccade start
               x.en.mm = x.cent * pixsize, y.en.mm = y.cent * pixsize, #current fixation = saccade end
               x.st.mm = ifelse(x.st.mm %>% is.na(), 0, x.st.mm), x.en.mm = ifelse(x.en.mm %>% is.na(), 0, x.en.mm), #set first fixation to fixation cross
               y.st.mm = ifelse(y.st.mm %>% is.na(), 0, y.st.mm), y.en.mm = ifelse(y.en.mm %>% is.na(), 0, y.en.mm), #set first fixation to fixation cross
               distfix.st = distfix(x.st.mm, y.st.mm, distance),
               distfix.en = distfix(x.en.mm, y.en.mm, distance)
        ) %>% filter(dur > 0) %>% 
        rowwise() %>% mutate(angle = angle(c(x.st.mm, y.st.mm, distfix.st), c(x.en.mm, y.en.mm, distfix.en)))
      
      fixations.trial.analysis.bins = fixations.trial %>% filter(end > min(bins), start < max(bins)) %>% 
        filter(start > 0) %>% #filter out fixations starting before stimulus onset (even if they extend into a stimulus ROI)
        mutate(start = ifelse(start < min(bins), min(bins), start),
               end = ifelse(end > max(bins), max(bins), end),
               dur = end - start) %>% filter(dur > 0)
      
      #output variables
      dwell = with(fixations.trial.analysis, sum(dur[inRoi] / sum(dur)))
      dwell.non = with(fixations.trial.analysis, sum(dur[inRoi.non] / sum(dur[inRoi | inRoi.non])))
      ms = with(fixations.trial.analysis, start[inRoi] %>% head(1)); if (is_empty(ms)) ms = ratingStart
      ms.non = with(fixations.trial.analysis, start[inRoi.non] %>% head(1)); if (is_empty(ms.non)) ms.non = ratingStart
      diagnosticFirst = ms < ms.non
      fixN = fixations.trial.analysis %>% nrow()
      mFixTime = fixations.trial.analysis %>% pull(dur) %>% mean()
      roiSwitch = with(fixations.trial.analysis, roi[roi != "no ROI"]) %>% as.numeric() %>% diff() %>% {. != 0} %>% sum()
      scanPath = with(fixations.trial.analysis, sum(angle, na.rm=T))
      
      dwell.bins = c()
      dwell.non.bins = c()
      for (b in 2:length(bins)) {
        bin = bins[(b-1):b]
        fix.bin = fixations.trial.analysis.bins %>% filter(end >= min(bin) & start <= max(bin)) %>% #only fixations that touch current bin
          mutate(start = ifelse(start <= min(bin), min(bin), start), 
                 end = ifelse(end >= max(bin), max(bin), end),
                 dur = end - start) #prune values to borders of bin to calculate duration correctly
        dwell.bins = with(fix.bin, sum(dur[inRoi] / sum(dur))) %>% c(dwell.bins, .)
        dwell.non.bins = with(fix.bin, sum(dur[inRoi.non] / sum(dur[inRoi | inRoi.non]))) %>% c(dwell.non.bins, .)
      }
      
      if (showTrialPlots) {
        picPath = list.files(paste0(path.rois, ".."), "png", full.names=T); picPath = picPath[picPath %>% grep(condition, .)]
        #picPath = list.files(paste0(path.rois, ".."), "jpg", full.names=T)[condition%/%10 * 2-1 + round((condition%%10-1)/6)]
        fixationCross = head(fixations.trial.analysis, 1) %>% mutate(start = 0, end = 0, x = screen.width/2, y = screen.height/2, roi="no ROI", dur=0)
        fixations.trial.plot = rbind(fixationCross, fixations.trial.analysis) %>% mutate(roi = factor(roi, levels=c("diagnostic", "non-diagnostic", "no ROI")))
        { fixations.trial.plot %>% ggplot(aes(x=x, y=y)) + 
            annotation_custom(grid::rasterGrob(
              png::readPNG(picPath),
              #jpeg::readJPEG(picPath),
              width=unit(1,"npc"), height=unit(1,"npc")), 
              xmin = screen.width/2-642/2, xmax = screen.width/2+642/2, 
              ymin = screen.height/2-676/2, ymax = screen.height/2+676/2) +
            geom_point(aes(x=screen.width/2, y=screen.height/2), color="blue", shape="+", size=5) + #fixation cross
            #geom_point(aes(x=x[2], y=y[2], size=dur[2]), shape=21) + #first fixation after fix cross accentuated
            geom_point(aes(size=dur, color=roi), alpha=.25) + geom_path() + #all fixations as transparent circles with fill
            #geom_point(aes(size=dur, color=inRoi), shape=21) + geom_path() + #all fixations as opaque circles without fill with path
            #geom_rect(mapping=aes(xmin=0, xmax=screen.width, ymin=0, ymax=screen.height), fill=NA, color="black") + #screen border
            #xlim(0, screen.width) + ylim(0, screen.height) + #full screen
            scale_x_continuous(limits=c(screen.width/2-642/2, screen.width/2+642/2), expand=c(0,0)) +
            scale_y_continuous(limits=c(screen.height/2-676/2, screen.height/2+676/2), expand=c(0,0)) +
            coord_fixed() + ggtitle(paste0(code, ": trial ", i.total, ", condition ", condition)) + 
            scale_size(range = c(3, 10)) +
            labs(size = "duration") +
            scale_color_manual(values=c("green3", "red", "grey")) +
            theme_bw() + theme(plot.title = element_text(hjust = 0.5)) } %>% print()
      }
      
      #assign output variables
      eye.vp[i, -2] = c(i.total, dwell, dwell.non, ms, ms.non, diagnosticFirst, fixN, mFixTime, roiSwitch, scanPath, dwell.bins, dwell.non.bins); eye.vp[i, 2] = condition
    }
    
    return(eye.vp)
  }
}



# Load Data ---------------------------------------------------------------
#behavior.valid = read_rds("behavior.valid.rds" %>% paste0(path.rds, .))
messages = "Messages.txt" %>% paste0(path.eye, .) %>% 
  read_delim(delim="\t", locale=locale(decimal_mark=","), na=".", show_col_types=F) %>% #read.table(filePath, skip=1, dec=",", sep="\t", na.strings=".")
  rename(subject = RECORDING_SESSION_LABEL, trial = TRIAL_LABEL, time = CURRENT_MSG_TIME, event = CURRENT_MSG_TEXT) %>% 
  separate("subject", c("subject", "block"), sep = "_") %>% 
  mutate(block = block %>% as.numeric(), 
         trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
         trial = trial + (block-1)*trials.N/2) %>% 
  separate(event, c("left", "right"), sep = " \\| ", remove = F) %>% 
  separate(left, c("distract.left", "target.left"), sep = " ") %>% 
  separate(right, c("distract.right", "target.right"), sep = " ") %>% 
  mutate(angry = if_else(distract.left %>% grepl("category1", .), "left", "right"),
         #angryCheck = if_else(distract.right %>% grepl("category1", .), "right", "left"),
         target.left = target.left %>% na_if("NA"), target.right = target.right %>% na_if("NA"),
         #targetCheck = if_else(target.right %>% is.na(), "left", "right"),
         target = case_when(subject %>% str_starts("b") ~ NA, #two targets for dual probe
                            target.left %>% is.na() ~ "right", 
                            target.right %>% is.na() ~ "left",
                            T ~ "error")
  ) #%>% filter(angry != angryCheck | target != targetCheck | target == "error")

#messages %>% count(subject) %>% filter(n != trials.N)
messages %>% count(subject) %>% filter(n != trials.N, subject %in% {behavior.overview %>% filter(noET %>% is.na() == F) %>% pull(subject)} == F)
#a07: no 2nd block (known from other variables)
#problem resolved: [a19: last 1 trial  of 1st block missing (but not in reaction times)]
#problem resolved: [b04: last 4 trials of 2nd block missing (but not in reaction times)]
#a19 & b04 amended by writing messages that occur without eye-tracking data (eyes were missing for so long such that the eye-tracker went into standby)
#View(messages %>% filter(subject %in% {messages %>% count(subject) %>% filter(n != trials.N) %>% pull(subject) %>% setdiff(behavior.overview %>% filter(noET %>% is.na() == F) %>% pull(subject))}), "messageCheck")
#all missings were manually compared against sequence files to assure trial number is correct

#messages %>% pull(time) %>% summary()
messages %>% pull(time) %>% quantile(c(0, 0.001, 0.002, 0.003, 0.004), na.rm=T) #very small percentage of weird message timing
messages %>% filter(time < 900)
#messages %>% filter(time < 900) %>% pull(time) %>% hist()
#Checked the raw data: Before problematic message, eyes had not been detected for longer time
# => tracker went into stand-by (?) and did not write further timestamps with missing data
# ==> message timing tells how long before stimulus onset, valid eye-tracking data are existent
#     (e.g., -300 => 300 ms AFTER stimulus onset)

fixations = "Fixations.txt" %>% paste0(path.eye, .) %>% 
  read_delim(delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F) %>% #read.table(filePath, skip=1, dec=",", sep="\t", na.strings=".")
  #fixations = fixations %>% select(-2) #drop 2nd column containing only "Trial: " #only works if sep-parameter is set to default
  rename(subject = X1, trial = X2, start = X3,end = X4, x = X5, y = X6) %>% 
  separate("subject", c("subject", "block"), sep = "_") %>% 
  mutate(block = block %>% as.numeric(), 
         trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
         trial = trial + (block-1)*trials.N/2,
         x = x %>% gsub(",", ".", .) %>% na_if("NA") %>% as.numeric(),
         y = y %>% gsub(",", ".", .) %>% na_if("NA") %>% as.numeric(),
         y = screen.height - y) %>% 
  #filter(subject %in% exclusions.eye.num == F) %>% #exclusion is done in baseline validation below
  left_join(messages %>% select(subject, block, trial, time) %>% rename(onset = time),
            by=c("subject", "block", "trial")) %>% 
  left_join(behavior, by=c("subject", "block", "trial")) %>% #get conditions & other variables
  arrange(subject, trial)

fixations.distractors = fixations %>% 
  mutate(start = start - onset, end = end - onset, #realign such that 0 = stim start
         start = ifelse(start < 0, 0, start), #discard fraction of fixation before distractor onset
         end = if_else(end > SOA, SOA, end), #discard fraction of fixation after distractor offset
         dur = end - start) %>% 
  filter(dur > 0) #note: this may drop some trials completely

#TODO remove trials with invalid responses?
#fixations.valid %>% filter(response %>% is.na()) %>% select(subject, trial) %>% unique()
#fixations.valid %>% filter(angry %>% is.na() | congruency %>% is.na()) %>% select(subject, trial, congruency, angry, response_dual) %>% unique() %>% print(n = nrow(.))


# Valid Fixations ---------------------------------------------------------
vpn.eye = fixations.distractors %>% pull(subject) %>% unique() %>% setdiff(exclusions.eye) %>% sort()

#subjects with 0 valid fixations during any distractor phase of the experiment
fixations %>% pull(subject) %>% unique() %>% setdiff(vpn.eye) %>% sort()

#determine valid trials by fixation time (invalid trials need not be evaluated by their baseline)
eye.valid.trial = fixations.distractors %>% 
  summarize(.by = c(subject, block, trial, SOA), valid = sum(dur)/first(SOA), n = n())

eye.valid.trial %>% mutate(block = block %>% as_factor()) %>% 
  ggplot(aes(x = valid, fill = block)) +
  geom_histogram(breaks=seq(0, 1, length.out=20+1), color="black") +
  geom_vline(xintercept = validFixTime.trial, color = "red", linetype = "dashed", linewidth = 2) +
  scale_fill_viridis_d() + myGgTheme

eye.valid.trial %>% group_by(subject) %>% 
  summarize(valid.n = sum(valid > validFixTime.trial),
            valid.p = valid.n / trials.N) %>% 
  summarize(n = sum(valid.p > validFixTime.subj),
            n_ex = sum(valid.p <= validFixTime.subj) + {fixations %>% pull(subject) %>% unique() %>% setdiff(vpn.eye) %>% length()},
            retention = n / (n + n_ex))

#check if missingness increases over time (e.g., fatigue => droopy eyes)
rToFishZ = function(r) { return(.5*log((1+r)/(1-r))) }
fishZtoR = function(Z) { return(ifelse(Z==Inf, 1, (exp(2*Z)-1)/(exp(2*Z)+1))) }
fnFishZ = function(r, fn=mean, prune=.999, ...) { 
  warning.txt = ""
  if (any(abs(r) >= 1, na.rm=T)) {
    warning.txt = "Correlation(s) contain(s) values at or beyond 1. This will bias results for summary functions." %>% paste0(warning.txt, .)
    if (prune %>% is.na()) warning = "Consider setting prune parameter." %>% paste(warning.txt, .)
    else {
      warning.txt = paste0("Pruning to ", prune, ".") %>% paste(warning.txt, .)
      
      r = r %>% prune(abs=prune) #if_else(abs(r) >= 1, prune * if_else(r > 0, 1, -1), r)
    }
    warning(warning.txt)
  }
  
  return(r %>% rToFishZ() %>% fn(...) %>% fishZtoR()) 
}

#eye.valid.trial %>% count(subject, block) %>% filter(n <= 3)
eye.trial_decline.simple = eye.valid.trial %>% 
  filter(subject %in% {eye.valid.trial %>% count(subject, block) %>% filter(n <= 3) %>% pull(subject) %>% c("b28")} == F) %>% 
  mutate(trial = trial - (block-1)*trials.N.block) %>% #restart trial count at block start (due to recalibration of eye tracker)
  summarize(.by = c(subject, block),
            #TODO better use binomial regression
            cor = cor(trial, valid),
            p = cor.test(trial, valid)$p.value, 
            cor.test = cor.test(trial, valid) %>% apa::cor_apa(r_ci=T, print=F)) %>% 
  arrange(p)
#eye.trial_decline.simple %>% filter(p < .05) %>% print(n=nrow(.)) #all significant correlations are negative
eye.trial_decline.simple %>% summarize(cor.m = fnFishZ(cor), cor.sd = fnFishZ(cor, fn=sd))

eye.valid.trial %>% 
  ggplot(aes(x = trial, y = valid, 
             color = subject, 
             group = interaction(subject, block))) +
  #geom_point(alpha=.1) +
  stat_smooth(method="lm") +
  myGgTheme

#beta mixed-effects regression
library(glmmTMB)
library(broom.mixed)

eye.trial_decline <- glmmTMB(
  valid ~ SOA * trial + (trial | subject),
  data = eye.valid.trial %>% 
    mutate(valid = (valid * (n() - 1) + 0.5) / n(), #correction for 0 and 1 values (not allowed)
           trial = trial - (block-1)*trials.N.block), #restart trial count at block start (due to recalibration of eye tracker
  family = beta_family(link = "logit")
)
tidy(eye.trial_decline, effects = "fixed", conf.int = TRUE)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
# fixations.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% select(subject, block, trial) %>% 
#   left_join(fixations, by=c("subject", "trial", "block"))
fixations.valid = fixations #take all fixations


# baseline validation --------------------------------------------------------------------
baselines.summary1 = validateBaselines(fixations.valid %>% filter(block==1), messages %>% filter(block==1), exclusions.eye, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="_1") %>% mutate(block=1) %>% select(subject, block, everything())
baselines.trial1 = baselines.trial; rm(baselines.trial)
baselines.summary2 = validateBaselines(fixations.valid %>% filter(block==2), messages %>% filter(block==2), exclusions.eye, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="_2") %>% mutate(block=2) %>% select(subject, block, everything())
baselines.trial2 = baselines.trial; rm(baselines.trial)

# baseline validation summary
baselines.summary = baselines.summary1 %>% bind_rows(baselines.summary2) %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread,
         block = block %>% as_factor())
rm(baselines.summary1, baselines.summary2)

baselines.summary %>% ggplot(aes(x = invalid, fill = block)) + 
  geom_histogram(boundary=outlierLimit.eye, binwidth=.025, color="black") + 
  geom_vline(xintercept = outlierLimit.eye, color = "red", linetype = "dashed", linewidth = 2) +
  scale_fill_viridis_d() + myGgTheme

baselines.summary %>% summarize(.by = block, totalN = n(), includedN = sum(included), includedP = mean(included))

baselines.summary %>% summarize(.by = subject, invalid = mean(invalid)) %>% arrange(desc(invalid))
baselines.summary %>% group_by(subject) %>% summarize(invalid = mean(invalid)) %>% 
  filter(invalid < outlierLimit.eye) %>% summarize(invalid.m = mean(invalid), invalid.sd = sd(invalid))


# ROI analysis ------------------------------------------------------------
fixations.valid = fixations.valid %>% 
  mutate(ROI.side = case_when(x <= screen.width/2 - indifference.zone ~ "left",
                              x >= screen.width/2 + indifference.zone ~ "right",
                              T ~ NA),
         ROI = if_else(ROI.side == angry, "angry", "neutral"))
  
# fixations.valid %>% ggplot(aes(x = x, y = y, group = interaction(subject, trial), color = ROI.side)) +
#   geom_rect(xmax=screen.width/2-indifference.zone, xmin=-Inf, ymin=-Inf, ymax=+Inf, fill="blue", alpha=.2) +
#   geom_rect(xmin=screen.width/2+indifference.zone, xmax=+Inf, ymin=-Inf, ymax=+Inf, fill="red", alpha=.2) +
#   geom_point(alpha=.1) + #geom_line(alpha=.1) +
#   myGgTheme
fixations.valid %>% ggplot(aes(x = x, fill = ROI.side)) +
  geom_histogram(boundary = screen.width/2, binwidth = indifference.zone/2, color = "black") +
  #geom_vline(xintercept = c(-indifference.zone, +indifference.zone) + screen.width/2, linetype = "dashed", color = "red") +
  myGgTheme

fixations.first = fixations.valid %>% 
  filter(start > 0, 
         #end < SOA, #100 ms SOA is too short for a saccade => also count saccades after face offset (despite strong confound with target)
  ) %>% 
  summarize(.by = c(subject, trial, paradigm, SOA, angry),
            firstDwell = ROI %>% na.omit() %>% first()) %>% 
  summarize(.by = -c(trial, firstDwell),
            firstDwell.angry = mean(firstDwell=="angry", na.rm=T)) %>% 
  arrange(subject, paradigm, SOA, angry) %>% mutate(SOA = SOA %>% as_factor())


# Statistical Analysis ----------------------------------------------------
fixations.first.missing = fixations.first %>% filter(firstDwell.angry %>% is.na()) %>% pull(subject) %>% unique() %>% sort()
fixations.first %>% pull(subject) %>% unique() %>% setdiff(fixations.first.missing) %>% length() %>% paste0("N = ", .)

fixations.first %>% 
  filter(subject %in% fixations.first.missing == F) %>% 
  mutate(firstDwell.angry = firstDwell.angry - .5) %>% #interpretability of intercept
  ez::ezANOVA(dv = firstDwell.angry, wid = subject,
              within = .(SOA, angry),
              between = paradigm,
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)

# paradigm x angry
fixations.first %>% 
  filter(subject %in% fixations.first.missing == F) %>% 
  filter(paradigm == "Dot Probe") %>% 
  mutate(firstDwell.angry = firstDwell.angry - .5) %>% #interpretability of intercept
  ez::ezANOVA(dv = firstDwell.angry, wid = subject,
              within = .(SOA, angry),
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)

fixations.first %>% 
  filter(subject %in% fixations.first.missing == F) %>% 
  filter(paradigm == "Dual Probe") %>% 
  mutate(firstDwell.angry = firstDwell.angry - .5) %>% #interpretability of intercept
  ez::ezANOVA(dv = firstDwell.angry, wid = subject,
              within = .(SOA, angry),
              type=2, detailed=T) %>% apa::anova_apa(force_sph_corr = T)


fixations.first %>% summarize(.by = -c(subject, firstDwell.angry),
                              firstDwell.angry.se = se(firstDwell.angry, na.rm=T),
                              firstDwell.angry = mean(firstDwell.angry, na.rm=T)) %>% 
  ggplot(aes(x = angry, y = firstDwell.angry, color = SOA)) +
  facet_wrap(~paradigm) +
  geom_hline(yintercept = .5, color = "black", linetype = "dashed") +
  geom_errorbar(aes(ymin = firstDwell.angry - firstDwell.angry.se, ymax = firstDwell.angry + firstDwell.angry.se)) +
  geom_point() + 
  myGgTheme
