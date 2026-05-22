###############################################################################
# Project "Attentional Bias Phenotypes" of Mario Reutter
# Script by Matthias Gamer
#
# Before running, convert idf-files to text (*_Samples.txt) using IDF converter contained in iTools_setup.exe
# 
# Eyetracking-Daten in Event-Listen umwandeln
# - Getrennte Listen fuer Fixation, Saccades und Messages im SR Research Format erstellen

###############################################################################
# Data formats:
# 
# MESSAGES:
# RECORDING_SESSION_LABEL	TRIAL_LABEL	CURRENT_MSG_TIME	CURRENT_MSG_TEXT
# vpja01a	Trial: 1	1063	Stimulus 5D.jpg
# vpja01a	Trial: 2	1050	Stimulus 3B.jpg
# 
# FIXATIONS:
# RECORDING_SESSION_LABEL	TRIAL_LABEL	CURRENT_FIX_START	CURRENT_FIX_END	CURRENT_FIX_X	CURRENT_FIX_Y
# vpsa01	Trial: 1	0	128	892,10	564,20
# vpsa01	Trial: 1	156	320	962,00	551,50
# vpsa01	Trial: 1	348	1172	913,30	572,00
# vpsa01	Trial: 1	1192	1428	892,20	537,30
#
# SACCADES:
# RECORDING_SESSION_LABEL	TRIAL_LABEL	CURRENT_SAC_START_TIME	CURRENT_SAC_END_TIME	CURRENT_SAC_START_X	CURRENT_SAC_END_X	CURRENT_SAC_START_Y	CURRENT_SAC_END_Y	CURRENT_SAC_CONTAINS_BLINK
# vp01l	Trial: 1	833	840	916,70	928,00	581,50	583,50	false
# vp01l	Trial: 1	1076	1083	924,20	912,40	586,60	583,80	false
# vp01l	Trial: 1	1639	1665	907,90	889,10	588,80	497,80	false
# vp01l	Trial: 1	1817	1840	880,20	806,20	492,60	489,60	false
###############################################################################

library(tidyverse)
#source("0 General.R")
distance = 650 #screen center to subject eyes in mm
screen.width.mm = 520 #screen width in mm
screen.width.px = 1920 #screen width in pixels
screen.height.px = 1200 #screen height in pixels (only for plot)

options(warn=2)  # Halt on warnings

helpers.dir = "eyeScoring/"
source(paste0(helpers.dir, "eventdetection_lowspeed.R"))
source(paste0(helpers.dir, "eventdetection_highspeed.R"))
source(paste0(helpers.dir, "targetCols.R"))

pathToCode = function(path, path.sep="/", file.ext="\\.") {
  first = path %>% gregexpr(path.sep, .) %>% sapply(max) %>% {. + 1}
  last = path %>% gregexpr(file.ext, .) %>% sapply(max) %>% {. - 1}
  return(path %>% substring(first, last))
}
codeToNum = function(code) code %>% gsub("\\D+", "", .) %>% as.integer()
toPhysString = function(x, width=2) formatC(x, width=width, format="d", flag="0") #add leading zeros (and convert to string)
toPhysFileName = function(x, width=2, prefix="vp", postfix=".txt") x %>% toPhysString() %>% paste0(prefix, ., postfix)
toCode = function(x, prefix="vp", postfix="") x %>% toPhysString() %>% paste0(prefix, ., postfix)

# One degree of visual angle amounts to how many pixels?
winpix <- distance * tan(pi/180) * screen.width.px/screen.width.mm

vpn <- list.files(path.eye.raw, pattern=".txt", full.names = T)
if (exists("allfixtab")) rm("allfixtab"); if (exists("allsactab")) rm("allsactab"); if (exists("allmsgtab")) rm("allmsgtab")
for (vp in vpn) {
  #vp = head(vpn, 1)
  #vp = sample(vpn, 1)
  
  #vp = vpn %>% Filter(\(x) {x %>% grepl("a05_2", .)}, .)
  #vp = vpn %>% Filter(\(x) {x %>% grepl("a19_2", .)}, .) #problem resolved
  #vp = vpn %>% Filter(\(x) {x %>% grepl("b04_2", .)}, .) #problem resolved
  
  code <- vp %>% pathToCode(file.ext = " ")
  print(code)
  
  
  ############################################
  # Read data for definition of trials
  daten <- scan(vp,what="character",skip=38,quiet=TRUE)
  
  # Find Messages
  imsg <- which(daten=="MSG")
  # Trial definition
  #msgtxt <- daten[imsg+5] #for calib
  msgtxt <- paste(daten[imsg+4], daten[imsg+5], daten[imsg+6], daten[imsg+7], daten[imsg+8]) #for main experiment... note to future self: don't use whitespaces in messages :)
  #if (length(msgtxt) != 16) { #for calib
  if (length(msgtxt) != trials.N/2) { #for main experiment
    print(paste0(code, ": ERROR Trial labels"))
  }
  
  ############################################
  # Read data for scoring events
  #daten <- read.table(vp,skip=38,colClasses=c("numeric","character",rep("numeric",11)),fill=TRUE)
  daten <- read.table(vp,skip=38,fill=TRUE) %>% tibble() %>% select(all_of(targetIndeces))
  names(daten) <- c("Time","Type","Trial","PupilL_x","pupilL_y","PupilL_mm","PupilR_x","pupilR_y","PupilR_mm", 
                    "PosL_x","PosL_y","PosR_x","PosR_y")
  
  # Find Messages
  imsg <- which(daten$Type=="MSG")
  
  # Process trialwise data (Timestamps are in nanoseconds)
  for (trial in 1:length(imsg)) {
    st <- daten$Time[imsg[trial]]-prestim*1000
    en <- daten$Time[imsg[trial]]+max(poststim)*1000
    
    tdaten <- daten[(daten$Time>=st) & (daten$Time<=en),]
    
    # Remove Messages
    tdaten <- tdaten[tdaten$Type=="SMP",]
    
    if (nrow(tdaten)>0) {
      # Create new timestamps in ms
      onset  <- round((daten$Time[imsg[trial]]-tdaten$Time[1])/1000,1)
      tdaten$Time <- round((tdaten$Time-tdaten$Time[1])/1000,1)
      
      # Interpolate missing data rows
      imissing <- as.numeric(which((apply(tdaten[,4:ncol(tdaten)],1,sum)==0) & 
                                     (c(1,diff(tdaten$PupilL_x))!=0) & 
                                     (c(diff(tdaten$PupilL_x),1)!=0)))
      imissing <- imissing[!(imissing %in% c(1,nrow(tdaten)))]
      
      tdaten[imissing,4:13] <- (tdaten[imissing-1,4:13]+tdaten[imissing+1,4:13])/2
      
      # Compute average of L and R eye (if one is 0 take the other)
      tdaten$Pos_x <- (tdaten$PosL_x+tdaten$PosR_x)/2
      tdaten$Pos_x[tdaten$PosL_x==0] <- tdaten$PosR_x[tdaten$PosL_x==0]
      tdaten$Pos_x[tdaten$PosR_x==0] <- tdaten$PosL_x[tdaten$PosR_x==0]
      tdaten$Pos_y <- (tdaten$PosL_y+tdaten$PosR_y)/2
      tdaten$Pos_y[tdaten$PosL_y==0] <- tdaten$PosR_y[tdaten$PosL_y==0]
      tdaten$Pos_y[tdaten$PosR_y==0] <- tdaten$PosL_y[tdaten$PosR_y==0]
      
      # Detect blinks and NA gaze position
      iblink <- which(tdaten$PupilL_mm==0 | tdaten$PupilR_mm==0)
      tdaten$blink <- 0; tdaten$blink[iblink] <- 1
      # Add 2 Samples before and after blink for accurate gaze detection
      tdaten$blink <- as.numeric((tdaten$blink + c(tdaten$blink[3:nrow(tdaten)],0,0)
                                  + c(0,0,tdaten$blink[1:(nrow(tdaten)-2)]))>0)
      tdaten$Pos_x[tdaten$blink==1] <- NA
      tdaten$Pos_y[tdaten$blink==1] <- NA
      
      # Determine fixations and saccades
      #timeline <- tdaten$Time
      #x <- tdaten$Pos_x
      #y <- tdaten$Pos_y
      #blk <- tdaten$blink
      #pixperdeg <- winpix
      #plotdata <- TRUE
      
      events <- eventdetection_highspeed(tdaten$Time, tdaten$Pos_x, tdaten$Pos_y, tdaten$blink, winpix, plotdata=FALSE, trial)
      #events <- eventdetection_highspeed(tdaten$Time, tdaten$Pos_x, tdaten$Pos_y, tdaten$blink, winpix, plotdata=TRUE, trial)
      #readline(paste("Trial:",trial))
      
      # Store fixations
      if (!is.na(events$fixations)[1]) {
        if (exists("allfixtab")) {
          allfixtab <- rbind(allfixtab,data.frame(#code=rep(paste(vp,"s",session,sep=""),nrow(events$fixations)),
                                                  code=rep(code,nrow(events$fixations)),
                                                  trial=rep(paste("Trial:",trial),nrow(events$fixations)),
                                                  events$fixations))
        } else {
          allfixtab <- data.frame(#code=rep(paste(vp,"s",session,sep=""),nrow(events$fixations)),
                                  code=rep(code,nrow(events$fixations)),
                                  trial=rep(paste("Trial:",trial),nrow(events$fixations)),
                                  events$fixations)
        }
      }
      
      # Store saccades
      if (!is.na(events$saccades)[1]) {
        if (exists("allsactab")) {
          allsactab <- rbind(allsactab,data.frame(code=rep(code,nrow(events$saccades)),
                                                  #code=rep(paste(vp,"s",session,sep=""),nrow(events$saccades)),
                                                  trial=rep(paste("Trial:",trial),nrow(events$saccades)),
                                                  events$saccades))
        } else {
          allsactab <- data.frame(code=rep(code,nrow(events$saccades)),
                                  #code=rep(paste(vp,"s",session,sep=""),nrow(events$saccades)),
                                  trial=rep(paste("Trial:",trial),nrow(events$saccades)),
                                  events$saccades)
        }
      }
      
      # Store messages
      if (exists("allmsgtab")) {
        allmsgtab <- rbind(allmsgtab,data.frame(code=code,
                                                #code=paste(vp,"s",session,sep=""),
                                                trial=paste("Trial:",trial),
                                                msgtime=onset,
                                                msgtxt=msgtxt[trial]))
      } else {
        allmsgtab <- data.frame(code=code,
                                #code=paste(vp,"s",session,sep=""),
                                trial=paste("Trial:",trial),
                                msgtime=onset,
                                msgtxt=msgtxt[trial])
      }
    } else {
      # Store messages without valid fixations
      if (exists("allmsgtab")) {
        allmsgtab <- rbind(allmsgtab,data.frame(code=code,
                                                #code=paste(vp,"s",session,sep=""),
                                                trial=paste("Trial:",trial),
                                                msgtime=NA,
                                                msgtxt=msgtxt[trial]))
      } else {
        allmsgtab <- data.frame(code=code,
                                #code=paste(vp,"s",session,sep=""),
                                trial=paste("Trial:",trial),
                                msgtime=NA,
                                msgtxt=msgtxt[trial])
      }
    }
  }
}

# Save results
if (!file.exists(path.eye)) dir.create(path.eye)
names(allmsgtab) <- c("RECORDING_SESSION_LABEL","TRIAL_LABEL","CURRENT_MSG_TIME","CURRENT_MSG_TEXT")
write.table(allmsgtab,paste(path.eye,"Messages.txt",sep=""),sep="\t",dec=",",quote=FALSE,row.names=FALSE)

names(allfixtab) <- c("RECORDING_SESSION_LABEL","TRIAL_LABEL",
                      "CURRENT_FIX_START","CURRENT_FIX_END","CURRENT_FIX_X","CURRENT_FIX_Y")
write.table(allfixtab,paste(path.eye,"Fixations.txt",sep=""),sep="\t",dec=",",quote=FALSE,row.names=FALSE)

# Code numeric to boolean
allsactab$blink <- as.logical(allsactab[,ncol(allsactab)])
allsactab <- allsactab[,c(1:8,10)]
names(allsactab) <- c("RECORDING_SESSION_LABEL","TRIAL_LABEL",
                      "CURRENT_SAC_START_TIME","CURRENT_SAC_END_TIME",
                      "CURRENT_SAC_START_X","CURRENT_SAC_END_X","CURRENT_SAC_START_Y","CURRENT_SAC_END_Y",
                      "CURRENT_SAC_CONTAINS_BLINK")
write.table(allsactab,paste(path.eye,"Saccades.txt",sep=""),sep="\t",dec=",",quote=FALSE,row.names=FALSE)


# Check Eye Calib ---------------------------------------------------------
read_tsv(paste0(path.eye %>% gsub("/Summary", "/test/Summary", .), "Fixations.txt"), locale = locale(decimal_mark = ",")) %>% 
  left_join(read_tsv(paste0(path.eye %>% gsub("/Summary", "/test/Summary", .), "Messages.txt"), locale = locale(decimal_mark = ",")), by = c("RECORDING_SESSION_LABEL", "TRIAL_LABEL")) %>% 
  separate(CURRENT_MSG_TEXT, c("distractL", "targetL", NA, "distractR", "targetR"), sep = " ") %>% 
  mutate(angry = if_else(distractL %>% grepl("category1", .), "left", "right"),
         targetSide = if_else(targetL=="NA", "right", "left"),
         congruency = if_else(angry==targetSide, "congruent", "incongruent"),
         targetKind = if_else(targetSide=="left", targetL, targetR)) %>% 
  #select(-(1:7)) %>% unique() #for checking
  mutate(start = CURRENT_FIX_START - CURRENT_MSG_TIME,
         start = if_else(start < 0, 0, start),
         end = CURRENT_FIX_END - CURRENT_MSG_TIME,
         end = if_else(end < 0, 0, end),
         dur = end - start) %>% 
  filter(dur > 0) %>% 
  rename(subject = RECORDING_SESSION_LABEL) %>% 
  ggplot(aes(x = CURRENT_FIX_X, y = CURRENT_FIX_Y, color = subject, size = dur)) +
  geom_rect(color="black", fill=NA, xmin = 0, xmax = screen.width.px, ymin = 0, ymax = screen.height.px, inherit.aes = F) +
  geom_path(aes(group = interaction(subject, TRIAL_LABEL)), size=1, alpha = .125) +
  geom_point(alpha = .25) + 
  scale_color_viridis_d() + myGgTheme


# Check Fixations ---------------------------------------------------------
read_tsv(paste0(path.eye, "Fixations.txt"), locale = locale(decimal_mark = ",")) %>% 
  left_join(read_tsv(paste0(path.eye, "Messages.txt"), locale = locale(decimal_mark = ",")), by = c("RECORDING_SESSION_LABEL", "TRIAL_LABEL")) %>% 
  separate(CURRENT_MSG_TEXT, c("distractL", "targetL", NA, "distractR", "targetR"), sep = " ") %>% 
  mutate(angry = if_else(distractL %>% grepl("category1", .), "left", "right"),
         targetSide = if_else(targetL=="NA", "right", "left"),
         congruency = if_else(angry==targetSide, "congruent", "incongruent"),
         targetKind = if_else(targetSide=="left", targetL, targetR)) %>% 
  #select(-(1:7)) %>% unique() #for checking
  mutate(start = CURRENT_FIX_START - CURRENT_MSG_TIME,
         start = if_else(start < 0, 0, start),
         end = CURRENT_FIX_END - CURRENT_MSG_TIME,
         end = if_else(end < 0, 0, end),
         dur = end - start) %>% 
  filter(dur > 0) %>% 
  rename(subject = RECORDING_SESSION_LABEL) %>% 
  ggplot(aes(x = CURRENT_FIX_X, y = CURRENT_FIX_Y, color = subject, size = dur)) +
  facet_grid(rows = vars(angry), cols = vars(targetSide), labeller = "label_both") +
  geom_rect(color="black", fill=NA, xmin = 0, xmax = screen.width.px, ymin = 0, ymax = screen.height.px, inherit.aes = F) +
  #geom_path(aes(group = interaction(subject, TRIAL_LABEL)), size=1, alpha = .125) +
  geom_point(alpha = 1/2^6) + guides(color = "none") +
  scale_color_viridis_d() + myGgTheme
