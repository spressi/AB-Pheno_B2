library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")

writeCorrectedMarkers = T #rewrite marker file for subjects with inverted markers (low voltage = signal instead of high)

files.eeg.markers = list.files(path.eeg.raw, pattern = ".vmrk", full.names = T) %>% 
  Filter(\(x) x %>% grepl("_2", .) == F, .) %>% #get rid of second file
  Filter(\(x) x %>% grepl("_original", .) == F, .) #get rid of original marker files (just backup for transparency)
eeg.markers.list = list()
for (file in files.eeg.markers) {
  #file = files.eeg.markers %>% sample(1) #for testing
  eeg.markers.list[[pathToCode(file)]] = file %>% read.csv(skip=11, header=F, col.names = c("marker", "value", "sample", "size", "channel"))
}
#tidy up
eeg.markers = eeg.markers.list %>% bind_rows(.id = "subject") %>% tibble() %>% 
  filter(marker %>% grepl("Stimulus", .)) %>% 
  mutate(value = value %>% gsub("S\\s*", "", .) %>% as.integer(),
         paradigm = if_else(subject %>% grepl("a", .), "Dot Probe", "Dual Probe"))
eeg.markers %>% count(subject) %>% filter(n != 1152) %>% arrange(n)


# correct for inverted markers --------------------------------------------
invertedMarkers = eeg.markers %>% count(subject) %>% filter(n > 1152) %>% pull(subject) #some subjects had inverted markers (low voltage = signal instead of high). Every reset of the pins thus resulted in an additional marker
#eeg.markers %>% filter(subject %in% invertedMarkers) #255 is marker stop => recode as 255 - value
eeg.markers = eeg.markers %>% 
  mutate(value = if_else(subject %in% invertedMarkers, 255 - value, value)) %>% 
  filter(value != 0) %>% #get rid of stop markers
  #check assumption: first marker is always "Mk2" (last is "Mk2304")
  #summarize(.by = subject, marker.min = first(marker), marker.max = last(marker)) %>% filter(marker.min != "Mk2=Stimulus" | marker.max != "Mk2304=Stimulus", subject %in% invertedMarkers)
  mutate(.by = subject, markerNum = 1:n()+1)
#check assumption: only subjects with inverted markers are affected by wrong marker order
#eeg.markers %>% mutate(markerCheck = paste0("Mk", markerNum, "=Stimulus")) %>% filter(markerCheck != marker) %>% pull(subject) %>% unique() %>% {. == invertedMarkers} %>% all()

#check inverted markers
#eeg.markers %>% filter(subject %in% invertedMarkers)



# correct for Dual Probe --------------------------------------------------
# adding 5 to markers for SOA==100 messed up everything for Dual Probe
# targets 1 & 2 will get code 12, which will be mapped to 17 for SOA==100 => can't distinguish between 12 for SOA==100 and 17 for SOA==500
# targets 9 & 4 will get code 94 (for angry left) - if this happens within SOA==100, we get 99, which will get confused with a correct response
# => get rid of target information in markers and just discern between angry face position (1 vs. 2) and soa (0 vs. 10)

#add trial numbers (note: don't do this before having corrected inverted markers!)
eeg.markers2 = eeg.markers %>% 
  filter(markerNum %% 2 == 0) %>% #first marker starts with 2 => we need all even markers
  mutate(.by = subject, trial = 1:n()) %>% 
  left_join(eeg.markers, .) %>% #insert trial numbers (but only for stimuli, not for responses!)
  left_join(behavior %>% select(subject, trial, SOA)) %>% #insert SOA from behavior (only for stimuli, not for responses!)
  
  mutate(value = case_when(paradigm == "Dual Probe" & SOA==100 ~ value-5, #revert the +5 from SOA==100
                           T ~ value),
  #) %>% filter(paradigm == "Dual Probe") %>% count(value) %>% arrange(value) %>% data.frame()
         value = case_when(paradigm == "Dual Probe" & trial %>% is.na() == F & value < 100 ~ 1, #angry left  := 1 (i.e., drop information on targets)
                           paradigm == "Dual Probe" & trial %>% is.na() == F & value > 100 ~ 2, #angry right := 2 (i.e., drop information on targets)
                           T ~ value), #leave responses and Dot Probe values untouched
         value = case_when(paradigm == "Dual Probe" & SOA==100 ~ value+5, #now add +5 for SOA==100
                           T ~ value)) %>% 
  #fill(trial, .direction = "down") %>% 
  
  #create output that can be written into a file
  mutate(marker = paste0("Mk", markerNum, "=Stimulus"), #this also overwrites correct markers but result has been checked to be the same
         output = paste(marker, paste0("S", value), sample, size, channel, sep = ","))
# check if Dot Probe has been affected
#eeg.markers2 %>% rename(value2 = value) %>% full_join(eeg.markers) %>% filter(value!=value2) %>% pull(subject) %>% unique() %>% Filter(\(x) x %>% str_starts("a"), .)
eeg.markers = eeg.markers2; rm(eeg.markers2)

eeg.markers %>% count(subject) %>% filter(n != 1152) %>% arrange(n)
#a07: only first block
#a13: last 2 markers (i.e., last trial) missing in EEG (verified, see code below)
#a30: last 2 markers (i.e., last trial) missing in EEG (verified, see code below)

# #check what happened to a13
# list.files(path.seq, pattern = "a13", full.names = T) %>% Filter(\(x) x %>% grepl("_0", .) == F, .) %>% #get rid of training
#   lapply(read_tsv) %>% bind_rows() %>% 
#   mutate(trial = 1:n(), condition = condition + if_else(soa == 100, 5, 0)) %>% 
#   select(trial, condition) %>% 
#   filter(condition != eeg.markers %>% filter(subject == "a13", value %in% c(99, 96) == F) %>% pull(value) %>% c(NA))
# #=> no mismatch across trials

# #check what happened to a30
# list.files(path.seq, pattern = "a30", full.names = T) %>% Filter(\(x) x %>% grepl("_0", .) == F, .) %>% #get rid of training
#   lapply(read_tsv) %>% bind_rows() %>%
#   mutate(trial = 1:n(), condition = condition + if_else(soa == 100, 5, 0)) %>%
#   select(trial, condition) %>%
#   filter(condition != eeg.markers %>% filter(subject == "a30", value %in% c(99, 96) == F) %>% pull(value) %>% c(NA))
# #=> no mismatch across trials


# Breaks ------------------------------------------------------------------
breakTime = 2000 #in samples
eeg.markers.breaks = eeg.markers %>% 
  mutate(.by = subject, timediff = lead(sample)-sample) %>% 
  summarize(.by = subject, 
            breaks = sum(timediff > breakTime, na.rm=T),
            breakIndex = if_else(breaks == 1, which.max(timediff), NA),
            n = n()
            )
eeg.markers.breaks %>% filter(breaks != 1 | breakIndex != trials.N | n != trials.N*2) #2 markers per trial (distractors & response) => break should be at trials.N * 2 / 2 = trials.N
#break always after 576th marker


# Write Files -------------------------------------------------------------
if (writeCorrectedMarkers) {
  markerFilesToWrite = c(invertedMarkers, eeg.markers %>% pull(subject) %>% unique() %>% Filter(\(x) x %>% str_starts("b"), .)) %>% unique() %>% sort()
  for (s in markerFilesToWrite) {
    #s = sample(markerFilesToWrite, 1) #for testing
    filename = files.eeg.markers %>% Filter(\(x) x %>% grepl(s, .), .)
    filename.copy = filename %>% gsub(".vmrk", "_original.vmrk", ., fixed = T)
    if (file.exists(filename.copy)) {
      message(paste0(s, ": Original file already exists. Skipping creation of adjusted marker file."))
      next
    }
    
    if (file.exists(filename.copy)==F) #careful! there is no option that prevents file.rename from overwriting an existing file => check yourself
      file.rename(filename, filename.copy) #this retains "last modified" as the original file creation date
    
    file = readLines(filename.copy)
    #file[12] #assert that line 12 is the last line that should remain unmodified: Mk1=New Segment
    file = c(file[1:12],
             eeg.markers %>% filter(subject == s) %>% pull(output))
    writeLines(file, filename)
    cat(paste0(s, ": Adjusted marker file created."))
  }
}
