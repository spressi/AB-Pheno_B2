library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")

writeCorrectedMarkers = F #rewrite marker file for subjects with inverted markers (low voltage = signal instead of high)
markers.n = trials.N*3 #TODO add extra markers (block & EOG calibration)
#TODO adjust computations by extra markers! save number in separate variable?

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

eeg.markers %>% count(subject) %>% filter(n != markers.n) %>% arrange(n)


# Breaks ------------------------------------------------------------------
#TODO will be obsolete when adding start & finish marker
breakTime = 2000 #in samples
eeg.markers.breaks = eeg.markers %>% 
  mutate(.by = subject, timediff = lead(sample)-sample) %>% 
  summarize(.by = subject, 
            breaks = sum(timediff > breakTime, na.rm=T),
            breakIndex = if_else(breaks == 1, which.max(timediff), NA),
            n = n()
            )
eeg.markers.breaks %>% filter(breaks != 1 | breakIndex != markers.n/2 | n != markers.n) #3 markers per trial (distractors, target, & response) => break should be at trials.N * 3 / 2
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
