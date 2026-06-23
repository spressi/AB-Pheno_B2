library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")

markers.n = trials.N*3 + # 3 triggers per trial (distractors, target, response)
  2*2 + #start & end triggers for each block
  (6 + 2)*2 #EOG calibration: 6 positions (2x middle, top, right, bottom, left) + start & end triggers; for each block

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
#eeg.markers %>% count(subject, value) %>% View()

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
