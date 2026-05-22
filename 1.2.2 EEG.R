library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")
#source("1.2.1 EEG Markers.R")

# Markers -----------------------------------------------------------------
files.eeg.markers = list.files(path.eeg.raw, pattern = ".mrk", full.names = T) %>% 
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

# assert correct markers
eeg.markers %>% filter(paradigm=="Dot Probe" & value %in% c(11, 12, 16, 17, 21, 22, 26, 27, 111, 112, 116, 117, 121, 122, 126, 127, 96, 99) == F | 
                         paradigm=="Dual Probe" & value %in% c(1, 2, 6, 7, 88, 99) == F)
eeg.markers %>% filter(paradigm=="Dot Probe") %>% count(value) %>% arrange(value)
eeg.markers %>% filter(paradigm=="Dual Probe") %>% count(value) %>% arrange(value)

# # assert balancing of angry left & right
# sequences %>% filter(subject %>% str_starts("b")) %>% 
#   filter(subject %in% {eeg.markers %>% pull(subject) %>% unique()}) %>% 
#   #pull(subject) %>% unique()
#   mutate(angry = if_else(distractor_left %>% grepl("AN", .), "left", "right")) %>% 
#   count(angry)

# assert correct timing
eeg.markers %>% 
  mutate(.by = subject, 
         trial = ceiling(1:n()/2),
         samplediff = sample-lag(sample)) %>% 
  filter(samplediff == min(samplediff, na.rm=T)) %>% 
  left_join(sequences %>% select(subject, trial, SOA, iti)) #premature response within the first 100 ms => shorter than a normal 100 ms trial :(

eeg.markers.long = eeg.markers %>% 
  mutate(.by = subject, trial = ceiling(1:n()/2)) %>% 
  mutate(.by = c(subject, trial), helper = 1:n()) %>% 
  mutate(kind = if_else(helper == 1, "stim", "response")) %>% select(-helper) %>% 
  pivot_wider(names_from = kind, values_from = c(value, sample), id_cols = c(subject, trial)) %>% 
  left_join(behavior %>% select(subject, block, trial)) %>% relocate(subject, block, trial) %>% #insert block
  left_join(sequences %>% select(subject, trial, SOA, iti)) %>% #insert SOA & iti
  mutate(stimToResp = (sample_response - sample_stim) / hz.eeg * 1000,
         rt = stimToResp - SOA) %>% 
  mutate(.by = c(subject, block), 
         respToNextStim = (lead(sample_stim) - sample_response) / hz.eeg * 1000,
         stimToNextStim = (lead(sample_stim) - sample_stim) / hz.eeg * 1000,
         stimToNextResp = (lead(sample_response) - sample_stim) / hz.eeg * 1000)
eeg.markers.long %>% select(-contains("sample")) %>% 
  filter(#.by = subject,
    stimToResp == min(stimToResp, na.rm=T) |
      stimToResp == max(stimToResp, na.rm=T) |
      #respToNextStim == min(respToNextStim, na.rm=T) |
      #respToNextStim == max(respToNextStim, na.rm=T) |
      stimToNextStim == min(stimToNextStim, na.rm=T) |
      stimToNextStim == max(stimToNextStim, na.rm=T)
  ) %>% arrange(rt) #%>% relocate(stimToResp, stimToNextStim)

eeg.markers.long %>% select(-contains("sample")) %>% 
  filter(value_response == 99) %>% 
  filter(stimToResp == min(stimToResp, na.rm=T) |
           stimToResp == max(stimToResp, na.rm=T) |
           stimToNextResp == min(stimToNextResp, na.rm=T)
  ) %>% arrange(rt) %>% relocate(stimToResp, stimToNextResp)
#max time to correct response: 2010 ms
#min time to NEXT correct response: 2024 ms
# => barely not overlapping!
  
# Impedances --------------------------------------------------------------
files.eeg.headers = list.files(path.eeg.raw, pattern = ".vhdr", full.names = T)
eeg.impedances.list = list()
for (file in files.eeg.headers) {
  #file = files.eeg.headers %>% sample(1) #for testing
  
  skip = 112 #start with skipping 112 lines (b06 for whatever reason)
  repeat { #do-while loop
    checkFile = file %>% 
      read_table(skip = skip, col_names = c("electrode", "impedance"), show_col_types = F, na = "???")
    if (checkFile %>% pull(electrode) %>% grepl("Fp1", .) %>% any()) { #check if Fp1 is contained (first electrode)
      break
    } else {
      skip = skip - 1
    }
  }
  
  eeg.impedances.list[[pathToCode(file)]] = checkFile
}
#tidy up
eeg.impedances = eeg.impedances.list %>% bind_rows(.id = "subject") %>% 
  mutate(electrode = electrode %>% gsub(":", "", .),
         time = if_else(subject %>% grepl("_2", .), "after", "before") %>% as_factor(),
         subject_session = subject,
         subject = subject %>% gsub("_2", "", .))


#sanity checks
eeg.impedances %>% count(electrode, name = "files") %>% count(files, name = "electrodes")
eeg.impedances %>% filter(impedance %>% is.na()) %>% pull(subject_session) %>% unique() #TODO check NAs

##impedance change before/after
# with(eeg.impedances %>% summarize(.by = c(subject, time),
#                              impedance = mean(impedance, na.rm=T)) %>% 
#   pivot_wider(id_cols = subject, names_from = time, values_from = impedance),
#   t.test(before, after, paired=T)
# ) %>% apa::t_apa(es_ci=T)
#just saving another file does not work to measure impedances after experiment :(
# => need to click on impedance measurement again; otherwise previous values will just get carried over

eeg.impedances.m = eeg.impedances %>% summarize(.by = c(subject, time),
                                                impedance = mean(impedance, na.rm=T)) %>% 
  pivot_wider(names_from = time, values_from = impedance) %>% filter(before != after)

with(eeg.impedances.m, t.test(after, before, paired=T)
) %>% apa::t_apa(es_ci=T)
eeg.impedances.m %>% summarize(before.m = mean(before), before.sd = sd(before), 
                               after.m = mean(after), after.sd = sd(after))
