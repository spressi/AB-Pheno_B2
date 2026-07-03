library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")

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
         paradigm = if_else(subject %>% grepl("a", .), "Discrimination", "Localization"))

eeg.markers %>% count(subject) %>% filter(n != markers.n) %>% arrange(n)
#a03: EEG recording started too late, first 4 EOG calibration markers (3 trials) missing => use 2nd EOG for both blocks?

#check missing markers
eeg.markers %>% count(subject, value) %>% 
  mutate(correct = case_when(value %in% c(88, 99) ~ T, #always assume correct here / check later
                             value == 211 & n == 4 ~ T, #EOG center
                             value >= 200 & n == 2 ~ T, #all other EOG
                             n %in% c(8, 16) ~ T, #could differentiate: a needs 8 for target markers, b 16 throughout
                             T ~ F)) %>% 
  filter(correct == F) %>% 
  mutate(meaning = case_when(value == 200 ~ "Block Start", value == 201 ~ "Block End",
                             value == 210 ~ "EOG Start", value == 211 ~ "EOG Center", value == 212 ~ "EOG Top", value == 213 ~ "EOG Right", value == 214 ~ "EOG Bottom", value == 215 ~ "EOG Left", value == 220 ~ "EOG End",
                             value %% 10 <= 2 ~ "distractor",
                             value %% 10 >= 4 ~ "target"))
#note: target markers will be skipped if a response is premature, cf.
#behavior %>% filter(expositionCheck %>% is.na())

#check missing trials (should stand out on other variables, too)
eeg.markers %>% filter(value %in% c(88, 99)) %>% count(subject, value) %>% pivot_wider(names_from = value, values_from = n) %>% mutate(sum = `88` + `99`) %>% filter(sum != trials.N)


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
  #filter(samplediff == min(samplediff, na.rm=T)) %>% 
  filter(samplediff < 50) %>% 
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
#TODO renew calculation of premature responses in 1.1 Behavior?

eeg.markers.long %>% select(-contains("sample")) %>% 
  filter(value_response == 99) %>% 
  filter(stimToResp == min(stimToResp, na.rm=T) |
           stimToResp == max(stimToResp, na.rm=T) |
           stimToNextResp == min(stimToNextResp, na.rm=T)
  ) %>% arrange(rt) %>% relocate(stimToResp, stimToNextResp)
#max time to correct response: TODO ms
#min time to NEXT correct response: TODO ms
# => check for overlap
  
# Impedances --------------------------------------------------------------
files.eeg.headers = list.files(path.eeg.raw, pattern = ".vhdr", full.names = T)
eeg.impedances.list = list()
for (file in files.eeg.headers) {
  #file = files.eeg.headers %>% sample(1) #for testing
  
  skip = 111 #start with skipping 111 lines (will be reduced for each subject)
  repeat { #do-while loop
    checkFile = file %>% 
      read_table(skip = skip, col_names = c("electrode", "impedance"), show_col_types = F, na = "???")
    if (checkFile %>% pull(electrode) %>% grepl("Fp1", .) %>% any()) { #check if Fp1 is contained (first electrode)
      break
    } else {
      skip = skip - 1
      if (skip < 0) stop(paste0("No start of impedance list found for ", file %>% pathToCode(), ". Check vhdr file if increasing skip parameter will solve the issue: \n", file))
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
eeg.impedances %>% filter(impedance %>% is.na()) %>% pull(subject_session) %>% unique()

##impedance change before/after
eeg.impedances.m = eeg.impedances %>% summarize(.by = c(subject, time),
                                                impedance = mean(impedance, na.rm=T)) %>% 
  pivot_wider(names_from = time, values_from = impedance)

with(eeg.impedances.m, t.test(after, before, paired=T)) %>% apa::t_apa(es_ci=T)
eeg.impedances.m %>% summarize(before.m = mean(before), before.sd = sd(before), 
                               after.m = mean(after), after.sd = sd(after))
