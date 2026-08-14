library(tidyverse)
#source("0 General.R")
#source("1.1 Behavior.R")

# Markers -----------------------------------------------------------------
markers.n = trials.N*3 + # 3 triggers per trial (distractors, target, response)
  2*2 + #start & end triggers for each block
  (6 + 2)*2 #EOG calibration: 6 positions (2x middle, top, right, bottom, left) + start & end triggers; for each block

files.eeg.markers = list.files(path.eeg.raw, pattern = ".mrk", full.names = T) %>% 
  Filter(\(x) x %>% str_detect("_2") == F, .) %>% #get rid of second file
  Filter(\(x) x %>% str_detect("_original") == F, .) #get rid of original marker files (just backup for transparency)
eeg.markers.list = list()
for (file in files.eeg.markers) {
  #file = files.eeg.markers %>% sample(1) #for testing
  eeg.markers.list[[pathToCode(file)]] = file %>% read.csv(skip=11, header=F, col.names = c("marker", "value", "sample", "size", "channel"))
}

#tidy up
eeg.markers = eeg.markers.list %>% bind_rows(.id = "subject") %>% tibble() %>% 
  filter(marker %>% str_detect("Stimulus")) %>% 
  mutate(value = value %>% str_replace("S\\s*", "") %>% as.integer(),
         kind = case_when(value >= 200 ~ "extra",
                          value %in% c(88, 99) ~ "response",
                          value %% 10 <= 2 ~ "distractor",
                          value %% 10 >= 4 ~ "target"),
         paradigm = if_else(subject %>% str_detect("a"), "Discrimination", "Localization"))

#trials
eeg.markers = eeg.markers %>% 
  filter(kind == "distractor") %>% #count trials using distractor triggers (which never turned out missing)
  mutate(.by = subject, trial = 1:n()) %>% 
  select(subject, sample, trial) %>% full_join(eeg.markers, ., join_by(subject, sample)) %>% 
  #group_by(subject) %>% fill(trial, .direction = "down") %>% ungroup() %>% 
  fill(trial, .direction = "down") %>% mutate(trial = if_else(kind == "extra", NA, trial))


# * Missing Markers -------------------------------------------------------
eeg.markers %>% count(subject) %>% filter(n != markers.n) %>% mutate(diff = n - markers.n) #%>% arrange(n)
#a03: first 4 EOG calibration markers (3 trials) missing => use 2nd EOG for both blocks?
#b04: first EOG start marker missing (no problem) + 1 response missing???
#b06: 2 targets missing + 1 response missing???
#b07: 2 targets missing

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
eeg.markers %>% filter(value %in% c(88, 99)) %>% count(subject, value) %>% pivot_wider(names_from = value, values_from = n) %>% mutate(sum = `88` + `99`) %>% filter(sum != trials.N) %>% mutate(diff = sum - trials.N)


# * Assertions ------------------------------------------------------------
# # assert balancing of angry left & right
# sequences %>% filter(subject %>% str_starts("b")) %>% 
#   filter(subject %in% {eeg.markers %>% pull(subject) %>% unique()}) %>% 
#   #pull(subject) %>% unique()
#   mutate(angry = if_else(distractor_left %>% str_detect("AN"), "left", "right")) %>% 
#   count(angry)

# assert correct timing
eeg.markers %>% 
  mutate(.by = subject, 
         trial = ceiling(1:n()/2),
         samplediff = sample-lag(sample)) %>% 
  #filter(samplediff == min(samplediff, na.rm=T)) %>% 
  filter(samplediff < 50) %>% 
  left_join(sequences %>% select(subject, trial, SOA, iti)) #premature response within the first 100 ms => shorter than a normal 100 ms trial :(


# * Long Format -----------------------------------------------------------
eeg.markers.wide = eeg.markers %>% 
  filter(trial %>% is.na() == F) %>% 
  pivot_wider(names_from = kind, values_from = c(value, sample), id_cols = c(subject, trial)) %>% 
  left_join(behavior %>% select(subject, block, trial), by = join_by(subject, trial)) %>% relocate(subject, block, trial) %>% #insert block
  left_join(sequences %>% select(subject, trial, SOA, iti), by = join_by(subject, trial)) %>% #insert SOA & iti
  mutate(stimToResp = (sample_response - sample_distractor) / hz.eeg * 1000,
         rt = stimToResp - SOA) %>% 
  mutate(.by = c(subject, block), 
         respToNextStim = (lead(sample_distractor) - sample_response) / hz.eeg * 1000,
         stimToNextStim = (lead(sample_distractor) - sample_distractor) / hz.eeg * 1000,
         stimToNextResp = (lead(sample_response) - sample_distractor) / hz.eeg * 1000)
eeg.markers.wide %>% select(-starts_with("value_"), -starts_with("sample_")) %>% 
  filter(#.by = subject,
    stimToResp == min(stimToResp, na.rm=T) |
      stimToResp == max(stimToResp, na.rm=T) |
      #respToNextStim == min(respToNextStim, na.rm=T) |
      #respToNextStim == max(respToNextStim, na.rm=T) |
      stimToNextStim == min(stimToNextStim, na.rm=T) |
      stimToNextStim == max(stimToNextStim, na.rm=T)
  ) %>% arrange(rt) #%>% relocate(stimToResp, stimToNextStim)
#TODO renew calculation of premature responses in 1.1 Behavior?

eeg.markers.wide %>% 
  filter(value_response == 99) %>% 
  select(-starts_with("value_"), -contains("sample")) %>% 
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
    if (checkFile %>% pull(electrode) %>% str_detect("Fp1") %>% any()) { #check if Fp1 is contained (first electrode)
      break
    } else {
      skip = skip - 1
      if (skip < 0) stop(paste0("No start of impedance list found for ", file %>% pathToCode(), ". Check vhdr file if increasing skip parameter will solve the issue: \n", file))
    }
  }
  
  eeg.impedances.list[[pathToCode(file)]] = checkFile %>% 
    mutate(impedance = if_else(impedance=="Out", Inf, impedance %>% as.integer()))
}
#tidy up
eeg.impedances = eeg.impedances.list %>% bind_rows(.id = "subject") %>% 
  mutate(electrode = electrode %>% str_replace(":", ""),
         time = if_else(subject %>% str_detect("_2"), "after", "before") %>% as_factor(),
         subject_session = subject,
         subject = subject %>% str_replace("_2", ""))


#sanity checks
eeg.impedances %>% count(electrode, name = "files") %>% count(files, name = "electrodes")
eeg.impedances %>% filter(impedance %>% is.na()) %>% pull(subject_session) %>% unique()
eeg.impedances %>% filter(impedance == Inf)

##impedance change before/after
eeg.impedances.m = eeg.impedances %>% 
  filter(impedance != Inf) %>% 
  summarize(.by = c(subject, time),
            impedance = mean(impedance, na.rm=T)) %>% 
  pivot_wider(names_from = time, values_from = impedance)

with(eeg.impedances.m, t.test(after, before, paired=T)) %>% apa::t_apa(es_ci=T)
eeg.impedances.m %>% summarize(before.m = mean(before), before.sd = sd(before), 
                               after.m = mean(after), after.sd = sd(after))
