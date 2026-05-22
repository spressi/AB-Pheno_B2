originals = list.files(path.eeg.raw, pattern = "_original.vmrk", full.names = T)
file.rename(originals, originals %>% gsub("_original", "", .))
