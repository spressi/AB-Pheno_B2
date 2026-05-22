library(tidyverse)
#source("0 General.R")

missing.max = .2 #maximum proportion of missing values to score questionnaire within subject

missingItems.vars = function(vars) return(rowSums(is.na(across(vars))) / ncol(across(vars)))
sumscore.vars = function(vars, missing.max=1) {
  ifelse(missingItems.vars(vars) > missing.max, NA,
         rowMeans(across(vars), na.rm=T) * ncol(across(vars)))
}

missingItems = function(identifier) return(rowSums(is.na(across(starts_with(identifier)))) / ncol(across(starts_with(identifier))))
sumscore = function(identifier, missing.max=1) {
  ifelse(missingItems(identifier) > missing.max, NA,
         rowMeans(across(starts_with(identifier)), na.rm=T) * ncol(across(starts_with(identifier))))
}


# Questionnaires ----------------------------------------------------------
que.raw = read_csv2(path.que) %>% 
  rename(subject = VP01s) %>% 
  select(subject, STARTED,
         starts_with("sias"), #SIAS (social anxiety)
         starts_with("FB02_"), #STAI trait (general anxiety)
         starts_with("SD")) #sociodemographics

#que.raw %>% filter(subject %>% duplicated()) #no restarts of questionnaires
#que.raw %>% filter(if_any(starts_with("sias") | starts_with("FB02_"), is.na))

que = que.raw %>% mutate(
  #items 5, 9, & 11 have already been reversed in SoSciSurvey
  sias = sumscore("sias", missing.max),
  
  #across(c(FB02_01, FB02_06, FB02_07, FB02_10, FB02_13, FB02_14, FB02_16, FB02_19), ~ 5 - .), #STAI trait reverse items
  #items have already been reversed in SoSciSurvey
  stai = sumscore("FB02_", missing.max),
  
  age = SD01 + 17, #dropdown menu, first entry is 18
  
  edu = SD02 %>% replace_values(5 ~ 0) %>% #no finished education first
    #as_factor() %>% #doesn't work...
    recode_values(0 ~ "kein Schulabschluss",
                  1 ~ "Hauptschulabschluss",
                  2 ~ "mittlere Reife",
                  3 ~ "Hochschulreife",
                  4 ~ "abgeschlossenes (Fach-)Hochschul-Studium"), 
  
  work = SD03 %>% 
    #as_factor() %>% #doesn't work...
    recode_values(1 ~ "arbeitssuchend",
                  2 ~ "Schüler:in",
                  3 ~ "in Ausbildung",
                  4 ~ "Student:in",
                  5 ~ "angestellt/selbstständig"), 
  
  gender = SD04 %>% recode_values(1 ~ "male",
                                  2 ~ "female",
                                  3 ~ "non-binary"), 
  
  hand = SD06 %>% recode_values(1 ~ "left-handed",
                                2 ~ "right-handed",
                                3 ~ "ambidextrous"), 
  comment = SD07_01
) %>% select(subject, sias, stai, age:comment)


# Samples Descriptives ----------------------------------------------------
# que %>% summarize(across(c(age, sias, stai),
#                          descriptives.list)) %>%
#   pivot_longer(everything()) %>% separate(name, c("variable", "name")) %>% pivot_wider() #separate may fail depending on variable name
que %>% pivot_longer(c(age, sias, stai)) %>% 
  summarize(.by = name, across(value, descriptives.list, .names = "{.fn}"))

que %>% checkContent(gender)

que %>% checkContent(edu)
que %>% checkContent(work)
que %>% checkContent(hand)

que %>% filter(comment %>% is.na() == F)
