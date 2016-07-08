# Plotting script of mate-fidelity
# Luke J. Eberhart-Phillips
# July 8, 2016

library(dplyr)
library(Rmisc)
library(extrafont)
library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(grid)
library(RColorBrewer)

# Find fonts from computer that are menlo or Menlo
font_import(pattern="[M/m]enlo") 
loadfonts(device = "win") # load these into R

# Madagascar data import, tidying, and preprations
# import files
Madagascar_BirdRef <- 
  read.csv("data/Madagascar/Madagascar_BirdRef_2002-2015.csv",
           header=T, stringsAsFactors = FALSE)
Madagascar_Nests <- 
  read.csv("data/Madagascar/Madagascar_Nests_2002-2015.csv",
           header=T, stringsAsFactors = FALSE)
Madagascar_Sexes <- 
  read.csv("data/Madagascar/Madagascar_Sexes_2003-2015.csv",
           header=T)

Madagascar_Sexes <- Madagascar_Sexes[,c("species", "ring", "mol_sex")]
Madagascar_Sexes$mol_sex <- toupper(Madagascar_Sexes$mol_sex)

Andava_Nests <- Madagascar_Nests[((grepl("Andava", 
                                         Madagascar_Nests$site)|
                                     grepl("Sirabe", 
                                           Madagascar_Nests$site))),]
Andava_BirdRef <- Madagascar_BirdRef[((grepl("Andava", 
                                             Madagascar_BirdRef$site)|
                                         grepl("Sirabe", 
                                               Madagascar_BirdRef$site))),]

# subset the data to extract only nests that hatched
#Andava_Nests_Hatched <- dplyr::filter(Andava_Nests, fate == "HATCH")

# subset the data to remove entries that have a "+" by the no_chicks
Andava_Nests <- Andava_Nests[!grepl("[+]",as.character(Andava_Nests$no_chicks)),]

# Make all species entries capitalized and consistent
Andava_Nests$species <- toupper(Andava_Nests$species)
Andava_Nests$species <- ifelse(Andava_Nests$species == "WP", "WFP",
                               ifelse(Andava_Nests$species == "KP", "KIP",
                                      Andava_Nests$species))

# Make all species entries capitalized and consistent
Andava_BirdRef$species <- toupper(Andava_BirdRef$species)
Andava_BirdRef$species <- ifelse(Andava_BirdRef$species == "WP", "WFP",
                               ifelse(Andava_BirdRef$species == "KP", "KIP",
                                      Andava_BirdRef$species))

# define no_chicks as.numeric
Andava_Nests$no_chicks <- as.numeric(as.character(Andava_Nests$no_chicks))

# insert zeros in situations where the fate of a nest was not unknown or hatched
Andava_Nests$no_chicks <- ifelse(Andava_Nests$fate == "UNKNOWN", NA,
                                ifelse(Andava_Nests$fate != "HATCH",
                                       0, Andava_Nests$no_chicks))

# make sure that there are no white spaces hidden in the ID info
Andava_Nests$year <- as.factor(gsub(" ", "", Andava_Nests$year, fixed = TRUE))
Andava_Nests$species <- as.factor(gsub(" ", "", Andava_Nests$species, fixed = TRUE))
Andava_Nests$nest <- as.factor(gsub(" ", "", Andava_Nests$nest, fixed = TRUE))
Andava_BirdRef$year <- as.factor(gsub(" ", "", Andava_BirdRef$year, fixed = TRUE))
Andava_BirdRef$species <- as.factor(gsub(" ", "", Andava_BirdRef$species, fixed = TRUE))
Andava_BirdRef$nest <- as.factor(gsub(" ", "", Andava_BirdRef$nest, fixed = TRUE))

# check the levels of the ID info to see if there are any obivous mistakes
levels(Andava_Nests$year)
levels(Andava_Nests$species)
levels(Andava_Nests$nest)
levels(Andava_BirdRef$year)
levels(Andava_BirdRef$species)
levels(Andava_BirdRef$nest)

# paste the year, site, and ID into one unique ID variable
Andava_Nests$ID <- paste(Andava_Nests$species, Andava_Nests$year, Andava_Nests$nest, 
                         sep = "_")
Andava_BirdRef$ID <- paste(Andava_BirdRef$species, Andava_BirdRef$year, Andava_BirdRef$nest, 
                           sep = "_")

# link the nest ID of the BirdRef file to the nest ID in the nest file.
Andava_Joined <- inner_join(Andava_Nests[,c("species","no_chicks", "clutch_size", "ID", "year")],
                           Andava_BirdRef[,c("parent1", "parent2", "ID")], by = "ID")

# Correct Sexing file
Madagascar_Sexes[which(Madagascar_Sexes$ring == "FH47594"),c("mol_sex")] <- "M"
Madagascar_Sexes[which(Madagascar_Sexes$ring == "FH47634"),c("mol_sex")] <- "M"
Madagascar_Sexes[which(Madagascar_Sexes$ring == "FH69169"),c("mol_sex")] <- "F"
Madagascar_Sexes[which(Madagascar_Sexes$ring == "FH69256"),c("mol_sex")] <- "F"


names(Andava_Joined)[6] <- "ring"
Andava_Joined_Sex <- inner_join(Andava_Joined, Madagascar_Sexes[,c("ring", "mol_sex")], by = "ring")
names(Andava_Joined_Sex)[6] <- "parent1"
names(Andava_Joined_Sex)[8] <- "sex_parent1"
names(Andava_Joined_Sex)[7] <- "ring"
Andava_Joined_Sex <- inner_join(Andava_Joined_Sex, Madagascar_Sexes[,c("ring", "mol_sex")], by = "ring")
names(Andava_Joined_Sex)[7] <- "parent2"
names(Andava_Joined_Sex)[9] <- "sex_parent2"
Andava_Joined_Sex <- Andava_Joined_Sex[,c("species", "no_chicks", "clutch_size",
                                          "ID", "year", "parent1", "parent2",
                                          "sex_parent1", "sex_parent2")]

same_sex <- Andava_Joined_Sex[which(Andava_Joined_Sex$sex_parent1 == Andava_Joined_Sex$sex_parent2),]
same_sex <- same_sex[!duplicated(same_sex),]
#write.csv(same_sex, "Data_files/Andavadoaka_families_with_homosexual_parents.csv",row.names = FALSE)

Andava_Joined_Sex <- anti_join(Andava_Joined_Sex, same_sex, by = "ID")
Andava_Joined_Sex <- Andava_Joined_Sex[!duplicated(Andava_Joined_Sex),]
Andava_Joined_Sex$female <- ifelse(Andava_Joined_Sex$sex_parent1 == "F", Andava_Joined_Sex$parent1,
                                   ifelse(Andava_Joined_Sex$sex_parent2 == "F", Andava_Joined_Sex$parent2, ""))
Andava_Joined_Sex$male <- ifelse(Andava_Joined_Sex$sex_parent1 == "M", Andava_Joined_Sex$parent1,
                                   ifelse(Andava_Joined_Sex$sex_parent2 == "M", Andava_Joined_Sex$parent2, ""))
Andava_Joined_Sex$pair <- paste(Andava_Joined_Sex$female, Andava_Joined_Sex$male, sep = "-")
Andava_Joined_Sex <- Andava_Joined_Sex[which(str_detect(Andava_Joined_Sex$male, "FH") & str_detect(Andava_Joined_Sex$female, "FH")),]
Andava_Joined_Sex$female <- as.factor(Andava_Joined_Sex$female)
Andava_Joined_Sex$male <- as.factor(Andava_Joined_Sex$male)
Andava_Joined_Sex$pair <- as.factor(Andava_Joined_Sex$pair)

df <- Andava_Joined_Sex[, c("species", "year", "ID", "female", "male", "pair")]
df <- df[!duplicated(df),]
# MP FH72484 is monogamous 3 times within 2015 season
females <- dcast(df, species + female  ~ year)
number_males_p_female <- aggregate(male ~ female, df, function(x) length(unique(x)))
number_attempts_p_female <- aggregate(male ~ female, df, function(x) length(duplicated(x)))
test2 <- df[df$female %in% as.vector(number_attempts_p_female[which(number_attempts_p_female$male > 1),1]),]
test2 <- test2[with(test2,order(female,year)), ]

females <- inner_join(females, number_males_p_female)
females[,c(3:8)] <- 
  lapply(females[,c(3:8)], as.numeric)

females$attempts <- rowSums(females[, c(3:8)])
females$years <- rowSums(females[, c(3:8)] > 0)

#females$mate_fidelity <- females$male / females$years
females_no_1 <- filter(females, male  != 1 | years != 1 | attempts != 1)
# Kip FH47071 monogamous between and within seasons
# Kip FH47194 polygamous between seasons
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(2,9)] <- c("focal", "mate")


males <- dcast(df, species + male  ~ year)
number_females_p_male <- aggregate(female ~ male, df, function(x) length(unique(x)))
males <- inner_join(males, number_females_p_male)
males[,c(3:8)] <- 
  lapply(males[,c(3:8)], as.numeric)
males$attempts <- rowSums(males[, c(3:8)])
males$years <- rowSums(males[, c(3:8)] > 0)
#males$mate_fidelity <- males$female / males$years
males_no_1 <- filter(males, female  != 1 | years != 1 | attempts != 1)
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(2,9)] <- c("focal", "mate")

mating_df_Mada <- rbind(females_no_1, males_no_1)
mating_df_Mada$status <- ifelse(mating_df_Mada$mate == 1 & mating_df_Mada$years == mating_df_Mada$attempts, "Monogamous between years",
                           ifelse(mating_df_Mada$mate == 1 & mating_df_Mada$years < mating_df_Mada$attempts, "Monogamous within years",
                                  ifelse(mating_df_Mada$mate > 1 & mating_df_Mada$years == mating_df_Mada$attempts, "Polygamous between years",
                                         ifelse(mating_df_Mada$mate > 1 & mating_df_Mada$years < mating_df_Mada$attempts, "Polygamous within years", "XXX"))))
mating_Mada <- mating_df_Mada[,c("species", "sex", "status", "focal")]

# Ceuta data import, tidying, and preprations
Ceuta_BirdRef <- 
  read.csv("data/Ceuta/Ceuta_BirdRef_2006-2012.csv",
           header=T, stringsAsFactors = FALSE)
Ceuta_Nests <- 
  read.csv("data/Ceuta/Ceuta_Nests_2006-2012.csv",
           header=T, stringsAsFactors = FALSE)
Ceuta_Captures <- 
  read.csv("data/Ceuta/Ceuta_Captures_2006-2012.csv",
           header=T)

Ceuta_Captures$sex <- as.factor(toupper(substr(Ceuta_Captures$age_sex, 2, 2)))
Ceuta_Captures <- Ceuta_Captures[which(Ceuta_Captures$sex == "F" |
                                         Ceuta_Captures$sex == "M"),]
Ceuta_Sexes <- dcast(Ceuta_Captures, ring  ~ sex)
Ceuta_Sexes$sex <- ifelse(Ceuta_Sexes$F > 0 & Ceuta_Sexes$M == 0, "F",
                          ifelse(Ceuta_Sexes$M > 0 & Ceuta_Sexes$F == 0, "M", "Conflict"))
Ceuta_Sexes <- filter(Ceuta_Sexes, sex != "Conflict") [, c(1,4)]

# subset the data to extract only nests that hatched
#Ceuta_Nests_Hatched <- dplyr::filter(Ceuta_Nests, fate == "HATCH")

# subset the data to remove entries that have a "+" by the no_chicks
Ceuta_Nests <- Ceuta_Nests[!grepl("[+]",as.character(Ceuta_Nests$no_chicks)),]

# define no_chicks as.numeric
Ceuta_Nests$no_chicks <- as.numeric(as.character(Ceuta_Nests$no_chicks))

# insert zeros in situations where the fate of a nest was not unknown or hatched
Ceuta_Nests$no_chicks <- ifelse(Ceuta_Nests$fate == "UNKNOWN", NA,
                                ifelse(Ceuta_Nests$fate != "HATCH",
                                       0, Ceuta_Nests$no_chicks))

# make sure that there are no white spaces hidden in the ID info
Ceuta_Nests$year <- as.factor(gsub(" ", "", Ceuta_Nests$year, fixed = TRUE))
Ceuta_Nests$site <- as.factor(gsub(" ", "", Ceuta_Nests$site, fixed = TRUE))
Ceuta_Nests$nest <- as.factor(gsub(" ", "", Ceuta_Nests$nest, fixed = TRUE))
Ceuta_BirdRef$year <- as.factor(gsub(" ", "", Ceuta_BirdRef$year, fixed = TRUE))
Ceuta_BirdRef$site <- as.factor(gsub(" ", "", Ceuta_BirdRef$site, fixed = TRUE))
Ceuta_BirdRef$nest <- as.factor(gsub(" ", "", Ceuta_BirdRef$nest, fixed = TRUE))

# check the levels of the ID info to see if there are any obivous mistakes
levels(Ceuta_Nests$year)
levels(Ceuta_Nests$site)
levels(Ceuta_Nests$nest)
levels(Ceuta_BirdRef$year)
levels(Ceuta_BirdRef$site)
levels(Ceuta_BirdRef$nest)

# paste the year, site, and ID into one unique ID variable
Ceuta_Nests$ID <- paste(Ceuta_Nests$year, Ceuta_Nests$site, Ceuta_Nests$nest, sep = "_")
Ceuta_BirdRef$ID <- paste(Ceuta_BirdRef$year, Ceuta_BirdRef$site, Ceuta_BirdRef$nest, sep = "_")

# link the nest ID of the BirdRef file to the nest ID in the nest file.
Ceuta_Joined <- inner_join(Ceuta_Nests[,c("no_chicks", "clutch_size", "ID", "year")],
                           Ceuta_BirdRef[,c("male", "female", "ID")], by = "ID")
Ceuta_Joined$species <- "SP"
Ceuta_Joined_Sex <- Ceuta_Joined[, c(7,1:6)]
# same_sex <- Ceuta_Joined_Sex[which(Ceuta_Joined_Sex$sex_parent1 == Ceuta_Joined_Sex$sex_parent2),]
# same_sex <- same_sex[!duplicated(same_sex),]
#Ceuta_Joined_Sex <- anti_join(Ceuta_Joined_Sex, same_sex, by = "ID")
Ceuta_Joined_Sex <- Ceuta_Joined_Sex[!duplicated(Ceuta_Joined_Sex),]
Ceuta_Joined_Sex$pair <- paste(Ceuta_Joined_Sex$female, Ceuta_Joined_Sex$male, sep = "-")
Ceuta_Joined_Sex <- Ceuta_Joined_Sex[which(!is.na(Ceuta_Joined_Sex$female) & !is.na(Ceuta_Joined_Sex$male)),]
Ceuta_Joined_Sex$female <- as.factor(Ceuta_Joined_Sex$female)
Ceuta_Joined_Sex$male <- as.factor(Ceuta_Joined_Sex$male)
Ceuta_Joined_Sex$pair <- as.factor(Ceuta_Joined_Sex$pair)
df <- Ceuta_Joined_Sex[, c("species", "year", "ID", "female", "male", "pair")]
df <- df[!duplicated(df),]

females <- dcast(df, species + female  ~ year)
number_males_p_female <- aggregate(male ~ female, df, function(x) length(unique(x)))
females <- inner_join(females, number_males_p_female)
females[,c(3:9)] <- 
  lapply(females[,c(3:9)], as.numeric)
females$attempts <- rowSums(females[, c(3:9)])
females$years <- rowSums(females[, c(3:9)] > 0)
females_no_1 <- filter(females, male  != 1 | years != 1 | attempts != 1)
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(2,10)] <- c("focal", "mate")

males <- dcast(df, species + male  ~ year)
number_females_p_male <- aggregate(female ~ male, df, function(x) length(unique(x)))
males <- inner_join(males, number_females_p_male)
males[,c(3:9)] <- 
  lapply(males[,c(3:9)], as.numeric)
males$attempts <- rowSums(males[, c(3:9)])
males$years <- rowSums(males[, c(3:9)] > 0)
males_no_1 <- filter(males, female  != 1 | years != 1 | attempts != 1)
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(2,10)] <- c("focal", "mate")

mating_df_SP <- rbind(females_no_1, males_no_1)
mating_df_SP$status <- ifelse(mating_df_SP$mate == 1 & mating_df_SP$years == mating_df_SP$attempts, "Monogamous between years",
                              ifelse(mating_df_SP$mate == 1 & mating_df_SP$years < mating_df_SP$attempts, "Monogamous within years",
                                     ifelse(mating_df_SP$mate > 1 & mating_df_SP$years == mating_df_SP$attempts, "Polygamous between years",
                                            ifelse(mating_df_SP$mate > 1 & mating_df_SP$years < mating_df_SP$attempts, "Polygamous within years", "XXX"))))

mating_SP <- mating_df_SP[,c("species", "sex", "status", "focal")]


# Tuzla data import, tidying, and preprations
Tuzla_BirdRef <- 
  read.csv("data/Tuzla/BirdRef_Tuzla_1996_to_1999.csv",
           header=T, stringsAsFactors = FALSE)
Tuzla_Nests <- 
  read.csv("data/Tuzla/Tuzla_Nests.csv",
           header=T, stringsAsFactors = FALSE)

# Capitalize all column headers to keep scripts consistent between populations.
# First create a function that capitalizes the first letter of each word in a
# vector, then apply it to the column names of each dataset.
capwords <- function(s, strict = FALSE) 
{
  s <- tolower(s)
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

names(Tuzla_BirdRef) <- capwords(names(Tuzla_BirdRef))
names(Tuzla_Nests) <- capwords(names(Tuzla_Nests))

# subset the data to extract only nests that hatched
#Tuzla_Nests_Hatched <- dplyr::filter(Tuzla_Nests, fate == "HATCH")

# subset the data to remove entries that have a "+" by the no_chicks
Tuzla_Nests <- Tuzla_Nests[!grepl("[+]",as.character(Tuzla_Nests$No_chicks)),]

# define no_chicks as.numeric
Tuzla_Nests$No_chicks <- as.numeric(as.character(Tuzla_Nests$No_chicks))

# insert zeros in situations where the fate of a nest was not unknown or hatched
Tuzla_Nests$No_chicks <- ifelse(!grepl("[H/h]", Tuzla_Nests$Fate),
                                0, Tuzla_Nests$No_chicks)

# make sure that there are no white spaces hidden in the ID info
Tuzla_Nests$Year <- as.factor(gsub(" ", "", Tuzla_Nests$Year, fixed = TRUE))
Tuzla_Nests$Site <- as.factor(gsub(" ", "", Tuzla_Nests$Site, fixed = TRUE))
Tuzla_Nests$Nest <- as.factor(gsub(" ", "", Tuzla_Nests$Nest, fixed = TRUE))
Tuzla_BirdRef$Year <- as.factor(gsub(" ", "", Tuzla_BirdRef$Year, fixed = TRUE))
Tuzla_BirdRef$Site <- as.factor(gsub(" ", "", Tuzla_BirdRef$Site, fixed = TRUE))
Tuzla_BirdRef$Nest <- as.factor(gsub(" ", "", Tuzla_BirdRef$Nest, fixed = TRUE))

# check the levels of the ID info to see if there are any obivous mistakes
levels(Tuzla_Nests$Year)
Tuzla_Nests <- Tuzla_Nests[which(Tuzla_Nests$Year == "1996" |
                                   Tuzla_Nests$Year == "1997" |
                                   Tuzla_Nests$Year == "1998" |
                                   Tuzla_Nests$Year == "1999" |
                                   Tuzla_Nests$Year == "2000" |
                                   Tuzla_Nests$Year == "2001"),]
levels(Tuzla_Nests$Site)
levels(Tuzla_Nests$Nest)
levels(Tuzla_BirdRef$Year)
levels(Tuzla_BirdRef$Site)
levels(Tuzla_BirdRef$Nest)

# paste the year, site, and ID into one unique ID variable
Tuzla_Nests$ID <- paste(Tuzla_Nests$Year, Tuzla_Nests$Site, Tuzla_Nests$Nest, sep = "_")
Tuzla_BirdRef$ID <- paste(Tuzla_BirdRef$Year, Tuzla_BirdRef$Site, Tuzla_BirdRef$Nest, sep = "_")

# link the nest ID of the BirdRef file to the nest ID in the nest file.
Tuzla_joined <- inner_join(Tuzla_Nests[,c("No_chicks", "Clutch_size", "ID", "Year")],
                           Tuzla_BirdRef[,c("Male", "Female", "ID")], by = "ID")

Tuzla_joined$species <- "KP_Tuzla"
Tuzla_joined[, c("ID", "Male", "Female", "species")] <- 
  lapply(Tuzla_joined[, c("ID", "Male", "Female", "species")], as.factor)

Tuzla_joined_Sex <- Tuzla_joined [, c(7,1:6)]
# names(Tuzla_joined_Sex)[6] <- "parent1"
# names(Tuzla_joined_Sex)[7] <- "parent1"
# Tuzla_joined_Sex$sex_parent1 <- "M"
# Tuzla_joined_Sex$sex_parent2 <- "F"
Tuzla_joined_Sex <- Tuzla_joined[,c("species", "No_chicks", "Clutch_size",
                                    "ID", "Year", "Male", "Female")]
names(Tuzla_joined_Sex) <- c("species", "no_chicks", "clutch_size",
                             "ID", "year", "male", "female")

same_sex <- Tuzla_joined_Sex[which(Tuzla_joined_Sex$sex_parent1 == Tuzla_joined_Sex$sex_parent2),]
same_sex <- same_sex[!duplicated(same_sex),]
#write.csv(same_sex, "Data_files/Andavadoaka_families_with_homosexual_parents.csv",row.names = FALSE)

Tuzla_joined_Sex <- anti_join(Tuzla_joined_Sex, same_sex, by = "ID")
Tuzla_joined_Sex <- Tuzla_joined_Sex[!duplicated(Tuzla_joined_Sex),]
# Tuzla_joined_Sex$female <- ifelse(Tuzla_joined_Sex$sex_parent1 == "F", Tuzla_joined_Sex$parent1,
#                                   ifelse(Tuzla_joined_Sex$sex_parent2 == "F", Tuzla_joined_Sex$parent2, ""))
# Tuzla_joined_Sex$male <- ifelse(Tuzla_joined_Sex$sex_parent1 == "M", Tuzla_joined_Sex$parent1,
#                                 ifelse(Tuzla_joined_Sex$sex_parent2 == "M", Tuzla_joined_Sex$parent2, ""))
Tuzla_joined_Sex$pair <- paste(Tuzla_joined_Sex$female, Tuzla_joined_Sex$male, sep = "-")
Tuzla_joined_Sex <- Tuzla_joined_Sex[which(!is.na(Tuzla_joined_Sex$female) & !is.na(Tuzla_joined_Sex$male)),]
Tuzla_joined_Sex$female <- as.factor(Tuzla_joined_Sex$female)
Tuzla_joined_Sex$male <- as.factor(Tuzla_joined_Sex$male)
Tuzla_joined_Sex$pair <- as.factor(Tuzla_joined_Sex$pair)

df <- Tuzla_joined_Sex[, c("species", "year", "ID", "female", "male", "pair")]
df <- df[!duplicated(df),]
# MP FH72484 is monogamous 3 times within 2015 season
females <- dcast(df, species + female  ~ year)
number_males_p_female <- aggregate(male ~ female, df, function(x) length(unique(x)))
number_attempts_p_female <- aggregate(male ~ female, df, function(x) length(duplicated(x)))
test2 <- df[df$female %in% as.vector(number_attempts_p_female[which(number_attempts_p_female$male > 1),1]),]
test2 <- test2[with(test2,order(female,year)), ]

females <- inner_join(females, number_males_p_female)
females[,c(3:6)] <- 
  lapply(females[,c(3:6)], as.numeric)

females$attempts <- rowSums(females[, c(3:6)])
females$years <- rowSums(females[, c(3:6)] > 0)

#females$mate_fidelity <- females$male / females$years
females_no_1 <- filter(females, male  != 1 | years != 1 | attempts != 1)
# Kip FH47071 monogamous between and within seasons
# Kip FH47194 polygamous between seasons
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(2,7)] <- c("focal", "mate")


males <- dcast(df, species + male  ~ year)
number_females_p_male <- aggregate(female ~ male, df, function(x) length(unique(x)))
males <- inner_join(males, number_females_p_male)
males[,c(3:6)] <- 
  lapply(males[,c(3:6)], as.numeric)
males$attempts <- rowSums(males[, c(3:6)])
males$years <- rowSums(males[, c(3:6)] > 0)
#males$mate_fidelity <- males$female / males$years
males_no_1 <- filter(males, female  != 1 | years != 1 | attempts != 1)
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(2,7)] <- c("focal", "mate")

mating_df_KP <- rbind(females_no_1, males_no_1)
mating_df_KP$status <- ifelse(mating_df_KP$mate == 1 & mating_df_KP$years == mating_df_KP$attempts, "Monogamous between years",
                              ifelse(mating_df_KP$mate == 1 & mating_df_KP$years < mating_df_KP$attempts, "Monogamous within years",
                                     ifelse(mating_df_KP$mate > 1 & mating_df_KP$years == mating_df_KP$attempts, "Polygamous between years",
                                            ifelse(mating_df_KP$mate > 1 & mating_df_KP$years < mating_df_KP$attempts, "Polygamous within years", "XXX"))))

mating_KP_Tuzla <- mating_df_KP[,c("species", "sex", "status", "focal")]

###############################################################################
# MAIO data import, tidying, and preprations
Maio_BirdRef <- 
  read.csv("data/Maio/Maio_BirdRef_2007-2015.csv",
           header=T, stringsAsFactors = FALSE)
Maio_Nests <- 
  read.csv("data/Maio/Maio_Nests_2007-2015.csv",
           header=T, stringsAsFactors = FALSE)
# subset the data to extract only nests that hatched
#Maio_Nests_Hatched <- dplyr::filter(Maio_Nests, fate == "HATCH")

# Create a function that replaces all "#N/A" entries from excel version of the
# data with "NA" for R recognition.
NA_excel_replace <- function(dataframe)
{
  dataframe <- data.frame(lapply(dataframe, as.character), 
                          stringsAsFactors=FALSE)
  NAs <- dataframe == "#N/A"
  is.na(dataframe)[NAs] <- TRUE
  return(dataframe)
}

# apply the function to replace "#N/A" entires in all datasets
Maio_Nests <- NA_excel_replace(Maio_Nests)
Maio_BirdRef <- NA_excel_replace(Maio_BirdRef)

# Define the columns of the dataframes as the correct type.  First make a
# function, then apply it to the datasets
Define_plover_columns <- function(dataframe, type)
{
  if(type == "Captures")
  {
    #dataframe <- dataframe[,-1]
    dataframe[,c("year", "date", "time")] <- 
      lapply(dataframe[,c("year", "date", "time")], as.integer)
    dataframe[,c("site", "nest", "sex", "ring", "code", "observer")] <- 
      lapply(dataframe[,c("site", "nest", "sex", "ring", "code", "observer")],
             as.factor)
    dataframe[,c("comments_field", "comments_stdfile")] <- 
      lapply(dataframe[,c("comments_field", "comments_stdfile")], 
             as.character)
    dataframe[,c("weight", "left_wing", "right_wing", "left_tarsus", 
                 "right_tarsus", "bill")] <- 
      lapply(dataframe[,c("weight", "left_wing", "right_wing", "left_tarsus", 
                          "right_tarsus", "bill")], as.numeric)
    return(dataframe)
  }
  if(type == "Resightings")
  {
    dataframe[,c("year", "date", "time", "easting", "northing")] <- 
      lapply(dataframe[,c("year", "date", "time", "easting", "northing")], 
             as.integer)
    dataframe[,c("site", "sex", "code", "observer", "habitat")] <- 
      lapply(dataframe[,c("site", "sex", "code", "observer", "habitat")], 
             as.factor)
    dataframe[,c("comments")] <- lapply(dataframe[,c("comments")], 
                                        as.character)
    dataframe[,c("distance", "degree")] <- 
      lapply(dataframe[,c("distance", "degree")], as.numeric)
    return(dataframe)
  }
  if(type == "Birdref")
  {
    dataframe[,c("year")] <- as.integer(dataframe[,c("year")])
    dataframe[,c("site", "nest", "male", "female", "chick1", "chick2", "chick3",
                 "field_sex_m", "field_sex_f", "mol_sex_m", "mol_sex_f", 
                 "captured_focalyear_m", "captured_focalyear_f")] <- 
      lapply(dataframe[,c("site", "nest", "male", "female", "chick1", 
                          "chick2", "chick3", "field_sex_m", "field_sex_f", 
                          "mol_sex_m", "mol_sex_f", "captured_focalyear_m", 
                          "captured_focalyear_f")], as.factor)
    dataframe[,c("comments_field", "comments_stdfile")] <- 
      lapply(dataframe[,c("comments_field", "comments_stdfile")], as.character)
    return(dataframe)
  }
  if(type == "Broodfates")
  {
    dataframe[,c("year", "date", "time", "easting", "northing", "chicks")] <- 
      lapply(dataframe[,c("year", "date", "time", "easting", "northing", 
                          "chicks")], as.integer)
    dataframe[,c("site", "brood", "observer", "habitat", "parents")] <- 
      lapply(dataframe[,c("site", "brood", "observer", "habitat", "parents")], 
             as.factor)
    dataframe[,c("comments")] <- lapply(dataframe[,c("comments")], as.character)
    dataframe[,c("distance", "degree")] <- 
      lapply(dataframe[,c("distance", "degree")], as.numeric)
    return(dataframe)
  }
  if(type == "Nests")
  {
    #dataframe <- dataframe[,-1]
    dataframe[,c("year", "found_date", "floating_date", "found_time",
                 "laying_date", "end_date", "clutch_size", "no_chicks",
                 "easting", "northing")] <- 
      lapply(dataframe[,c("year", "found_date", "floating_date", "found_time", 
                          "laying_date", "end_date", "clutch_size", "no_chicks", 
                          "easting", "northing")], as.integer)
    dataframe[,c("site", "nest", "observer", "fate", "float1", "float2", 
                 "float3")] <- lapply(dataframe[,c("site", "nest", "observer", 
                                                   "fate", "float1", "float2", 
                                                   "float3")], as.factor)
    dataframe[,c("comments_field", "comments_stdfile")] <- 
      lapply(dataframe[,c("comments_field", "comments_stdfile")], as.character)
    dataframe[,c("length1", "width1", "length2", "width2", "length3", 
                 "width3")] <- lapply(dataframe[,c("length1", "width1", 
                                                   "length2", "width2", 
                                                   "length3", "width3")], 
                                      as.numeric)
    return(dataframe)
  }
  if(type == "Sexes")
  {
    dataframe[,c("mol_sex", "ring")] <- lapply(dataframe[,c("mol_sex", 
                                                            "ring")], as.character)
    return(dataframe)
  }
}

Maio_Nests <- Define_plover_columns(dataframe = Maio_Nests, type = "Nests")
Maio_BirdRef <- 
  Define_plover_columns(dataframe = Maio_BirdRef, type = "Birdref")

# subset the data to remove entries that have a "+" by the no_chicks
Maio_Nests <- Maio_Nests[!grepl("[+]",as.character(Maio_Nests$no_chicks)),]

# insert zeros in situations where the fate of a nest was not unknown or hatched
Maio_Nests$no_chicks <- ifelse(Maio_Nests$fate == "UNKNOWN", NA,
                               ifelse(Maio_Nests$fate != "HATCH",
                                      0, Maio_Nests$no_chicks))

# make sure that there are no white spaces hidden in the ID info
Maio_Nests$year <- as.factor(gsub(" ", "", Maio_Nests$year, fixed = TRUE))
Maio_Nests$site <- as.factor(gsub(" ", "", Maio_Nests$site, fixed = TRUE))
Maio_Nests$nest <- as.factor(gsub(" ", "", Maio_Nests$nest, fixed = TRUE))
Maio_BirdRef$year <- as.factor(gsub(" ", "", Maio_BirdRef$year, fixed = TRUE))
Maio_BirdRef$site <- as.factor(gsub(" ", "", Maio_BirdRef$site, fixed = TRUE))
Maio_BirdRef$nest <- as.factor(gsub(" ", "", Maio_BirdRef$nest, fixed = TRUE))

# check the levels of the ID info to see if there are any obivous mistakes
levels(Maio_Nests$year)
levels(Maio_Nests$site)
levels(Maio_Nests$nest)
levels(Maio_BirdRef$year)
levels(Maio_BirdRef$site)
levels(Maio_BirdRef$nest)

# paste the year, site, and ID into one unique ID variable
Maio_Nests$ID <- paste(Maio_Nests$year, Maio_Nests$site, Maio_Nests$nest, sep = "_")
Maio_BirdRef$ID <- paste(Maio_BirdRef$year, Maio_BirdRef$site, Maio_BirdRef$nest, sep = "_")

# link the nest ID of the BirdRef file to the nest ID in the nest file.
Maio_Joined <- inner_join(Maio_Nests[,c("no_chicks", "clutch_size", "ID", "year")],
                          Maio_BirdRef[,c("male", "female", "ID")], by = "ID")

Maio_Joined$species <- "KP_Maio"
Maio_Joined_Sex <- Maio_Joined[, c(7,1:6)]
same_sex <- Maio_Joined_Sex[which(Maio_Joined_Sex$sex_parent1 == Maio_Joined_Sex$sex_parent2),]
same_sex <- same_sex[!duplicated(same_sex),]
Maio_Joined_Sex <- anti_join(Maio_Joined_Sex, same_sex, by = "ID")
Maio_Joined_Sex <- Maio_Joined_Sex[!duplicated(Maio_Joined_Sex),]
Maio_Joined_Sex$pair <- paste(Maio_Joined_Sex$female, Maio_Joined_Sex$male, sep = "-")
Maio_Joined_Sex <- Maio_Joined_Sex[which(!is.na(Maio_Joined_Sex$female) & !is.na(Maio_Joined_Sex$male)),]
Maio_Joined_Sex$female <- as.factor(Maio_Joined_Sex$female)
Maio_Joined_Sex$male <- as.factor(Maio_Joined_Sex$male)
Maio_Joined_Sex$pair <- as.factor(Maio_Joined_Sex$pair)
df <- Maio_Joined_Sex[, c("species", "year", "ID", "female", "male", "pair")]
df <- df[!duplicated(df),]
females <- dcast(df, species + female  ~ year)
number_males_p_female <- aggregate(male ~ female, df, function(x) length(unique(x)))
number_attempts_p_female <- aggregate(male ~ female, df, function(x) length(duplicated(x)))
females <- inner_join(females, number_males_p_female)
females[,c(3:11)] <- 
  lapply(females[,c(3:11)], as.numeric)
females$attempts <- rowSums(females[, c(3:11)])
females$years <- rowSums(females[, c(3:11)] > 0)
females_no_1 <- filter(females, male  != 1 | years != 1 | attempts != 1)
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(2,12)] <- c("focal", "mate")

males <- dcast(df, species + male  ~ year)
number_females_p_male <- aggregate(female ~ male, df, function(x) length(unique(x)))
males <- inner_join(males, number_females_p_male)
males[,c(3:11)] <- 
  lapply(males[,c(3:11)], as.numeric)
males$attempts <- rowSums(males[, c(3:11)])
males$years <- rowSums(males[, c(3:11)] > 0)
#males$mate_fidelity <- males$female / males$years
males_no_1 <- filter(males, female  != 1 | years != 1 | attempts != 1)
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(2,12)] <- c("focal", "mate")

mating_df_KP_Maio <- rbind(females_no_1, males_no_1)
mating_df_KP_Maio$status <- ifelse(mating_df_KP_Maio$mate == 1 & mating_df_KP_Maio$years == mating_df_KP_Maio$attempts, "Monogamous between years",
                                   ifelse(mating_df_KP_Maio$mate == 1 & mating_df_KP_Maio$years < mating_df_KP_Maio$attempts, "Monogamous within years",
                                          ifelse(mating_df_KP_Maio$mate > 1 & mating_df_KP_Maio$years == mating_df_KP_Maio$attempts, "Polygamous between years",
                                                 ifelse(mating_df_KP_Maio$mate > 1 & mating_df_KP_Maio$years < mating_df_KP_Maio$attempts, "Polygamous within years", "XXX"))))

mating_KP_Maio <- mating_df_KP_Maio[,c("species", "sex", "status", "focal")]

###############################################################################
# COMPARATIVE plotting and analysis
total_mating_df <- rbind(mating_Mada, mating_SP, mating_KP_Maio, mating_KP_Tuzla)

total_mating_df$species <- ifelse(total_mating_df$species == "SP", "Snowy", 
                                  ifelse(total_mating_df$species == "KP_Tuzla", "Kentish-Tuzla",
                                         ifelse(total_mating_df$species == "MP", "Madagascar",
                                                ifelse(total_mating_df$species == "WFP", "White-fronted",
                                                       ifelse(total_mating_df$species == "KP_Maio", "Kentish-Maio", "Kittlitz's")))))

total_mating_df$status <- factor(total_mating_df$status, 
                                 levels = c("Monogamous within years",
                                            "Monogamous between years",
                                            "Polygamous between years",
                                            "Polygamous within years"))
total_mating_df$species <- factor(total_mating_df$species, 
                                  levels = c("Snowy",
                                             "Kentish-Tuzla",
                                             "Madagascar",
                                             "Kentish-Maio",
                                             "White-fronted",
                                             "Kittlitz's"))
total_mating_df$status_simple <- ifelse(total_mating_df$status == "Monogamous within years" |
                                          total_mating_df$status == "Monogamous between years",
                                 "Monogamous", "Polygamous")


sample_sizes_sex <- aggregate(focal ~ species + sex, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes <- aggregate(focal ~ species, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes_status_simple <- aggregate(focal ~ species + status_simple, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes_status_simple <- left_join(sample_sizes_status_simple, sample_sizes, by = c("species"))
sample_sizes_status_simple$prop <- sample_sizes_status_simple$focal.x/sample_sizes_status_simple$focal.y

sample_sizes_sex <- aggregate(focal ~ species + sex, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes <- aggregate(focal ~ species, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes_status<- aggregate(focal ~ species + status, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes_status <- left_join(sample_sizes_status, sample_sizes, by = c("species"))
sample_sizes_status$prop <- sample_sizes_status$focal.x/sample_sizes_status$focal.y
colnames(sample_sizes_status) <- c("Species", "Mating_status", "N_yes", "N_no", "Proportion")
# write.table(sample_sizes_status, 
#             file = "/Users/Luke/Documents/Academic_Projects/PhD/Charadrius_Mating_System_Chapter/Book_chapter_R_analysis/Data_files/Social_mating_system_summary_LEP.txt",
#             sep = "\t",
#             row.names = FALSE)

sample_sizes_status_simple_sex <- aggregate(focal ~ species + status_simple + sex, data = total_mating_df, FUN = function(x){NROW(x)})
sample_sizes_status_simple_sex <- left_join(sample_sizes_status_simple_sex, sample_sizes_sex, by = c("species", "sex"))
sample_sizes_status_simple_sex$prop <- sample_sizes_status_simple_sex$focal.x/sample_sizes_status_simple_sex$focal.y

polygamy_prop <- filter(sample_sizes_status_simple, status_simple == "Polygamous")[,c(1,5)]
colnames(polygamy_prop) <- c("population", "prop_poly")
polygamy_prop$population <- 
  factor(polygamy_prop$population ,
         levels = c("Snowy",
                    "Kentish (Tuzla)",
                    "Madagascar",
                    "Kentish (Maio)",
                    "White-fronted",
                    "Kittlitz's"))

polygamy_prop_sex <- filter(sample_sizes_status_simple_sex, status_simple == "Polygamous")[,c(1,3,6)]
colnames(polygamy_prop_sex) <- c("population", "sex", "prop_poly")
polygamy_prop_sex$population <- 
  factor(polygamy_prop_sex$population ,
         levels = c("Snowy",
                    "Kentish (Tuzla)",
                    "Madagascar",
                    "Kentish (Maio)",
                    "White-fronted",
                    "Kittlitz's"))


custom_pal <- c("#7b3294", "#9E6BB1", "#53B16B", "#008837")

matefidelity_plot_by_sex <- 
  ggplot() +
  geom_bar(position = "fill", alpha = 0.75, data = total_mating_df, aes(x = species, fill = status)) +
  geom_text(aes(y = 1.05, x = species, label = focal, family = "Arial"), data = sample_sizes_sex, size = 6) +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        legend.text=element_text(size = 11),
        legend.title=element_blank(),
        # legend.key.height=unit(0.8,"line"),
        # legend.key.width=unit(0.8,"line"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 13), 
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 13), 
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size=14, face = "bold"),
        strip.text.y = element_text(size=14, face = "bold")) +
  ylab("Proportion of individuals") +
  scale_fill_manual(values = custom_pal) +
  facet_grid(sex ~ .) +
  scale_y_continuous(limits = c(0, 1.05))

custom_pal <- 
  c(brewer.pal(8, "Set1")[c(2)], brewer.pal(8, "Set1")[c(1)])
matefidelity_plot <- 
  ggplot() +
  #coord_flip() +
  geom_bar(position = "fill", alpha = 0.75, data = total_mating_df, aes(x = species, fill = status_simple)) +
  # geom_text(aes(y = 1.05, x = species, label = c("Snowy", "Kentish (Tuzla)", "Madagascar", "Kentish (Maio)", "White-fronted", "Kittlitz's"), 
  #               family = "Arial"), data = sample_sizes, size = 3) +
  #geom_text(aes(y = 1.05, x = species, label = paste("n = ", focal, sep = ""), family = "Candara"), data = sample_sizes, size = 4) +
  theme_bw() +
  theme(text=element_text(family="Menlo"),
        legend.text=element_text(size = 11),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        axis.title.y = element_text(size=12, vjust=-0.2),
        axis.text.y  = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0.8,1.2,0.9,1), "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size=11),
        panel.margin = unit(0.75, "lines")) +
  # theme(text = element_text(family = "Arial"),
  #       legend.text=element_text(size = 11),
  #       legend.title=element_blank(),
  #       legend.position="top",
  #       legend.key.height=unit(0.8,"line"),
  #       legend.key.width=unit(0.8,"line"),
  #       axis.title.y = element_blank(),
  #       axis.text.y  = element_text(size = 11), 
  #       axis.title.x = element_text(size = 12),
  #       axis.text.x  = element_blank(), 
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       # strip.text.x = element_text(size=14, face = "bold"),
  #       # strip.text.y = element_text(size=14, face = "bold"),
  #       plot.margin = unit(c(0,1,1,1), "cm"),
  #       panel.border=element_blank(),
  #       axis.ticks.y = element_blank(),
  #       axis.ticks.x = element_blank()) +
  ylab("Proportion of individuals") +
  scale_fill_manual(values = custom_pal) +
  # scale_y_continuous(limits = c(0, 1.05), breaks = NULL) +
  # scale_x_discrete(breaks = NULL) +
  guides(fill = guide_legend(ncol = 2, bycol = TRUE))
matefidelity_plot

ggplot2::ggsave(matefidelity_plot, 
                filename = "matefidelity_plot.jpg", 
                path = "figs/",
                width = 10,
                height = 4.5, units = "in",
                dpi = 300)

# extract the female column, add a sex column.  extract the male colum, add a 
# sex column.  Stack these two dataframes.
Sex <- rep("Female", nrow(Andava_Joined_Sex))
Ring <- Andava_Joined_Sex$female
females <- data.frame(Ring, Sex)
Sex <- rep("Male", nrow(Andava_Joined_Sex))
Ring <- Andava_Joined_Sex$male
males <- data.frame(Ring, Sex)
Individuals <- rbind(males, females)

# Replicate each row by 2 then cbind the stacked dataframe from the previous
# step
Andava_Joined_Sex <- cbind(Andava_Joined_Sex[rep(row.names(Andava_Joined_Sex), 2), c("no_chicks", "clutch_size", "ID", "year", "species")],
                      Individuals)

# Change the order of the sex levels, so that females are first (for the plot)
Andava_Joined_Sex$Sex <- factor(Andava_Joined_Sex$Sex, levels = c("Female", "Male"))

# subset the data to remove entries that have a NA in the Ring column
Andava_Joined_Sex <- Andava_Joined_Sex[!is.na(Andava_Joined_Sex$Ring),]

# subset the data to remove entries that have a NA in the Ring column
Andava_Joined_Sex <- Andava_Joined_Sex[!is.na(Andava_Joined_Sex$no_chicks),]

# group data according to Year, Sex, then Ring
Andava_Joined_Sex <- group_by(Andava_Joined_Sex, species, year, Sex, Ring)

# sum the total chicks produced per bird each year
Andava_Joined_Sex_sum <- ungroup(dplyr::summarise(Andava_Joined_Sex, total_chicks_p_year = sum(as.numeric(no_chicks))))

# group data according to Sex then Ring
Andava_Joined_Sex_sum <- group_by(Andava_Joined_Sex_sum, species, Sex, Ring)

# calculate avg total chicks produced per bird each year
Andava_Joined_Sex_sum_avg <- ungroup(dplyr::summarise(Andava_Joined_Sex_sum, avg_chicks_p_year = mean(as.numeric(total_chicks_p_year))))

# summarize the avg annual no_chicks by sex
fecundity_summary <- Rmisc::summarySE(Andava_Joined_Sex_sum_avg, measurevar = "avg_chicks_p_year", groupvars = c("species", "Sex"))

# subset data by species
KiP_Andava_Joined_Sex_sum_avg <- Andava_Joined_Sex_sum_avg[which(Andava_Joined_Sex_sum_avg$species == "KIP"),]
MP_Andava_Joined_Sex_sum_avg <- Andava_Joined_Sex_sum_avg[which(Andava_Joined_Sex_sum_avg$species == "MP"),]
WfP_Andava_Joined_Sex_sum_avg <- Andava_Joined_Sex_sum_avg[which(Andava_Joined_Sex_sum_avg$species == "WFP"),]

# plot distributions of fecundity by sex
cbPalette <- c("#CC0000", "#0066CC")
KiP_Sex_fecund_plot <- ggplot(aes(y = avg_chicks_p_year, x = Sex, fill = Sex), data = KiP_Andava_Joined_Sex_sum_avg) + 
  theme_bw() +
  geom_boxplot() +
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = "none",
        #legend.text=element_text(size=30),
        #legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=30, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        #plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  #panel.background = element_rect(colour = "black", size=1)) +
  #scale_colour_brewer(palette = "Set1") + 
  scale_fill_manual(values=cbPalette) +
  xlab("") + 
  ylab("Avg. annual fecundity per individual") #+
#scale_y_continuous(limits = c(0, 1)) +
#scale_fill_hue(l=55)
KiP_Sex_fecund_plot

cbPalette <- c("#CC0000", "#0066CC")
MP_Sex_fecund_plot <- ggplot(aes(y = avg_chicks_p_year, x = Sex, fill = Sex), data = MP_Andava_Joined_Sex_sum_avg) + 
  theme_bw() +
  geom_boxplot() +
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = "none",
        #legend.text=element_text(size=30),
        #legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=30, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        #plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  #panel.background = element_rect(colour = "black", size=1)) +
  #scale_colour_brewer(palette = "Set1") + 
  scale_fill_manual(values=cbPalette) +
  xlab("") + 
  ylab("Avg. annual fecundity per individual") #+
#scale_y_continuous(limits = c(0, 1)) +
#scale_fill_hue(l=55)
MP_Sex_fecund_plot

cbPalette <- c("#CC0000", "#0066CC")
WfP_Sex_fecund_plot <- ggplot(aes(y = avg_chicks_p_year, x = Sex, fill = Sex), data = WfP_Andava_Joined_Sex_sum_avg) + 
  theme_bw() +
  geom_boxplot() +
  theme(text=element_text(size=16, family="Candara"), # set the font as Candara
        legend.position = "none",
        #legend.text=element_text(size=30),
        #legend.title=element_text(size=30),
        axis.title.x = element_text(size=35, vjust=-0.1),
        axis.text.x  = element_text(size=30), 
        axis.title.y = element_text(size=30, vjust=1.2),
        axis.text.y  = element_text(size=30), 
        #plot.title = element_text(vjust=1.2, size=20, face="bold"),
        panel.grid.major = element_blank()) +
  #panel.background = element_rect(colour = "black", size=1)) +
  #scale_colour_brewer(palette = "Set1") + 
  scale_fill_manual(values=cbPalette) +
  xlab("") + 
  ylab("Avg. annual fecundity per individual") #+
#scale_y_continuous(limits = c(0, 1)) +
#scale_fill_hue(l=55)
WfP_Sex_fecund_plot
