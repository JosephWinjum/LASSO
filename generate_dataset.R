library(readxl)
library(data.table)
library(pwr)
source("../code/functions.R")
data <- data.table(read_excel("ssrp_data.xlsx", na="N/A"))

################################################################################
# Match author data using same procedure as in merge_data.R
authors <- data.table(read_excel("authors.xlsx", na="N/A"))
setkey(authors, study_id, study_type, author_order)

# `author_data.csv` includes info information on all authors (citations etc)
author_data <- data.table(read.csv("../data/author_data.csv", sep = ";",
                                   fileEncoding = "macroman"))
setkey(author_data, author_id) # one row per author

# Add variables for number of authors
n_authors.o <- unique(authors[study_type=="o", .(study_id, n_authors)])
names(n_authors.o) <- c("id", "n_authors.o")
n_authors.r <- unique(authors[study_type=="r", .(study_id, n_authors)])
names(n_authors.r) <- c("id", "n_authors.r")
data <- merge(data, merge(n_authors.o, n_authors.r, by="id"), by="id")

# For each study, find the relevant authors
for (id in data$id) {

    # loop over replication and original studies
    for (type in c("o", "r")) {

        study_authors <- authors[study_id == id & study_type == type, author_id]
        if (length(study_authors) > 0 & study_authors[1] != "") {
            # Instead of saving individual author info in the data file do the aggregation directly

            # Average citations
            data[id, paste0("author_citations_avg.", type) :=
                 mean(author_data[study_authors, citations], na.rm=TRUE)]

            # Citations of most cited author
            data[id, paste0("author_citations_max.", type) :=
                 max(author_data[study_authors, citations])]

            # Gender ratio (percent male)
            data[id, paste0("authors_male.", type) :=
                 mean(author_data[study_authors, gender] == "M", na.rm=TRUE)]

            # Highest seniority (highest seniority on project)
            seniority_order <- c("Professor", "Associate Professor",
                                 "Assistant Professor", "Researcher",
                                 "Lecturer", "Assistant", "Other")
            data[id, paste0("highest_seniority.", type) :=
                 seniority_order[seniority_order %in%
                                 author_data[study_authors, position]][1]]
        }
    }
}

################################################################################
# Generate some variables

# Same country and language
data[, same_country := ifelse(experiment_country.o == experiment_country.r, 1, 0)]
data[, same_language := ifelse(experiment_language.o == experiment_language.r, 1, 0)]
data[, same_online := ifelse(online.o == online.r, 1, 0)]
data[, same_subjects := ifelse(subjects.o == subjects.r, 1, 0)]

# US original/replication
data[, us_lab.o := ifelse(experiment_country.o == "United States", 1, 0)]
data[, us_lab.r := ifelse(experiment_country.r == "United States", 1, 0)]

# Original power
data$power.o <- vector("numeric", length = nrow(data))
for (i in 1:nrow(data)) {
    set(data, i = i, j = "power.o",
        value = pwr.r.test(n = data[i, n_groups.o],
                           r = data[i, effect_size.o],
                           power = NULL, sig.level = 0.05)$power)

    set(data, i = i, j = "es_80power1",
        value = pwr.r.test(n = newdata$n_groups_planned1.r[i],
                           r = NULL, sig.level = 0.05, power = 0.8)$r)

    set(data, i = i, j = "es_80power2",
        value = pwr.r.test(n = newdata$n_groups_planned2.r[i],
                           r = NULL, sig.level = 0.05, power = 0.8)$r)
}

################################################################################
# We need to create factors of all categorical variables and make sure that the
# levels are the same as in the original data set for the predictions to work
olddata <- readRDS("../data/data.rds")
# remove parts of data set that we don't use so that factor levels are correct
olddata <- olddata[aggregated == TRUE | is.na(aggregated)]
olddata <- olddata[drop == FALSE]

# Since we don't have any replication results yet we fill with random
if (all(is.na(data$replicated))) data$replicated <- rep_len(c(0,1), nrow(data))
if (all(is.na(data$relative_es))) data$relative_es <- rep_len(c(0,1), nrow(data))

# ssrp.12 is really "Education" but we dont have any training data on this
# discipline, so "Cognitive" is closest
data$discipline[data$discipline == "Education"] <- "Cognitive"

# Creating factors with the correct levels (taken from original data)
data$replicated <- factor(data$replicated)
levels(data$replicated) <- c("not_replicated", "replicated")
data$discipline <- factor(data$discipline, levels = levels(olddata$discipline))
data$effect_type.o <- factor(data$effect_type.o, levels = levels(as.factor(olddata$effect_type.o)))
data$compensation.o <- factor(data$compensation.o, levels = levels(as.factor(olddata$compensation.o)))
data$compensation.r <- factor(data$compensation.r, levels = levels(as.factor(olddata$compensation.r)))

lvls <- c("Assistant", "Assistant Professor", "Associate Professor", "Professor", "Researcher")
data$highest_seniority.r <- factor(data$highest_seniority.r, levels = lvls)
data$highest_seniority.o <- factor(data$highest_seniority.o, levels = lvls)

# These variables are categorical but not used in model, so we don't save them as factors
# For example subjects.o has the incorrect levels so will create NA's
# data$subjects.o <- factor(data$subjects.o, levels = levels(as.factor(olddata$subjects.o)))
# data$subjects.r <- factor(data$subjects.r, levels = levels(as.factor(olddata$subjects.r)))
# data$experiment_country.o <- factor(data$experiment_country.o, levels = levels(as.factor(olddata$experiment_country.o)))
# data$experiment_country.r <- factor(data$experiment_country.r, levels = levels(as.factor(olddata$experiment_country.r)))
# data$experiment_language.o <- factor(data$experiment_language.o, levels = levels(as.factor(olddata$experiment_language.o)))
# data$experiment_language.r <- factor(data$experiment_language.r, levels = levels(as.factor(olddata$experiment_language.r)))

################################################################################
# Save data as csv and RDS

write.csv(data, "newdata.csv", row.names = FALSE, na = "",
          fileEncoding = "UTF-8")
saveRDS(data, "newdata.rds")
