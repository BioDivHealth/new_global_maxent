## ----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\#
## Extract the taxonomic information and list of sinonyms for a given specie
## ----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~#
#


clean_synonyms <- function(x) {
                x %>%
                    as.character() %>% # be sure everything is character
                    gsub('c\\(|\\)|"|\\\\', "", .) %>% # strip c(  )  "  and escaped quotes
                    strsplit(split = ",") %>% # break on commas
                    unlist() %>%
                    trimws() %>% # remove leading / trailing spaces
                    discard(~ .x %in% c("", "NA")) %>% # drop empties / literal "NA"
                    unique() %>% # deduplicate
                    paste(collapse = ", ")
            }

clean_synonyms2 <- function(x) {
    x %>%
        as.character() %>% # be sure everything is character
        gsub('c\\(|\\)|"|\\\\', "", .) %>% # strip c(  )  "  and escaped quotes
        strsplit(split = "[,;]") %>% # break on commas
        unlist() %>%
        trimws() %>% # remove leading / trailing spaces
        discard(~ .x %in% c("", "NA")) %>% # drop empties / literal "NA"
        unique() %>% # deduplicate
        paste(collapse = "; ") %>% 
        strsplit(split = "[,;]") %>% 
        unlist() %>%
        trimws() %>% # remove leading / trailing spaces
        discard(~ .x %in% c("", "NA")) %>% # drop empties / literal "NA"
        unique() %>% 
        paste(collapse = "; ")
    # deduplicate
}


clean__genus_synonyms <- function(x) {
                x %>%
                    as.character() %>% # be sure everything is character
                    strsplit(split = ";") %>% # break on semicolons
                    unlist() %>%
                    trimws() %>% # remove leading / trailing spaces
                    discard(~ .x %in% c("", "NA")) %>% # drop empties / literal "NA"
                    unique() %>% # deduplicate
                    paste(collapse = "; ")
            }


get_synonyms <- function(genus, species, infra = NULL,
                         subpopulation = NULL, key = NULL, ...) {
    rsp <- tryCatch(
        rredlist::rl_species(
            genus = genus, species = species,
            infra = infra, subpopulation = subpopulation,
            key = key, ...
        ),
        error = function(e) NULL # network/API failure
    )
    if (is.null(rsp)) {
        return(NULL)
    }
    syns <- rsp$taxon$synonyms
    if (is.null(syns) || !is.data.frame(syns) || nrow(syns) == 0) {
        return(NULL)
    }

    syns
}

safe_rl_species_latest <- function(genus, species,
                                   infra = NULL, subpopulation = NULL,
                                   key = NULL, parse = TRUE) {
    assess <- rredlist::rl_species(genus, species,
        infra = infra
    )$assessments

    assess$year_published <- as.numeric(assess$year_published)

    if (any(assess$latest, na.rm = TRUE)) {
        assess <- subset(assess, assess$latest)
    }

    # Fallback when nothing is flagged
    if (nrow(as.data.frame(assess)) == 0) {
        warning("No assessment marked 'latest'; using most recent year")
        assess <- assess[order(assess$year_published, decreasing = TRUE), , drop = FALSE]
    }

    rredlist::rl_assessment(
        id = assess$assessment_id[1],
        key = key,
        parse = parse, ...
    )
}


rl_species_latest <- function(genus, species, infra = NULL,
                              subpopulation = NULL,
                              key = NULL, parse = TRUE) {
    tmp <- rl_species(genus, species,
        infra = infra,
        subpopulation = subpopulation
    )$assessments
    if (any(tmp$latest, na.rm = TRUE)) {
        tmp <- subset(tmp, tmp$latest)
    }
    if (length(tmp) == 0) {
        stop("No assessments found for this species.")
    }
    tmp$year_published <- as.numeric(as.character(tmp$year_published))
    ord <- order(tmp$year_published, decreasing = TRUE)
    tmp <- tmp[ord, , drop = FALSE]
    rl_assessment(id = tmp$assessment_id[1], key = key, parse = parse)
}


# IUCN_get_data function
IUCN_get_data <- function(g, b, genus, species, infra = NULL, n_times = 3) {
    if (length(g) == 0 & length(b) == 0) {
        IUCN_id <- NA
        IUCN_name <- NA
        IUCN_Phylum <- NA
        IUCN_Class <- NA
        IUCN_Order <- NA
        IUCN_Family <- NA

        IUCN_Category <- NA
        IUCN_Present <- "No"
        IUCN_latest <- NA
        IUCN_date <- NA
        IUCN_Category <- NA

        IUCN_N_syn <- NA
        IUCN_status <- NA
        IUCN_syn <- NA
    } else {
        # if the name of the species correspond to a synonim in the IUCN-red list,
        # use the accepted name to retrieve the species information

        if (length(b) != 0 & length(g) == 0) {
            g <- NULL
            t_11 <- 1

            while (is.null(g) && t_11 <= n_times) {
                acc <- b %>% filter(status == "ACCEPTED")
                # replace NA in acc$infra_name with NULL
                infra <- NULL
                if (!is.na(acc$infra_name[1])) {
                    infra <- acc$infra_name[1]
                }
                try(g <- rl_species(genus = acc$genus_name[1], species = acc$species_name[1], infra = NULL)) # get the taxonomic information for the accepted name)
                # ry(g <- rl_search(name = b$result$accepted_name[1],silent=TRUE)$result)
                t_11 <- t_11 + 1
            }
        }

        if (length(g) != 0) {
            IUCN_id <- unique(g$taxon$sis_id)
            IUCN_name <- unique(g$taxon$scientific_name)

            if (is.null(g$taxon$phylum_name)) {
                IUCN_Phylum <- NA
            } else {
                IUCN_Phylum <- unique(g$taxon$phylum_name)
            }

            if (is.null(g$taxon$class_name)) {
                IUCN_Class <- NA
            } else {
                IUCN_Class <- unique(g$taxon$class_name)
            }

            if (is.null(g$taxon$order_name)) {
                IUCN_Order <- NA
            } else {
                IUCN_Order <- unique(g$taxon$order_name)
            }

            if (is.null(g$taxon$family_name)) {
                IUCN_Family <- NA
            } else {
                IUCN_Family <- unique(g$taxon$family_name)
            }

            if (!is.null(g$assessments$latest)) {
                latest <- g$assessments %>% filter(latest == TRUE)
                IUCN_Category <- latest$red_list_category_code
                IUCN_latest <- latest$latest
                IUCN_date <- latest$year_published
            }
            IUCN_Present <- "Yes"
        }

        if (length(b) != 0) {
            syns <- b %>%
                # filter(status == "ACCEPTED" | status=="ADD") %>%
                mutate(
                    full_synonym_name = ifelse(
                        is.na(infra_name),
                        paste(genus_name, tolower(species_name)),
                        paste(genus_name, tolower(species_name), tolower(infra_name))
                    )
                )
            z_c <- syns$full_synonym_name %>%
                strsplit(split = " ") %>%
                sapply(length) == 2 # check if the synonym is a species or sub-species

            IUCN_N_syn <- length(z_c[z_c == TRUE]) # number of sp sinonims, Subsp excluded
            IUCN_status <- paste(syns$status, collapse = ";")
            IUCN_syn <- (paste(syns$full_synonym_name[sapply(strsplit(syns$full_synonym_name, split = " "), length) == 2],
                collapse = ";"
            )) # combine the names into a single string with all the sinonims
        } else {
            IUCN_status <- NA
            IUCN_N_syn <- NA
            IUCN_syn <- NA
        }
    }

    # Combine the IUCN information
    cat("IUCN Present:", IUCN_Present, "\n")
    cat("IUCN ID:", IUCN_id, "\n")
    cat("IUCN Name:", IUCN_name, "\n")
    cat("IUCN Category:", IUCN_Category, "\n")
    cat("IUCN Latest:", IUCN_latest, "\n")
    cat("IUCN Date:", IUCN_date, "\n")
    cat("IUCN Synonyms:", IUCN_syn, "\n")
    cat("IUCN Number of Synonyms:", IUCN_N_syn, "\n")
    cat("IUCN Status:", IUCN_status, "\n")
    cat("IUCN Phylum:", IUCN_Phylum, "\n")
    cat("IUCN Class:", IUCN_Class, "\n")
    cat("IUCN Order:", IUCN_Order, "\n")
    cat("IUCN Family:", IUCN_Family, "\n")
    # ------------------------------------------------------------------
    # Replace zero-length objects with NA ------------------------------
    # ------------------------------------------------------------------
    vars <- c("IUCN_Present", "IUCN_id", "IUCN_name", "IUCN_latest",
              "IUCN_date", "IUCN_Category", "IUCN_N_syn", "IUCN_syn",
              "IUCN_status", "IUCN_Phylum", "IUCN_Class",
              "IUCN_Order", "IUCN_Family")
    
    for (v in vars) {
        if (exists(v, inherits = TRUE) && length(get(v)) == 0) {
            assign(v, NA, inherits = TRUE)
        }
    }
    IUCN_data <- unique(data.frame(
        IUCN_Present, IUCN_id, IUCN_name, IUCN_latest, IUCN_date,
        IUCN_Category, IUCN_N_syn, IUCN_syn, IUCN_status,
        IUCN_Phylum, IUCN_Class,
        IUCN_Order, IUCN_Family
    ))

    return(IUCN_data)
}

# ITIS_get_data function
ITIS_get_species_data <- function(spp.x, n_times = 3) {
    # b.3.1 Get TSN a reference number that we are going to need to gather information from ITIS----
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    TSN <- NULL
    ITIS_species_in_genus <- NA
    t_4 <- 1

    while (is.null(TSN) && t_4 <= n_times) {
        try(TSN <- get_tsn_(spp.x, searchtype = "scientific", silent = TRUE))
        t_4 <- t_4 + 1
        print(TSN)
    }

    if (is.list(TSN) & is.null(TSN[[1]])) {
        TSN[[1]] <- data.frame()
    }
    if(is.null(TSN)){
        TSN = list(data.frame(c(NULL)))
        #TSN[[1]] = data.frame(list(c(NA)), list(c(NA)))
    }
    if (nrow(TSN[[1]]) == 0) {
        try(TSN.two <- get_tsn(spp.x, searchtype = "scientific", silent = TRUE))
        if(!exists("TSN.two")){
            TSN.two = NA
            attr(TSN.two) = "does not exists"
        }
        if (attr(TSN.two, "match") == "found") {
            tsn2 <- itis_acceptname(TSN.two[[1]], silent = TRUE)
            ## 1. Add a new row (tibble is currently empty)
            TSN[[1]] <- TSN[[1]] |>
                add_row(
                    tsn            = TSN.two[[1]],
                    scientificName = spp.x,
                    commonNames    = NA,
                    nameUsage      = "valid"
                )
        } else if (attr(TSN.two, "match") == "not found") {
            TSN <- list(c(NULL))
        } else if(attr(TSN.two, "match") == "does not exist"){
            TSN <- list(c(NULL))
        }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    if (is.null(TSN[[1]])) {
        ITIS_Present <- "No"
        ITIS_id <- NA
        ITIS_name <- NA
        ITIS_is_valid <- NA
        ITIS_Phylum <- NA
        ITIS_Class <- NA
        ITIS_Order <- NA
        ITIS_Family <- NA

        ITIS_accept <- NA
        ITIS_syn <- NA
        ITIS_N_syn <- NA
    } else {
        tsn <- TSN[[1]]

        tsn_n <- tsn[tsn$scientificName == spp.x & tsn$nameUsage == "valid", ]

        if (length(tsn_n$tsn) == 0) { # if the original species name is not present in the returned list of synonims we igore the results

            ITIS_Present <- "Wrong"
            ITIS_id <- NA
            ITIS_name <- NA
            ITIS_is_valid <- NA
            ITIS_Phylum <- NA
            ITIS_Class <- NA
            ITIS_Order <- NA
            ITIS_Family <- NA

            ITIS_accept <- NA
            ITIS_syn <- NA
            ITIS_N_syn <- NA
        } else {
            ITIS_Present <- "Yes"

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            AceptName <- NULL
            t_4 <- 1

            while (is.null(AceptName) && t_4 <= n_times) {
                try(AceptName <- itis_acceptname(tsn_n$tsn, silent = TRUE))
                t_4 <- t_4 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

            if (AceptName$submittedtsn == AceptName$acceptedtsn) { # is the submitted Or_name name accepted by ITIS
                ITIS_name <- spp.x
                ITIS_accept <- "yes"
            } else if (!is.na(AceptName$acceptedname)) {
                ITIS_name <- AceptName$acceptedname # if the name is not accepted, we use the accepted name
                ITIS_accept <- "yes"
            } else if (is.na(AceptName$acceptedname)) {
                ITIS_name <- "No"
                ITIS_accept <- "No"
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            Syn <- NULL
            t_5 <- 1

            while (is.null(Syn) && t_5 <= n_times) {
                try(Syn <- synonyms(tsn$tsn, db = "itis"))
                t_5 <- t_5 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

            # Syn_l is a list of possible codes with their respective names
            # We need to unify the list and extract the names

            Syn <- rbindlist(Syn, fill = TRUE)
            # count the number of FALSE and TRUE       #strsplit(, " "), length
            z <- as.character(is.na(Syn$syn_name))
            # if the length equals cero, there is no syn names

            if (length(z[z == "FALSE"]) == 0) {
                ITIS_N_syn <- NA
                ITIS_syn <- NA
            } else {
                z_b <- as.character(sapply(strsplit(Syn$syn_name, split = " "), length) == 2)

                ITIS_N_syn <- length(z_b[z_b == TRUE]) # number of sp sinonims, Subsp excluded
                ITIS_syn <- paste(Syn$syn_name[sapply(strsplit(Syn$syn_name, split = " "), length) == 2],
                    collapse = ";"
                ) # combine the names into a single string with all the sinonims
            }

            ITIS_id <- tsn_n$tsn
            if (is.null(ITIS_name) | is.na(ITIS_name)) {
                ITIS_name <- tsn_n$scientificName
            }
            ITIS_is_valid <- tsn_n$nameUsage

            # Get the upstream taxonomic information
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            j <- NULL
            t_13 <- 1
            t_31 <- 1
            t_32 <- 1
            while (is.null(j) && t_31 <= n_times) {
                try(j <- itis_hierarchy(tsn_n$tsn, what = "full", silent = TRUE))
                t_31 <- t_31 + 1
                if (dim(j)[1]>0){
                j = j %>% filter(tsn!=tsn_n$tsn)}
            }
            while (dim(j)[1] == 0 && t_32 <= n_times) {
                try(j <- itis_hierarchy(AceptName$acceptedtsn, what = "full", silent = TRUE))
                t_32 <- t_32 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            if (dim(j)[1] > 0) {
                ITIS_Phylum <- toupper(as.character(j[j$rankname == "phylum", 4])) # Change from lower to upper chase
                ITIS_Class <- toupper(as.character(j[j$rankname == "class", 4]))
                ITIS_Order <- toupper(as.character(j[j$rankname == "order", 4]))
                ITIS_Family <- toupper(as.character(j[j$rankname == "family", 4]))
            }
        }
    }

    # Unify the taxonomic information from ITIS
    ITIS_data <- data.frame(
        ITIS_Present, ITIS_is_valid, ITIS_id, ITIS_name,
        ITIS_Phylum, ITIS_N_syn, ITIS_syn,
        ITIS_Class, ITIS_Order, ITIS_Family, 
        ITIS_species_in_genus
    )

    return(ITIS_data)
}

ITIS_get_genus_data <- function(genus, n_times = 3) {
    TSN <- NULL
    t_4 <- 1

    while (is.null(TSN) && t_4 <= n_times) {
        try(TSN <- get_tsn_(genus, searchtype = "scientific", silent = TRUE))
        t_4 <- t_4 + 1
        print(TSN)
    }

    if (is.list(TSN) & is.null(TSN[[1]])) {
        TSN[[1]] <- data.frame()
    }

    if (nrow(TSN[[1]]) == 0) {
        try(TSN.two <- get_tsn(genus, searchtype = "scientific", silent = TRUE))
        if (attr(TSN.two, "match") == "found") {
            tsn2 <- itis_acceptname(TSN.two[[1]], silent = TRUE)
            ## 1. Add a new row (tibble is currently empty)
            TSN[[1]] <- TSN[[1]] |>
                add_row(
                    tsn            = TSN.two[[1]],
                    scientificName = genus,
                    commonNames    = NA,
                    nameUsage      = "valid"
                )
        } else if (attr(TSN.two, "match") == "not found") {
            TSN <- list(c(NULL))
        }
    }

if (is.null(TSN[[1]])) {
        ITIS_Present <- "No"
        ITIS_id <- NA
        ITIS_name <- NA
        ITIS_is_valid <- NA
        ITIS_Phylum <- NA
        ITIS_Class <- NA
        ITIS_Order <- NA
        ITIS_Family <- NA

        ITIS_accept <- NA
        ITIS_syn <- NA
        ITIS_N_syn <- NA
        ITIS_species_in_genus <- NA
    } else {
        tsn <- TSN[[1]]

        tsn_n <- tsn[tsn$scientificName == genus & tsn$nameUsage == "valid", ]
        if (dim(tsn_n)[1]==0){
          tsn_n <- tsn[grepl(genus, tsn$scientificName) & tsn$nameUsage == "valid", ]
        }
        if (length(tsn_n$tsn) == 0) { # if the original species name is not present in the returned list of synonims we igore the results

            ITIS_Present <- "Wrong"
            ITIS_id <- NA
            ITIS_name <- NA
            ITIS_is_valid <- NA
            ITIS_Phylum <- NA
            ITIS_Class <- NA
            ITIS_Order <- NA
            ITIS_Family <- NA

            ITIS_accept <- NA
            ITIS_syn <- NA
            ITIS_N_syn <- NA
            ITIS_species_in_genus <- NA
        } else {
            ITIS_Present <- "Yes"
            ITIS_species_in_genus <- NA
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            AceptName <- NULL
            t_4 <- 1

            while (is.null(AceptName) && t_4 <= n_times) {
                try(AceptName <- itis_acceptname(tsn_n$tsn, silent = TRUE))
                t_4 <- t_4 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

            if (AceptName$submittedtsn == AceptName$acceptedtsn) { # is the submitted Or_name name accepted by ITIS
                ITIS_name <- genus
                ITIS_accept <- "yes"
            } else if (!is.na(AceptName$acceptedname)) {
                ITIS_name <- AceptName$acceptedname # if the name is not accepted, we use the accepted name
                ITIS_accept <- "yes"
            } else if (is.na(AceptName$acceptedname)) {
                ITIS_name <- "No"
                ITIS_accept <- "No"
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            Syn <- NULL
            t_5 <- 1

            while (is.null(Syn) && t_5 <= n_times) {
                try(Syn <- synonyms(tsn_n$tsn, db = "itis"))
                t_5 <- t_5 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

            # Syn_l is a list of possible codes with their respective names
            # We need to unify the list and extract the names

            Syn <- rbindlist(Syn, fill = TRUE)
            # count the number of FALSE and TRUE       #strsplit(, " "), length
            z <- as.character(is.na(Syn$syn_name))
            # if the length equals cero, there is no syn names

            if (length(z[z == "FALSE"]) == 0) {
                ITIS_N_syn <- NA
                ITIS_syn <- NA
            } else {
                z_b <- as.character(sapply(strsplit(Syn$syn_name, split = " "), length) == 1)
                Syn$syn_full <- paste(Syn$syn_name, Syn$syn_author)
                ITIS_N_syn <- length(z_b[z_b == TRUE]) # number of sp sinonims, Subsp excluded
                ITIS_syn <-  paste(Syn$syn_full, collapse = "; ")  # combine the names into a single string with all the sinonims
            }

            ITIS_id <- tsn_n$tsn
            if (is.null(ITIS_name) | is.na(ITIS_name)) {
                ITIS_name <- tsn_n$scientificName
            }
            ITIS_is_valid <- tsn_n$nameUsage

            # Get the upstream taxonomic information
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            j <- NULL
            t_13 <- 1
            t_31 <- 1
            t_32 <- 1
            while (is.null(j) && t_31 <= n_times) {
                try(j <- itis_hierarchy(tsn_n$tsn, what = "full", silent = TRUE))
                t_31 <- t_31 + 1
            }
            while (dim(j)[1] == 0 && t_32 <= n_times) {
                try(j <- itis_hierarchy(AceptName$acceptedtsn, what = "full", silent = TRUE))
                t_32 <- t_32 + 1
            }
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            if (dim(j)[1] > 0) {
                ITIS_Phylum <- toupper(as.character(j[j$rankname == "phylum", 4])) # Change from lower to upper chase
                ITIS_Class <- toupper(as.character(j[j$rankname == "class", 4]))
                ITIS_Order <- toupper(as.character(j[j$rankname == "order", 4]))
                ITIS_Family <- toupper(as.character(j[j$rankname == "family", 4]))

                if ( "species" %in% j$rankname){
                    ITIS_species_in_genus <- j %>%
                        filter(rankname == "species") %>%
                        pull(taxonname) %>%
                        paste(collapse = "; ")
                }
            }
        }
    }

    # Unify the taxonomic information from ITIS
    ITIS_data <- data.frame(
        ITIS_Present, ITIS_is_valid, ITIS_id, ITIS_name,
        ITIS_Phylum, ITIS_N_syn, ITIS_syn,
        ITIS_Class, ITIS_Order, ITIS_Family,
        ITIS_species_in_genus
    )

    return(ITIS_data)
}

# GBIF_get_data function
GBIF_get_data <- function(spp.x, species, n_times = 3) {
    
    GBIF_Present <- NA
    GBIF_id <- NA
    GBIF_name <- NA
    GBIF_Phylum <- NA
    GBIF_Class <- NA
    GBIF_Order <- NA
    GBIF_Family <- NA
    
    GBIF_Status <- NA
    GBIF_syn <- NA
    GBIF_N_syn <- NA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    key_1 <- NULL
    t_6 <- 1

    while (is.null(key_1) && t_6 <= n_times) {
        try(key_1 <- get_gbifid_(sci = spp.x)[[1]], silent = TRUE) # get the taxon key
        t_6 <- t_6 + 1
    }
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # For binomial species
    if (!is.na(species)){
    if (length(key_1) == 0) {
        GBIF_Present <- "No"
        GBIF_id <- NA

        GBIF_name <- NA

        GBIF_Phylum <- NA
        GBIF_Class <- NA
        GBIF_Order <- NA
        GBIF_Family <- NA

        GBIF_syn <- NA
        GBIF_N_syn <- NA
    } else {
        GBIF_id <- paste(key_1$specieskey[key_1$status == "ACCEPTED" & key_1$matchtype == "EXACT"], collapse = "-")
        if(GBIF_id == ""){
            max_conf = max(key_1$confidence)
            GBIF_id <- paste(key_1$usagekey[key_1$status == "ACCEPTED" & key_1$matchtype == "FUZZY" & key_1$confidence==max_conf], collapse = "-")
        }
        GBIF_name <- paste(unique(key_1$species[key_1$specieskey == GBIF_id]), collapse = "-")
        GBIF_Status = paste(unique(key_1$status[key_1$specieskey == GBIF_id]), collapse = "-")
        GBIF_Phylum <- toupper(paste(unique(key_1$phylum[key_1$specieskey == GBIF_id]), collapse = "-"))
        GBIF_Class <- toupper(paste(unique(key_1$class[key_1$specieskey == GBIF_id]), collapse = "-"))
        GBIF_Order <- toupper(paste(unique(key_1$order[key_1$specieskey == GBIF_id]), collapse = "-"))
        GBIF_Family <- toupper(paste(unique(key_1$family[key_1$specieskey == GBIF_id]), collapse = "-"))

        if (is.na(GBIF_name)) {
            GBIF_Present <- "No"
            GBIF_id <- NA

            GBIF_name <- NA

            GBIF_Phylum <- NA
            GBIF_Class <- NA
            GBIF_Order <- NA
            GBIF_Family <- NA
            
            GBIF_Status <- NA
            GBIF_syn <- NA
            GBIF_N_syn <- NA
        }

        if (GBIF_name == "" & !"SYNONYM" %in% key_1$status) {
            GBIF_Present <- "No"
            GBIF_id <- NA

            GBIF_name <- NA
         
            GBIF_Phylum <- NA
            GBIF_Class <- NA
            GBIF_Order <- NA
            GBIF_Family <- NA

            GBIF_Status <- NA
            GBIF_syn <- NA
            GBIF_N_syn <- NA
        } else {
            GBIF_Present <- "Yes"
            GBIF_syn <- paste(unique(key_1$species[key_1$status == "SYNONYM"]), collapse = ";")
            GBIF_syn_multi <- unique(key_1$species[key_1$status == "SYNONYM"])
            if(nchar(GBIF_id) == 0){
                GBIF_id <- paste(key_1$usagekey[key_1$status == "SYNONYM"], collapse = "-")
                GBIF_id_multi <- key_1$usagekey[key_1$status == "SYNONYM"]
            }
            if (length(key_1[key_1$status == "SYNONYM", 1]) == 0) {
                GBIF_N_syn <- NA
            } else {
                GBIF_N_syn <- length(key_1[key_1$status == "SYNONYM", 1])
            }
            if (!is.na(GBIF_syn) & GBIF_syn!="" & GBIF_name=="") {
             try(tax <- get_gbifid_(sci = GBIF_syn_multi[1])[[1]], silent = TRUE)
                if(dim(tax)[1]>0){
                    GBIF_name = paste(unique(tax$species[tax$status=="ACCEPTED" & tax$matchtype=="EXACT"]), collapse = "; ")
                    GBIF_Status <- paste(unique(tax$status[key_1$usagekey == GBIF_id_multi[1]]), collapse = "-")
                    GBIF_Phylum <- toupper(paste(unique(tax$phylum[tax$usagekey == GBIF_id_multi[1]]), collapse = "-"))
                    GBIF_Class <- toupper(paste(unique(tax$class[tax$usagekey == GBIF_id_multi[1]]), collapse = "-"))
                    GBIF_Order <- toupper(paste(unique(tax$order[tax$usagekey == GBIF_id_multi[1]]), collapse = "-"))
                    GBIF_Family <- toupper(paste(unique(tax$family[tax$usagekey == GBIF_id_multi[1]]), collapse = "-"))
                }
            }
        }
    }
    # For genus level
    } else {
        if (length(key_1) == 0) {
            GBIF_Present <- "No"
            GBIF_id <- NA
            
            GBIF_name <- NA
            
            GBIF_Phylum <- NA
            GBIF_Class <- NA
            GBIF_Order <- NA
            GBIF_Family <- NA
            
            GBIF_syn <- NA
            GBIF_N_syn <- NA
        } else {
            GBIF_Present <- "Yes"
            GBIF_id <- paste(key_1$usagekey[key_1$status == "ACCEPTED" & key_1$matchtype == "HIGHERRANK"], collapse = "-")
            GBIF_name <- paste(unique(key_1$genus[key_1$usagekey == GBIF_id]), collapse = "-")
            GBIF_Phylum <- toupper(paste(unique(key_1$phylum[key_1$usagekey == GBIF_id]), collapse = "-"))
            GBIF_Class <- toupper(paste(unique(key_1$class[key_1$usagekey == GBIF_id]), collapse = "-"))
            GBIF_Order <- toupper(paste(unique(key_1$order[key_1$usagekey == GBIF_id]), collapse = "-"))
            GBIF_Family <- toupper(paste(unique(key_1$family[key_1$usagekey == GBIF_id]), collapse = "-"))
            GBIF_syn = NA
            GBIF_N_syn = NA
        }
    }
    GBif_data <- data.frame(
        GBIF_Present, GBIF_id, GBIF_name,
        GBIF_N_syn, GBIF_syn, GBIF_Phylum, GBIF_Status,
        GBIF_Class, GBIF_Order, GBIF_Family
    )

    return(GBif_data)
}

# Resolve capitalization
CapSp <- function(x) {
    s <- strsplit(x, " ") %>% unlist()

    if (length(s) > 1) {
        paste(paste0(toupper(substring(s[1], 1, 1)), substring(s[1], 2)),
            paste(tolower(s[-1]), collapse = " "),
            sep = " "
        )
    } else {
        paste0(toupper(substring(s[1], 1, 1)), substring(s[1], 2))
    }
}

## ----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\#

retrieve_syns_new <- function(spp_name, # [Character] The species name from which to collect taxonomic information
                          n_times = 3, # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                          Gbif = FALSE # [Logical] Should we check Gbif for a taxonomic macthing of the species
) {
    # 0. Load the packages
    list.of.packages <- c(
        "tidyr", "rredlist", "taxize", "data.table", "stringr",
        "rgbif", "raster", "data.table", "dplyr"
    )

    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    if (length(new.packages)) install.packages(new.packages)

    lapply(list.of.packages, require, character.only = TRUE)
    rm(list.of.packages, new.packages)
    rm(spp.x, correct_name, y.d, res, res_max, accepted, genus, species, infra, elems,
       IUCN_data, ITIS_data, Tax_dat, Spp_syn, correct_name, iucn_unique, itis_unique)
    
    # a. Check the species names, if the name is binomial it runs the query to collect the taxonomic information
    # a.1 Removes special characters and fix the capitalization of the name
    spp.x <- stringr::str_trim(spp_name) %>%
        gsub(pattern = "[[:punct:]]", replacement = " ") %>%
        stringr::str_trim("both") %>%
        gsub(pattern = "  ", replacement = "") # Remove white spaces at the start and end of the string

    spp.x <- CapSp(spp.x)

    spp.x_copy = spp.x # save a copy of the original name
    # a.2 Check if the name is related with a species/sub-specie or other taxonomic class (Class, Order, Family)
    #     by analyzing the number of terms of the character string
    correct_name <- NULL
    t_11 <- 1

    while (is.null(correct_name) && t_11 <= n_times) {
        res <- tryCatch(
            {
                res <- gna_verifier(names = spp.x, all_matches = TRUE) #|>
                # dplyr::select(matchedCanonicalSimple) |>
                # dplyr::pull()
            },
            error = function(e) NULL # on error return NULL
        )
        print(res)
        if (!is.null(res) && length(res) > 0) {
            correct_name <- res
        }

        t_11 <- t_11 + 1
    }
    rm(t_11)

    # Summarize the results of the name checking and correction
    if (!is.null(correct_name)) {
        if (nrow(correct_name) != 0) {
            y.d <- cbind(spp.x, correct_name[, colnames(correct_name) %in%
                c(
                    "dataSourceTitleShort", "cardinalityScore",
                    "sortScore", "matchedCanonicalFull",
                    "taxonomicStatus", "curatedDataScore"
                )])

            names(y.d)[1] <- "or_name"
            if ("Accepted" %in% y.d$taxonomicStatus & !"Synonym" %in% y.d$taxonomicStatus) {
                accepted <- y.d %>% filter(taxonomicStatus == "Accepted")
                if (length(unique(accepted$matchedCanonicalFull > 1))) {
                    # If there are more than one accepted name, we select the one with the highest sortScore
                    res_max <- accepted %>% slice_max(sortScore, n = 1)
                    spp.x <- unique(res_max$matchedCanonicalFull)
                } else if (length(unique(accepted$matchedCanonicalFull)) == 1) {
                    # If there is only one accepted name, we use it
                    spp.x <- unique(accepted$matchedCanonicalFull)
                }
            } else if ("Synonym" %in% y.d$taxonomicStatus) {
                y2d = y.d %>% filter(taxonomicStatus != "N/A" )
                res_max <- y2d %>% slice_max(sortScore, n = 1)
                spp.x <- unique(res_max$matchedCanonicalFull)
            } else if (unique(y.d$taxonomicStatus=="N/A")){
                res_max <- y.d %>% slice_max(sortScore, n = 1)
                spp.x <- unique(res_max$matchedCanonicalFull)
            }

            ifelse(y.d$matchedCanonicalFull == y.d$or_name,
                y.d$Status <- "Correct",
                y.d$Status <- "Incorrect"
            )
        }
    } else {
        y.d <- data.frame(
            or_name = spp.x,
            y.d = NA,
            Status = "Not_found",
            data_source_title = NA,
            score = NA
        )
    }

    # b.Use the corrected or original name to look for taxonomic data----
    #   b.1. Get the basic data from the IUCN red list----
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Check if this is a genus-level name (single word) or species-level name (two or more words)
    if(is.na(spp.x) | is.null(spp.x)){
        spp.x = spp.x_copy
    }
    
    elems <- spp.x %>%
        strsplit(split = " ") %>%
        unlist()
    genus <- elems[1]
    species <- elems[2]
    infra <- NULL
    if (length(elems) == 3) {
        # If the name is a sub-species, we remove the sub-species name
        infra <- elems[3]
    }

    if (!is.na(species)) {
        g <- NULL
        t_11 <- 1
        # Get the taxonomic information from the IUCN red list
        while (is.null(g) && t_11 <= n_times) {
            try(g <- rl_species(genus = genus, species = species, infra = infra), silent = TRUE)
            t_11 <- t_11 + 1
        }
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        b <- NULL
        t_2 <- 1

        while (is.null(b) && t_2 <= n_times) {
            ## 1. try latest assessment -------------------------------------
            b <- tryCatch(
                rredlist::rl_species_latest(
                    genus = genus,
                    species = species,
                    infra = infra
                )$taxon$synonyms,
                error = function(e) NULL
            )

            ## 2. fall back if nothing returned -----------------------------
            if (is.null(b)) {
                b <- get_synonyms(genus = genus, species = species, infra = infra)
                IUCN_latest <- NA
                IUCN_date <- NA
                IUCN_Category <- NA
            }
            print(b)

            t_2 <- t_2 + 1
        }
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        # b.2 Get IUCN data using the dedicated function----
        IUCN_data <- IUCN_get_data(g = g, b = b, genus = genus, species = species, infra = infra, n_times = n_times)

        # b.3. Get ITIS data using the dedicated function----
        ITIS_data <- ITIS_get_species_data(spp.x = spp.x, n_times = n_times)

        # Should we retrieve synonim information from GBIF?
        if (Gbif == TRUE) {
            GBif_data <- GBIF_get_data(spp.x = spp.x,species = species, n_times = n_times)
        }

        # C. Return the Taxonomic information
        if (Gbif == TRUE) {
            Tax_dat <- cbind(Or_name = spp.x, IUCN_data, ITIS_data, GBif_data)

            Spp_syn <- c(spp.x, Tax_dat[, colnames(Tax_dat) %in% c(
                "IUCN_name", "IUCN_syn", "ITIS_name",
                "ITIS_syn", "GBIF_name", "GBIF_syn"
            )]) %>%
                paste(collapse = ";") %>%
                strsplit(split = ";") %>%
                unlist()

            Spp_syn <- Spp_syn[-c(Spp_syn %>% grep(pattern = "NA"))]
            Spp_syn <- Spp_syn[!duplicated(Spp_syn)]
            
            # initialise
            correct_name <- NA_character_
            
            # Extract unique names (remove NA, keep first if multiple)
            iucn_unique <- na.omit(unique(IUCN_data$IUCN_name))
            itis_unique <- na.omit(unique(ITIS_data$ITIS_name))
            gbif_unique <- na.omit(unique(GBif_data$GBIF_name))
            
            # Apply robust logic
            if (length(iucn_unique) > 0 && length(itis_unique) > 0 && identical(iucn_unique[1], itis_unique[1])) {
                correct_name <- iucn_unique[1]
            } else if (length(itis_unique) > 0 && (length(iucn_unique) == 0)) {
                correct_name <- itis_unique[1]
            } else if (length(iucn_unique) > 0) {
                correct_name <- iucn_unique[1]
            } else if (length(gbif_unique) > 0) {
                correct_name <- gbif_unique[1]
            }
             
            return(list(
                Submitted_name = spp_name,
                Spp_syn = Spp_syn,
                taxon_level = "species",
                correct_name = correct_name,
                IUCN_spp = IUCN_data$IUCN_name,
                Or_name = spp.x,
                TaxDat = Tax_dat
            ))
        } else {
            print("No GBIF data")
            if (all(!is.na(IUCN_data$IUCN_latest))) {
                IUCN_data2 <- IUCN_data %>% filter(IUCN_latest == TRUE)
            } else {
                IUCN_data2 <- IUCN_data
            }

            Tax_dat <- cbind(Or_name = spp.x, IUCN_data2, ITIS_data)
            

            ## --- usage inside your pipeline --------------------------------------------
            Spp_syn <- c(
                spp.x,
                Tax_dat[, colnames(Tax_dat) %in% c(
                    "IUCN_name", "IUCN_syn",
                    "ITIS_name", "ITIS_syn",
                    "GBIF_name", "GBIF_syn"
                )]
            ) %>%
                clean_synonyms()

            # initialise
            correct_name <- NA_character_

            iucn_unique <- unique(IUCN_data2$IUCN_name)
            itis_unique <- unique(ITIS_data$ITIS_name)
            gbif_unique <- unique(GBif_data$GBIF_name)
            # both names present and identical
            if (!is.na(iucn_unique) && !is.na(itis_unique) &&
                identical(iucn_unique, itis_unique)) {
                correct_name <- iucn_unique

                # only ITIS name present
            } else if (is.na(iucn_unique) && !is.na(itis_unique)) {
                correct_name <- itis_unique
                # only IUCN name present
            } else if (!is.na(iucn_unique) && is.na(itis_unique)) {
                correct_name <- iucn_unique
            } else if (!is.na(gbif_unique) && is.na(itis_unique) && is.na(iucn_unique)) {
                correct_name <- gbif_unique
            } 
            

            return(list(
                Submitted_name = spp_name,
                correct_name = correct_name,
                taxon_level = "species",
                Spp_syn = Spp_syn,
                #ITIS_species_in_genus = if ("ITIS_species_in_genus" %in% names(ITIS_data)) ITIS_data$ITIS_species_in_genus else NA,
                IUCN_spp = IUCN_data2$IUCN_name,
                Or_name = spp.x,
                TaxDat = Tax_dat
            ))
        }
    } else {
        print("Genus-level name")
       
        # There is no IUCN data for genus-level names
        IUCN_Present <- "No"
        IUCN_id <- NA
        IUCN_name <- NA
        IUCN_Phylum <- NA
        IUCN_Class <- NA
        IUCN_Order <- NA
        IUCN_Family <- NA
        IUCN_Category <- NA
        IUCN_latest <- NA
        IUCN_date <- NA
        IUCN_Category <- NA
        IUCN_N_syn <- NA
        IUCN_status <- NA
        IUCN_syn <- NA
        IUCN_data <- unique(data.frame(
        IUCN_Present, IUCN_id, IUCN_name, IUCN_latest, IUCN_date,
        IUCN_Category, IUCN_N_syn, IUCN_syn, IUCN_status,
        IUCN_Phylum, IUCN_Class,
        IUCN_Order, IUCN_Family))
        
        #Extract the ITIS data for the genus
        ITIS_data <- ITIS_get_genus_data(genus = genus, n_times = n_times)
        
        # GBIF DATA
        if (Gbif == TRUE) {
            GBif_data <- GBIF_get_data(spp.x = spp.x, species = species, n_times = n_times)
        } else {
            GBif_data <- data.frame(
                GBIF_Present = "No",
                GBIF_id = NA,
                GBIF_name = NA,
                GBIF_N_syn = NA,
                GBIF_syn = NA,
                GBIF_Phylum = NA,
                GBIF_Class = NA,
                GBIF_Order = NA,
                GBIF_Family = NA
            )
        }
        
        
        Tax_dat <- cbind(Or_name = spp.x, IUCN_data, ITIS_data, GBif_data)
            
        ## --- usage inside your pipeline --------------------------------------------
        Spp_syn <- c(
                spp.x,
                Tax_dat[, colnames(Tax_dat) %in% c(
                    "IUCN_name", "IUCN_syn",
                    "ITIS_name", "ITIS_syn",
                    "GBIF_name", "GBIF_syn"
                )]
            ) %>%
                clean__genus_synonyms()

            # initialise
            correct_name <- NA_character_

            # Extract unique names (remove NA, keep first if multiple)
            iucn_unique <- na.omit(unique(IUCN_data$IUCN_name))
            itis_unique <- na.omit(unique(ITIS_data$ITIS_name))
            gbif_unique <- na.omit(unique(GBif_data$GBIF_name))
            
            # Apply robust logic
            if (length(iucn_unique) > 0 && length(itis_unique) > 0 && identical(iucn_unique[1], itis_unique[1])) {
                correct_name <- iucn_unique[1]
            } else if (length(itis_unique) > 0 && (length(iucn_unique) == 0)) {
                correct_name <- itis_unique[1]
            } else if (length(iucn_unique) > 0) {
                correct_name <- iucn_unique[1]
            } else if (length(gbif_unique) > 0) {
                correct_name <- gbif_unique[1]
            }

            return(list(
                Submitted_name = spp_name,
                correct_name = correct_name,
                taxon_level = "genus",
                Spp_syn = Spp_syn,
                #ITIS_species_in_genus = if ("ITIS_species_in_genus" %in% names(ITIS_data)) ITIS_data$ITIS_species_in_genus else NA,
                IUCN_spp = IUCN_data$IUCN_name,
                Or_name = spp.x,
                TaxDat = Tax_dat
            ))

    }
    # Message
    # print("Taxonomic search done")
}
