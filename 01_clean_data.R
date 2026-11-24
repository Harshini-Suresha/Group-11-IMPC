#Used for datasets to build normalised relational tables and export the cleaned datasets plus the audit logs
library(data.table)
library(stringr)

#Input folder and Output folder for reproductiblity 
RAW_DIR   <- "group11_data"
OUT_DIR   <- "clean_data"
dir.create(OUT_DIR, showWarnings = FALSE)
#This code searches the group11_data.zip file then it unzips the file and moves it to the RAW_DIR folder.
if (!dir.exists(RAW_DIR) && file.exists("group11_data.zip")) unzip("group11_data.zip") 

#Custom logging function that helped us monitor long running R scripts on the HPC
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# 0. Helpers
# characters that were randomly injected into values

# trim removes leading/trailing whitespace.
NOISE_CHARS <- "[!@#$%^&*()_+.\\- ]"

## norm_key removes punctuation/spacing and lowercases values for tolerant matching.
norm_key  <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

## Store cleaning-log rows in an environment so helper functions can append to it.
trim      <- function(x) gsub("^\\s+|\\s+$", "", x)
## log_change compares raw and cleaned vectors and records only values that changed.
na_if_lit <- function(x) { x <- trim(x); x[x %in% c("", "NA", "N/A", "NULL", "na")] <- NA; x }

# A row is considered unchanged when both values are NA or both non-missing values match.
# ch identifies changed records; the log stores ID, field, before/after values, and rule.
log_env <- new.env(); log_env$rows <- list()

## log_change compares raw and cleaned vectors and records only values that changed.
log_change <- function(analysis_id, field, raw, clean, rule) {
  same <- (is.na(raw) & is.na(clean)) | (!is.na(raw) & !is.na(clean) & raw == clean)
  ch <- !same
  if (any(ch)) log_env$rows[[length(log_env$rows) + 1]] <- data.table( ## Matching uses norm_key, so punctuation/case differences do not prevent a match.
    analysis_id = analysis_id[ch], field_name = field,                  # Original missing values remain missing rather than being converted to a vocabulary label.
    raw_value = raw[ch], clean_value = clean[ch], rule = rule)
}

# vocab_clean maps values to the canonical spelling in a supplied controlled vocabulary.
vocab_clean <- function(x, vocab) {
  k <- norm_key(x); vk <- norm_key(vocab)
  out <- vocab[match(k, vk)]
  out[is.na(x)] <- NA
  out
}

# most frequent spelling within a grouping where the ties are broken by `prefer` regex
consensus <- function(value, group, prefer = NULL) {
  dt <- data.table(value, group)[!is.na(value) & !is.na(group), .N, by = .(group, value)]
  if (!is.null(prefer)) dt[, ok := grepl(prefer, value)] else dt[, ok := TRUE]
  setorder(dt, group, -ok, -N)
  dt[, .SD[1], by = group][, .(group, consensus = value)]
}

# Blank lines are skipped, but their expected positional field is recorded as an issue.
# 1. Import raw analysis files 
# Split each line into a candidate key and value at the first comma.
# Remove surrounding quotes/whitespace from the value while preserving NA text.
FIELDS <- c("gene_accession_id", "gene_symbol", "allele_symbol", "strain_name",
            "parameter_stable_id", "p_value", "parameter_name", "mouse_life_stage",
            "lacz_blood", "lacz_brain", "lacz_colon", "lacz_eye", "lacz_heart",
            "lacz_kidney", "lacz_liver", "lacz_lung", "lacz_skeletal_muscle",
            "lacz_skin", "lacz_vascular_system", "analysis_id")
# Some exports concatenate adjacent labels so it records the issue and uses the final valid label.

LACZ <- grep("^lacz_", FIELDS, value = TRUE)

files <- sort(list.files(RAW_DIR, pattern = "\\.csv$", full.names = TRUE))

# If a label is invalid, fall back to the expected field at that line position.
msg("Importing", length(files), "raw analysis files")

parse_file <- function(f) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  
  # Assign the parsed value to its field. Return the record plus source and issue metadata.
  rec <- setNames(rep(NA_character_, length(FIELDS)), FIELDS)
  issues <- character(0)
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (trim(ln) == "") {                     # blank line -> positional field missing
      if (i <= length(FIELDS)) issues <- c(issues, paste0("blank_line:", FIELDS[i]))
      next
    }
    key_raw <- sub(",.*$", "", ln)
    val_raw <- if (grepl(",", ln)) sub("^[^,]*,", "", ln) else NA_character_
    val <- if (is.na(val_raw)) NA_character_ else gsub('^"|"$', "", trim(val_raw))
    if (!is.na(val_raw) && trim(val_raw) == "NA") val <- "NA"   # unquoted NA
    keys <- strsplit(gsub('^"|"$', "", key_raw), '""', fixed = TRUE)[[1]]
    keys <- keys[keys != ""]
    if (length(keys) > 1) {                     # merged labels, e.g. "lacz_blood""lacz_brain"
      issues <- c(issues, paste0("merged_labels:", paste(keys, collapse = "+")))
      if (keys[length(keys)] %in% FIELDS) { rec[keys[length(keys)]] <- val; next }
      keys <- ""
    }
    key <- if (length(keys) == 1) keys else ""
    if (!key %in% FIELDS) {                     # malformed label -> use position
      if (i <= length(FIELDS)) {
        issues <- c(issues, paste0("label_repaired_by_position:", FIELDS[i]))
        key <- FIELDS[i]
      } else next
    }
    rec[key] <- val
  }
  c(rec, source_file = basename(f), import_issue = paste(issues, collapse = ";")) # Return the parsed record with its source filename and import issue metadata.
}

# Parse all input files and combine their records into one data.table.
raw <- rbindlist(lapply(files, function(f) as.list(parse_file(f))))

# Convert empty import issue descriptions into proper missing values.
raw[import_issue == "", import_issue := NA]

# Print the number of parsed records and records with structural issues.
msg("Parsed", nrow(raw), "records;", sum(!is.na(raw$import_issue)), "with structural issues")

# Export the parsed records as a tab-separated file.
fwrite(raw, file.path(OUT_DIR, "raw_records_parsed.tsv"), sep = "\t")

# Create a separate copy for subsequent data cleaning.
d <- copy(raw)
# Generate an analysis ID by removing the optional numeric suffix and CSV extension
d[, file_analysis_id := sub("(_\\d+)?\\.csv$", "", source_file)]
# Create a row identifier using the original source filename.
d[, row_uid := sprintf("%s", source_file)]

# 2. Field-level cleaning
# 2.1 analysis_id : file name is authoritative (SOP: 'an' + 10 digits)
## Store the original analysis IDs before making changes.
aid_raw <- d$analysis_id
# Convert literal missing-value strings into proper NA values.
d[, analysis_id := na_if_lit(analysis_id)]
d[is.na(analysis_id) | analysis_id != file_analysis_id, analysis_id_note :=
    fifelse(is.na(analysis_id), "missing_in_file_restored_from_filename",
            "mismatch_with_filename_restored_from_filename")]
d[, analysis_id := file_analysis_id] # Replace all analysis IDs with the authoritative IDs extracted from filenames.
# Log the original and corrected analysis IDs for data-cleaning traceability.
log_change(d$analysis_id, "analysis_id", aid_raw, d$analysis_id, "restore_from_filename")

# 2.2 gene_accession_id
# Store the original gene accession IDs before cleaning.
acc_raw <- d$gene_accession_id
# Convert accession IDs to uppercase, remove unwanted characters and standardise missing values.
d[, gene_accession_id := toupper(gsub("[^A-Za-z0-9:]", "", na_if_lit(gene_accession_id)))]

# Validate accession IDs against the expected MGI format (MGI: followed by digits).
# Set invalid or non-matching accession IDs to NA.
d[!grepl("^MGI:[0-9]+$", gene_accession_id), gene_accession_id := NA]

# Log the original and cleaned accession IDs for traceability.
log_change(d$analysis_id, "gene_accession_id", acc_raw, d$gene_accession_id, "strip_noise_validate_MGI")

# 2.3 gene_symbol : evidence-ranked spelling per MGI accession 
#     order: (1) MGI symbol in the Disease Ontology report when its
#     normalised key agrees, (2) the spelling that equals an allele-symbol prefix,
#     (3) the most frequent spelling with standard nomenclature form.

# Store the original gene symbols before cleaning.
gs_raw <- d$gene_symbol
# Remove unwanted special characters and spaces while preserving legal '.' and '-' characters.
gs <- gsub("[!@#$%^&*()_+ ]", "", na_if_lit(gs_raw))     # '.' and '-' are legal in symbols
d[, gs_tmp := gs] # Create a temporary column containing the cleaned gene symbols.

# Read the Disease Ontology report to obtain reference gene symbols.
mgi_sym <- local({
  dl <- readLines("Disease_information.txt", warn = FALSE); dl <- dl[!grepl("^#", dl) & trim(dl) != ""]
  x <- fread(text = dl, sep = "\t", header = FALSE, colClasses = "character", quote = "", fill = TRUE)
  unique(x[, .(acc = trim(V8), sym = trim(V7))])[, .SD[1], by = acc]
})

# Extract allele-symbol prefixes by removing everything from '<' onwards.
# Clean unwanted special characters while preserving legal symbol punctuation.
al_prefix <- trim(sub("<.*$", "", gsub("[!@#$%^&*()_+ ]", "", na_if_lit(d$allele_symbol))))

# Count the frequency of each cleaned gene-symbol spelling for each MGI accession.
cand <- d[!is.na(gs_tmp) & !is.na(gene_accession_id), .N, by = .(acc = gene_accession_id, v = gs_tmp)]

# Identify gene symbols that match an allele-symbol prefix for the same MGI accession.
pref <- unique(data.table(acc = d$gene_accession_id, v = al_prefix)[!is.na(v)])[, is_prefix := TRUE]
# Identify gene symbols that match an allele-symbol prefix for the same MGI accession.
cand <- merge(cand, pref, by = c("acc", "v"), all.x = TRUE)
# Mark candidates without an allele-prefix match as FALSE.
cand[is.na(is_prefix), is_prefix := FALSE]

# Identify gene symbols following standard nomenclature patterns.
cand[, nomen := grepl("^[A-Z][a-z0-9]*([.-][A-Za-z0-9]+)*$|^[A-Z0-9]+[a-z0-9]*$", v)]

# Rank candidates by allele-prefix match, standard nomenclature and frequency.
# Higher-priority candidates are placed first.
setorder(cand, acc, -is_prefix, -nomen, -N)

# Select the highest-ranked gene-symbol spelling for each MGI accession.
best <- cand[, .SD[1], by = acc]

# Add the reference gene symbols from the Disease Ontology report.
best <- merge(best, mgi_sym, by = "acc", all.x = TRUE)

# Prefer the Disease Ontology spelling when its normalised key matches the selected candidate.
best[!is.na(sym) & norm_key(sym) == norm_key(v), v := sym]

# Map the selected gene symbol back to every record using its MGI accession.
cg <- best$v[match(d$gene_accession_id, best$acc)]

# Replace gene symbols with the selected spelling, falling back to cleaned original symbols.
# Remove the temporary working column after assignment.
d[, gene_symbol := fifelse(!is.na(cg), cg, gs_tmp)][, gs_tmp := NULL]

# Recover gene symbols that are missing across all records using allele-symbol prefixes.

# Count allele-prefix occurrences for each MGI accession and rank by frequency.
ap <- data.table(acc = d$gene_accession_id, v = al_prefix)[!is.na(v) & !is.na(acc), .N, by = .(acc, v)][order(acc, -N)][, .SD[1], by = acc]

# Identify records with missing gene symbols whose MGI accession has an available allele prefix.
need <- is.na(d$gene_symbol) & d$gene_accession_id %in% ap$acc

# Restore missing gene symbols using the most frequent allele prefix for that accession.
d[which(need), gene_symbol := ap$v[match(gene_accession_id, ap$acc)]]

# Report the number of recovered gene symbols and the number still missing.
msg("Gene symbols recovered from allele prefix:", sum(need), "records;", sum(is.na(d$gene_symbol)), "still missing")

# Log the original and final gene symbols to maintain a record of cleaning changes.
log_change(d$analysis_id, "gene_symbol", gs_raw, d$gene_symbol, "evidence_ranked_symbol_per_MGI_accession")
# 2.4 allele_symbol : remove noise, then consensus within (accession, normalised key)
al_raw <- d$allele_symbol
al <- gsub("[!@#$%^&*_+ ]", "", na_if_lit(al_raw))
# allele prefix must be the gene symbol: repair when normalised keys agree
apfx <- sub("<.*$", "", al); fix <- !is.na(al) & grepl("<", al) & norm_key(apfx) == norm_key(d$gene_symbol) & apfx != d$gene_symbol
al[which(fix)] <- paste0(d$gene_symbol[which(fix)], substring(al[which(fix)], nchar(apfx[which(fix)]) + 1))
d[, al_tmp := al][, al_key := paste(gene_accession_id, norm_key(al_tmp))]
cons_al <- consensus(d$al_tmp, fifelse(is.na(d$al_tmp), NA_character_, d$al_key),
                     prefer = "^[A-Za-z0-9.-]+<[^<>]+>$")
ca <- cons_al$consensus[match(d$al_key, cons_al$group)]
d[, allele_symbol := fifelse(!is.na(ca) & !is.na(al_tmp), ca, al_tmp)]
d[, c("al_tmp", "al_key") := NULL]
# validation: brackets balanced and allele prefix equals gene symbol
d[, allele_valid := !is.na(allele_symbol) &
     grepl("^[^<>]+<[^<>]+>$", allele_symbol) &
     sub("<.*$", "", allele_symbol) == gene_symbol]
# where the allele is missing but the gene has exactly one observed allele, impute it
one_allele <- d[!is.na(allele_symbol), .(n = uniqueN(allele_symbol), a = allele_symbol[1]), by = gene_accession_id][n == 1]
d[, a := one_allele$a[match(gene_accession_id, one_allele$gene_accession_id)]]
d[is.na(allele_symbol) & !is.na(a), `:=`(allele_symbol = a, allele_note = "imputed_single_allele_for_gene")]
d[, a := NULL]
log_change(d$analysis_id, "allele_symbol", al_raw, d$allele_symbol,
           "noise_strip_consensus_or_single_allele_imputation")

# 2.5 strain / life stage / lacZ: Standardise values using controlled vocabularies defined in the SOP.

# Define the accepted mouse strain names according to the SOP.
STRAINS <- c("B6Brd;B6Dnk;B6N-Tyr<c-Brd>", "C57BL/6Brd-Tyr<c-Brd>", "C57BL/6Dnk",
             "C57BL/6JTyr;C57BL/6N", "C57BL/6N", "C57BL/6NCrl", "C57BL/6NJ",
             "C57BL/6NTac", "C57BL/6NTacDen")

# Define the accepted mouse developmental and life stages.
STAGES <- c("E9.5", "E12.5", "E15.5", "E18.5", "Early adult", "Middle aged adult", "Late adult")

# Define the accepted lacZ expression categories.
LACZ_V <- c("No data", "ambiguous", "expression", "imageOnly", "no expression", "tissue not available")

# Iterate through strain, life-stage and lacZ fields.
for (f in c("strain_name", "mouse_life_stage", LACZ)) {
  
  # Select the appropriate controlled vocabulary for the current field.
  vocab <- switch(f, strain_name = STRAINS, mouse_life_stage = STAGES, LACZ_V)
  
  # Store the original values before cleaning.
  r <- d[[f]]
  
  # Convert literal missing-value strings to NA and match values against the accepted vocabulary.
  v <- vocab_clean(na_if_lit(r), vocab)
  
  # Replace the original field with the standardised values.
  set(d, j = f, value = v)
  
  # Log the original and cleaned values for data-cleaning traceability.
  log_change(d$analysis_id, f, r, v, "controlled_vocabulary_match")
}

# 2.6 p_value: Convert p-values to numeric format and identify invalid values.

# Store the original p-values before cleaning.
pv_raw <- d$p_value
# Convert p-values from character to numeric after handling literal missing values.
# Suppress warnings generated by values that cannot be converted to numbers.
d[, p_value := suppressWarnings(as.numeric(na_if_lit(p_value)))]
# Identify missing and out-of-range p-values.
# Valid p-values must fall between 0 and 1 (inclusive).
d[, p_value_flag := fcase(
  is.na(p_value), "missing",
  p_value < 0 | p_value > 1, "out_of_range_[0,1]",
  default = NA_character_
)]
# Log the original and converted p-values for data-cleaning traceability.
log_change(d$analysis_id, "p_value", pv_raw, as.character(d$p_value), "trim_to_numeric")


# 3. Reference data: Import and prepare parameter and procedure reference tables.

# Import the IMPC parameter description reference file, treating all columns as character data.
pdesc <- fread("IMPC_parameter_description.txt", colClasses = "character", strip.white = TRUE)

# Import the IMPC procedure reference file, treating all columns as character data.
proc <- fread("IMPC_procedure.txt", colClasses = "character", strip.white = TRUE)

# Remove leading and trailing whitespace from all columns in the parameter reference table.
pdesc[, names(pdesc) := lapply(.SD, trim)]

# Remove leading and trailing whitespace from all columns in the procedure reference table.
proc[, names(proc) := lapply(.SD, trim)]

# Record the original number of parameter reference rows before removing duplicates.
n0 <- nrow(pdesc)

# Remove exact duplicate rows from the parameter reference table.
pdesc <- unique(pdesc)

# Report the number of exact duplicate parameter records removed.
msg("Parameter reference: removed", n0 - nrow(pdesc), "exact duplicates")

# Check that parameter original IDs are unique in both reference tables.
# Stop execution if duplicate IDs are detected.
stopifnot(!anyDuplicated(pdesc$impcParameterOrigId), !anyDuplicated(proc$impcParameterOrigId))

# Merge parameter descriptions with procedure information using the common original parameter ID.
# Add suffixes to distinguish columns originating from each reference table.
proc_ref <- merge(pdesc, proc, by = "impcParameterOrigId", suffixes = c("_param", "_proc"))

# Verify that the merge retained all parameter description records.
# Stop execution if the merged row count differs from the original parameter table.
stopifnot(nrow(proc_ref) == nrow(pdesc))

# Prepare the procedure reference table by standardising uninformative descriptions.

# Replace generic procedure descriptions with empty strings.
proc_ref[description_proc %in% c("Housing_and_husbandry", "experimental_design"), description_proc := ""]

# Create a procedure-level table containing one row per unique procedure name.
# Select the first non-empty description and extract procedure codes from parameter IDs.
procedure <- proc_ref[, .(
  description = {
    x <- unique(description_proc)
    x <- x[x != ""]
    if (length(x)) x[1] else NA_character_
  },
  procedure_code = paste(unique(sub("^IMPC_([A-Z]+)_.*$", "\\1", parameterId)), collapse = ",")
), by = .(procedure_name = name_proc)]

# Sort the procedure table alphabetically by procedure name.
setorder(procedure, procedure_name)
# Generate a sequential numeric identifier for each procedure.
procedure[, procedure_id := .I]
# Create a reference table containing one row per unique stable parameter ID.
# Retain the parameter name, description, associated procedure and mandatory status.
ref_params <- proc_ref[, .(
  parameter_name_ref = name_param[1],
  parameter_description = description_param[1],
  procedure_name = name_proc[1],
  is_mandatory = as.integer(toupper(isMandatory[1]) == "TRUE")
), by = .(parameter_stable_id = parameterId)]

# Preserve individual IMPC parameter versions using their original parameter IDs.
# Retain both original and stable IDs, descriptions and mandatory status.
param_version <- proc_ref[, .(
  impc_parameter_orig_id = as.integer(impcParameterOrigId),
  parameter_stable_id = parameterId,
  name = name_param,
  description = description_param,
  is_mandatory = as.integer(toupper(isMandatory) == "TRUE")
)]


# 3.1 parameter_stable_id repair

# Store original parameter IDs before cleaning.
pid_raw <- d$parameter_stable_id

# Convert literal missing-value strings into proper NA values.
pid <- na_if_lit(pid_raw)

# Create a lookup table mapping normalised reference IDs to their official spellings.
ref_key <- setNames(ref_params$parameter_stable_id, norm_key(ref_params$parameter_stable_id))

# Normalise parameter IDs to allow matching despite formatting differences.
k <- norm_key(pid)

# Match dataset parameter IDs against the official IMPC reference IDs.
pid_ref <- unname(ref_key[k])

# For IDs not found in the reference, identify consistent spellings within groups.
# Use a regular expression to prioritise IDs following the expected legacy ID structure.
d[, pid_tmp := toupper(pid)][, pid_k := k]
cons_pid <- consensus(d$pid_tmp[is.na(pid_ref)], d$pid_k[is.na(pid_ref)],
                      prefer = "^[A-Z]+(-[A-Z])*(-[A-Z]+)?_[A-Z0-9]+_[0-9]{3}_[0-9]{3}$")

# Prefer official reference IDs; otherwise use the consensus spelling for non-reference IDs.
d[, pid_clean := fifelse(!is.na(pid_ref), pid_ref,
                         cons_pid$consensus[match(pid_k, cons_pid$group)])]

# Preserve originally missing parameter IDs instead of introducing inferred values.
d[is.na(pid), pid_clean := NA]

# Record whether each parameter ID came from the IMPC reference or consensus-based correction.
d[, parameter_id_source := fifelse(!is.na(pid_ref), "IMPC_reference",
                                   fifelse(!is.na(pid_clean), "consensus_spelling", NA_character_))]


# 3.2 parameter_name repair and recovery of missing IDs

# Store original parameter names before cleaning.
pn_raw <- d$parameter_name

# Remove leading/trailing whitespace and collapse multiple internal spaces into single spaces.
d[, parameter_name := gsub("\\s+", " ", na_if_lit(parameter_name))]

# Build a lookup table linking normalised parameter names to their stable IDs.
name2id <- ref_params[, .(ids = list(unique(parameter_stable_id))), by = .(nk = norm_key(parameter_name_ref))]

# Retain only parameter names associated with exactly one unique stable ID to avoid ambiguous matches.
name2id <- name2id[lengths(ids) == 1, .(nk, id = unlist(ids))]

# Identify records with missing parameter IDs but available parameter names.
miss <- is.na(d$pid_clean) & !is.na(d$parameter_name)

# Recover parameter IDs by matching normalised parameter names against the unique reference lookup.
rec <- name2id$id[match(norm_key(d$parameter_name[miss]), name2id$nk)]

# Restore missing IDs where a unique reference match is available and record their recovery source.
d[which(miss)[!is.na(rec)], `:=`(pid_clean = rec[!is.na(rec)], parameter_id_source = "recovered_from_unique_parameter_name")]

# Store the final cleaned parameter IDs and remove temporary working columns.
d[, parameter_stable_id := pid_clean][, c("pid_tmp", "pid_k", "pid_clean") := NULL]

# Log original and cleaned parameter IDs for data-cleaning traceability.
log_change(d$analysis_id, "parameter_stable_id", pid_raw, d$parameter_stable_id,
           "reference_key_match_or_consensus_or_unique_name_recovery")

# canonical parameter name: reference name where available, else consensus experimental name
cons_pn <- consensus(d$parameter_name, d$parameter_stable_id)
pn_canon <- merge(cons_pn, ref_params[, .(parameter_stable_id, parameter_name_ref)],
                  by.x = "group", by.y = "parameter_stable_id", all.x = TRUE)
pn_canon[, canonical := fifelse(!is.na(parameter_name_ref), parameter_name_ref, consensus)]
d[, parameter_name_clean := pn_canon$canonical[match(parameter_stable_id, pn_canon$group)]]
d[is.na(parameter_name_clean), parameter_name_clean := parameter_name]
log_change(d$analysis_id, "parameter_name", pn_raw, d$parameter_name_clean, "canonical_parameter_name")


# 4. Duplicate analyses (files suffixed _1)

# Define the biological and experimental fields used to compare duplicate analyses.
cmp_cols <- c("gene_accession_id", "gene_symbol", "allele_symbol", "strain_name",
              "parameter_stable_id", "p_value", "mouse_life_stage", LACZ)

# Initialise a counter to track the number of cleaning changes per record.
d[, n_changes := 0L]

# Combine all recorded cleaning changes into a single audit table.
lg <- rbindlist(log_env$rows)

# Count the number of logged changes for each analysis ID.
chg <- lg[, .N, by = analysis_id]

# Identify analyses appearing more than once in the dataset.
dups <- d[analysis_id %in% d[duplicated(analysis_id), analysis_id]]

# Compare duplicate records across selected fields to determine whether they are identical after cleaning.
dup_report <- dups[, .(files = paste(source_file, collapse = ";"),
                       identical_after_cleaning = uniqueN(.SD) == 1),
                   by = analysis_id, .SDcols = cmp_cols]

# Report the number of duplicate analyses and how many are identical after cleaning.
msg("Duplicate analysis ids:", nrow(dup_report), "; identical after cleaning:",
    sum(dup_report$identical_after_cleaning))

# Export the duplicate analysis report for inspection.
fwrite(dup_report, file.path(OUT_DIR, "duplicate_analysis_report.tsv"), sep = "\t")

# Keep the original unsuffixed file as the primary record and use duplicate files to recover missing values.

# Identify duplicate files using their numeric filename suffix.
d[, is_duplicate_file := grepl("_\\d+\\.csv$", source_file)]

# Initialise a list to store conflicting values between original and duplicate records.
conflicts <- list()

# Process only duplicate analyses that are not identical after cleaning.
for (a in dup_report[identical_after_cleaning == FALSE, analysis_id]) {
  
  # Identify the original unsuffixed record and its corresponding duplicate record.
  io <- which(d$analysis_id == a & !d$is_duplicate_file)
  id <- which(d$analysis_id == a & d$is_duplicate_file)
  
  # Compare each selected field between the original and duplicate records.
  for (cc in cmp_cols) {
    vo <- d[[cc]][io]
    vd <- d[[cc]][id[1]]
    
    # If the original value is missing but the duplicate has a value, recover it from the duplicate.
    if (is.na(vo) && !is.na(vd)) {
      set(d, io, cc, vd)
      
      # Record the recovered value in the cleaning audit log.
      log_env$rows[[length(log_env$rows) + 1]] <- data.table(analysis_id = a, field_name = cc, raw_value = NA_character_,
                                                             clean_value = as.character(vd), rule = "filled_from_duplicate_source_file")
      
      # If both records contain different non-missing values, record the conflict without overwriting either value.
    } else if (!is.na(vo) && !is.na(vd) && vo != vd) {
      conflicts[[length(conflicts) + 1]] <- data.table(analysis_id = a, field = cc, original = as.character(vo), duplicate = as.character(vd))
    }
  }
}

# Report the number of duplicate analyses processed and conflicting values detected.
msg("Duplicate pairs merged; value conflicts:", length(conflicts))

# Export conflicting values separately for manual investigation.
if (length(conflicts)) fwrite(rbindlist(conflicts), file.path(OUT_DIR, "duplicate_value_conflicts.tsv"), sep = "\t")

# Preserve duplicate filename and analysis ID information for reference.
dup_files <- d[is_duplicate_file == TRUE, .(source_file, analysis_id)]

# Add the duplicate comparison result to the duplicate file tracking table.
dup_files <- merge(dup_files, dup_report[, .(analysis_id, identical_after_cleaning)], by = "analysis_id")

# Remove suffixed duplicate files, retaining only the original unsuffixed records.
d <- d[is_duplicate_file == FALSE]

# Verify that every remaining analysis ID is unique.
stopifnot(!anyDuplicated(d$analysis_id))

# 5.1 Disease Ontology report (MGI)

# Read the Disease Ontology reference file, preserving each line as text.
dl <- readLines("Disease_information.txt", warn = FALSE)

# Remove comment lines beginning with '#' and blank lines.
dl <- dl[!grepl("^#", dl) & trim(dl) != ""]

# Read the remaining tab-separated data.
# All columns are imported as character to preserve IDs and their formatting.
dis_all <- fread(text = dl, sep = "\t", header = FALSE,
                 colClasses = "character", quote = "", fill = TRUE)

# Assign meaningful column names to the eight columns.
setnames(dis_all, c("human_symbol", "human_name", "hgnc_id", "do_name", "doid",
                    "mouse_genotype_ids", "mouse_symbol", "mgi_gene_id"))

# Remove leading and trailing whitespace from every column.
dis <- dis_all[, lapply(.SD, trim)]

# Count rows before removing exact duplicate records.
n0 <- nrow(dis)

# Remove rows that are identical across all columns.
dis <- unique(dis)

# Report how many exact duplicates were removed.
msg("Disease file:", n0, "rows,", n0 - nrow(dis), "exact duplicates removed")

# Validate that disease, human gene and mouse gene IDs follow their expected formats.
stopifnot(all(grepl("^DOID:[0-9]+$", dis$doid)),
          all(grepl("^HGNC:[0-9]+$", dis$hgnc_id)),
          all(grepl("^MGI:[0-9]+$", dis$mgi_gene_id)))

# Create a unique human gene table containing each HGNC ID and its annotations.
human_gene <- unique(dis[, .(hgnc_id, human_symbol, human_name)])

# Check that each HGNC ID occurs only once in the human gene table.
stopifnot(!anyDuplicated(human_gene$hgnc_id))

# Create a disease table with one disease name for each Disease Ontology ID.
disease <- dis[, .(disease_name = do_name[1]), by = doid]

# Create unique links between human genes and diseases.
human_gene_disease <- unique(dis[, .(hgnc_id, doid)])

# Create unique links between mouse genes and their human orthologues.
mouse_human_ortholog <- unique(dis[, .(mgi_gene_id, hgnc_id)])

# Split multiple mouse genotype IDs into individual rows.
# Each genotype ID is linked to its corresponding human gene and disease.
disease_model_genotype <- dis[mouse_genotype_ids != "",
                              .(mgi_genotype_id = unlist(
                                strsplit(mouse_genotype_ids, "|", fixed = TRUE))),
                              by = .(hgnc_id, doid)]

# Remove whitespace from genotype IDs and eliminate duplicate records.
disease_model_genotype <- unique(
  disease_model_genotype[, mgi_genotype_id := trim(mgi_genotype_id)]
)


# 5.2 Mouse genes: IMPC genes plus disease-report genes

# Extract unique mouse genes from the IMPC dataset, excluding missing gene IDs.
g_impc <- unique(d[!is.na(gene_accession_id),
                   .(mgi_gene_id = gene_accession_id, gene_symbol)])
# Ensure that each MGI gene ID occurs only once in the IMPC gene table.
stopifnot(!anyDuplicated(g_impc$mgi_gene_id))
# Extract mouse gene IDs and symbols from the Disease Ontology report.
g_dis <- unique(dis[, .(mgi_gene_id, gene_symbol = mouse_symbol)])
# Combine IMPC genes with disease-report genes that are not already in IMPC.
# in_impc_data = 1 means the gene came from IMPC; 0 means it came only from the disease report.
gene <- rbind(g_impc[, .(mgi_gene_id, gene_symbol, in_impc_data = 1L)],
              g_dis[!mgi_gene_id %in% g_impc$mgi_gene_id,
                    .(mgi_gene_id, gene_symbol, in_impc_data = 0L)])
# Keep one record per MGI gene ID.
gene <- gene[, .SD[1], by = mgi_gene_id]
# Match genes present in both sources to compare their gene symbols.
sym_check <- merge(g_impc, g_dis, by = "mgi_gene_id",
                   suffixes = c("_impc", "_mgi"))
# Report how many genes were compared and how many have different symbols.
msg("Gene symbols cross-checked with MGI disease report:", nrow(sym_check),
    "; disagreements:",
    sum(sym_check$gene_symbol_impc != sym_check$gene_symbol_mgi))
# Export genes whose IMPC and MGI disease-report symbols disagree.
fwrite(sym_check[gene_symbol_impc != gene_symbol_mgi],
       file.path(OUT_DIR, "gene_symbol_disagreements.tsv"), sep = "\t")
# Extract unique allele records, linking each allele symbol to its mouse gene ID.
allele <- unique(d[!is.na(allele_symbol),
                   .(allele_symbol, mgi_gene_id = gene_accession_id, allele_valid)])
# Collapse repeated allele records and mark an allele as valid if any record says it is valid.
allele <- allele[, .(allele_valid = as.integer(any(allele_valid))),
                 by = .(allele_symbol, mgi_gene_id)]
# Sort alleles alphabetically and assign sequential numeric IDs.
setorder(allele, allele_symbol)
allele[, allele_id := .I]
# Create a strain lookup table with a numeric ID for each controlled-vocabulary strain.
strain <- data.table(strain_id = seq_along(STRAINS), strain_name = STRAINS)
# Create a life-stage lookup table and flag embryonic stages.
life_stage <- data.table(
  life_stage_id = seq_along(STAGES),
  life_stage_name = STAGES,
  is_embryonic = as.integer(grepl("^E", STAGES))
)
# Create a tissue lookup table by removing the 'lacz_' prefix from tissue names.
lacz_tissue <- data.table(
  tissue_id = seq_along(LACZ),
  tissue_name = sub("^lacz_", "", LACZ)
)
# Create a lookup table assigning numeric IDs to LacZ call labels.
lacz_call <- data.table(call_id = seq_along(LACZ_V), call_label = LACZ_V)


# 5.3 Parameters: every ID used in the experiments, enriched with the reference
exp_ids <- unique(d[!is.na(parameter_stable_id), .(parameter_stable_id)])
parameter <- merge(exp_ids, ref_params, by = "parameter_stable_id", all.x = TRUE)
parameter <- rbind(parameter, ref_params[!parameter_stable_id %in% exp_ids$parameter_stable_id], fill = TRUE)
parameter[, parameter_name := pn_canon$canonical[match(parameter_stable_id, pn_canon$group)]]
parameter[is.na(parameter_name), parameter_name := parameter_name_ref]
parameter[, in_impc_reference := as.integer(!is.na(procedure_name))]
# infer procedure for centre-specific IDs from the procedure code (e.g. JAX_OFD_* -> Open Field)
code2proc <- procedure[!grepl(",", procedure_code), .(procedure_code, procedure_name)]
parameter[, code := fifelse(grepl("^[A-Z-]+_[A-Z]{3}_[0-9]{3}_[0-9]{3}$", parameter_stable_id),
                            sub("^[A-Z-]+_([A-Z]{3})_.*$", "\\1", parameter_stable_id), NA_character_)]
parameter[, procedure_link := fifelse(in_impc_reference == 1, "IMPC_reference", NA_character_)]
parameter[is.na(procedure_name) & code %in% code2proc$procedure_code,
          `:=`(procedure_name = code2proc$procedure_name[match(code, code2proc$procedure_code)],
               procedure_link = "inferred_from_procedure_code")]
parameter[, procedure_id := procedure$procedure_id[match(procedure_name, procedure$procedure_name)]]
parameter <- parameter[, .(parameter_stable_id, parameter_name, parameter_description,
                           procedure_id, procedure_link, is_mandatory, in_impc_reference)]
stopifnot(!anyDuplicated(parameter$parameter_stable_id))

# 5.4 Parameter groupings (many-to-many). Rules use the parameter name,
#     description and procedure; documented in parameter_group.rule_description.
# Join parameter records with procedure names using procedure_id.
# all.x = TRUE retains parameters even when no matching procedure is found.
pg_src <- merge(parameter, procedure[, .(procedure_id, procedure_name)], by = "procedure_id", all.x = TRUE)
# Combine parameter name and description into lowercase text for keyword matching.
pg_src[, txt := tolower(paste(parameter_name, parameter_description, sep = " | "))]
# Store the procedure name in lowercase; replace missing procedure names with "".
pg_src[, pr  := tolower(fifelse(is.na(procedure_name), "", procedure_name))]
# Define parameter groups, their descriptions, classification rules and rule expressions.
GROUPS <- list(
  list("Weight", "Body, organ and tissue weights",
       "name/description contains 'weight' or 'mass'; procedures Body Weight, Heart Weight, Organ Weight",
       quote(grepl("weight|\\bmass\\b", txt) | pr %in% c("body weight", "heart weight", "organ weight"))),
  list("Images", "Parameters that are images or derived from imaging",
       "name/description contains 'image'/'imaging'/'x-ray'/'microct'/'opt'; procedures X-ray, MicroCT*, OPT*",
       quote(grepl("imag|x-ray|xray|micro ?ct|\\bopt\\b", txt) | grepl("^x-ray|^microct|^opt", pr))),
  list("Brain", "Brain and central nervous system",
       "name/description contains brain, cerebr*, cerebell*, hippocamp*, cortex, ventricle(s) of brain, neur*, spinal cord; procedure Auditory Brain Stem Response",
       quote(grepl("brain|cerebr|cerebell|hippocamp|cortex|neur|spinal cord|hypothalam|thalam|pituitar|medulla", txt) |
               pr == "auditory brain stem response")),
  list("Blood & clinical chemistry", "Haematology, clinical chemistry and circulating analytes",
       "procedures Hematology, Clinical Chemistry, Insulin Blood Level; or name contains blood/plasma/serum",
       quote(pr %in% c("hematology", "clinical chemistry", "insulin blood level") | grepl("blood|plasma|serum", txt))),
  list("Cardiovascular", "Heart and vascular function / morphology",
       "procedures Echo, Electrocardiogram (ECG), Heart Weight; or name contains heart/cardi*/aort*/ventric*/atri*/vascular",
       quote(pr %in% c("echo", "electrocardiogram (ecg)", "heart weight") |
               grepl("heart|cardi|aort|vascular|\\bqrs\\b|\\bqt|\\bpr interval|ejection|stroke volume", txt))),
  list("Eye & vision", "Eye morphology and vision",
       "procedure Eye Morphology; or name contains eye/retina/lens/cornea/iris/pupil/optic",
       quote(pr == "eye morphology" | grepl("\\beye|retina|\\blens\\b|cornea|\\biris\\b|pupil|optic|vitreous|conjunctiv", txt))),
  list("Skeleton & bone", "Skeletal morphology and bone density",
       "procedure X-ray; or name contains bone/skull/femur/tibia/vertebra/rib/digit/BMD/skeletal",
       quote(pr == "x-ray" | grepl("bone|skull|femur|tibia|vertebra|\\bribs?\\b|digit|\\bbmd\\b|skelet|pelvis|scapula|clavicle|mandible|teeth|tooth", txt))),
  list("Behaviour & neurological", "Behavioural and neurological screens",
       "procedures Open Field, Light-Dark Test, Fear Conditioning, Y-maze, Acoustic Startle/PPI, Grip Strength, Combined SHIRPA, Formalin Test",
       quote(pr %in% c("open field", "light-dark test", "fear conditioning", "y-maze",
                       "acoustic startle and pre-pulse inhibition (ppi)", "grip strength",
                       "combined shirpa and dysmorphology", "formalin test") |
               grepl("startle|locomot|freezing|grip|anxiety|rotarod|behavio", txt))),
  list("Metabolism", "Energy balance, glucose handling and body composition",
       "procedures Indirect Calorimetry, IPGTT, Insulin Blood Level, Body Composition (DEXA); or name contains glucose/insulin/fat/lean/oxygen consumption/respiratory exchange",
       quote(pr %in% c("indirect calorimetry", "intraperitoneal glucose tolerance test (ipgtt)",
                       "insulin blood level", "body composition (dexa lean/fat)") |
               grepl("glucose|insulin|\\bfat\\b|lean|oxygen consumption|respiratory exchange|energy expenditure|triglycer|cholesterol", txt))),
  list("Immune system", "Immunophenotyping and immune cell populations",
       "procedure Immunophenotyping; or name contains T cell/B cell/NK/CD4/CD8/CD45/lymph*/monocyte/neutrophil/macrophage/granulocyte",
       quote(pr == "immunophenotyping" | grepl("t cell|b cell|\\bnk\\b|cd4|cd8|cd45|cd19|lymph|monocyte|neutrophil|macrophage|granulocyte|eosinophil|basophil|dendritic|immun", txt))),
  list("Embryo & development", "Embryonic / placental morphology and viability",
       "procedures Gross Morphology Embryo/Placenta*, Viability*, Embryo LacZ, MicroCT/OPT of embryos",
       quote(grepl("embryo|placenta|viability|microct|^opt", pr) | grepl("embryo|placent|yolk sac|somite", txt)))
)
parameter_group <- rbindlist(lapply(seq_along(GROUPS), function(i)
  data.table(group_id = i, group_name = GROUPS[[i]][[1]], group_description = GROUPS[[i]][[2]],
             rule_description = GROUPS[[i]][[3]])))
parameter_group_member <- rbindlist(lapply(seq_along(GROUPS), function(i)
  pg_src[eval(GROUPS[[i]][[4]]), .(group_id = i, parameter_stable_id)]))
msg("Parameter group sizes:"); print(parameter_group_member[, .N, by = group_id][parameter_group, on = "group_id"][, .(group_name, N)])

# 5.5 Analysis fact table + lacZ (long format)
analysis <- d[, .(analysis_id, mgi_gene_id = gene_accession_id,
                  allele_id = allele$allele_id[match(paste(allele_symbol, gene_accession_id),
                                                     paste(allele$allele_symbol, allele$mgi_gene_id))],
                  strain_id = strain$strain_id[match(strain_name, strain$strain_name)],
                  life_stage_id = life_stage$life_stage_id[match(mouse_life_stage, life_stage$life_stage_name)],
                  parameter_stable_id, p_value, p_value_flag, parameter_id_source,
                  source_file, import_issue,
                  analysis_id_note)]
# analyses with identical content (gene, allele, strain, stage, parameter, p-value)
# but different ids are retained; later ones point to the first of their group
analysis[, dup_key := fifelse(is.na(p_value) | is.na(parameter_stable_id), NA_character_,
                              paste(mgi_gene_id, allele_id, strain_id, life_stage_id, parameter_stable_id, p_value))]
setorder(analysis, analysis_id)
analysis[!is.na(dup_key), first_id := analysis_id[1], by = dup_key]
analysis[, content_duplicate_of := fifelse(!is.na(first_id) & first_id != analysis_id, first_id, NA_character_)]
analysis[, c("dup_key", "first_id") := NULL]
msg("Content-identical repeat analyses flagged:", sum(!is.na(analysis$content_duplicate_of)))
lacz <- melt(d[, c("analysis_id", LACZ), with = FALSE], id.vars = "analysis_id",
             variable.name = "tissue", value.name = "call", variable.factor = FALSE)
lacz <- lacz[!is.na(call), .(analysis_id,
                             tissue_id = lacz_tissue$tissue_id[match(sub("^lacz_", "", tissue), lacz_tissue$tissue_name)],
                             call_id = lacz_call$call_id[match(call, lacz_call$call_label)])]

# 5.6 Bespoke task: Ensembl homology (wide BioMart export -> long, species-generic)

# Import the Ensembl BioMart file.
# All columns are character to preserve gene IDs and other text values.
hml <- fread("bespokeTask_hml.txt", sep = "\t",
             colClasses = "character", quote = "", na.strings = NULL)

# Remove leading and trailing whitespace from all columns.
hml[, names(hml) := lapply(.SD, trim)]

# Record the number of rows before restructuring.
n_hml_rows <- nrow(hml)

# Rename the first column to identify the mouse Ensembl gene ID.
setnames(hml, 1, "mouse_ensembl_gene_id")

# Rename the mouse gene name column.
setnames(hml, "Gene name", "mouse_gene_name")

# Define a function that extracts homology information for one species.
species_cols <- function(sp) {
  
  # Find the five columns corresponding to the selected species:
  # Ensembl gene ID, gene name, last common ancestor,
  # homology type and percentage identity.
  cols <- grep(paste0("^", sp, " gene stable ID$|^", sp, " gene name$|^Last common ancestor with ", sp,
                      "$|^", sp, " homology type$|^%id. target ", sp),
               names(hml), value = TRUE)
  # Confirm that exactly five expected columns were identified.
  stopifnot(length(cols) == 5)
  # Extract the mouse gene ID and the five species-specific columns.
  x <- hml[, c("mouse_ensembl_gene_id", cols), with = FALSE]
  # Rename columns consistently, regardless of the target species.
  setnames(x, c("mouse_ensembl_gene_id", "homolog_ensembl_gene_id",
                "homolog_gene_name", "last_common_ancestor",
                "homology_type", "pct_identity"))
  # Add a column identifying the target species.
  x[, species := sp]
  # Retain only records with a non-empty homologue Ensembl ID.
  x[homolog_ensembl_gene_id != ""]
}
# Identify all target species from columns ending in 'gene stable ID'.
SPECIES <- sub(" gene stable ID$", "",
               grep(" gene stable ID$", names(hml), value = TRUE))
# Exclude the generic 'Gene' column, which is not a target species.
SPECIES <- setdiff(SPECIES, "Gene")
# Apply the extraction function to every target species and combine the results.
hom_long <- rbindlist(lapply(SPECIES, species_cols))
# Record the number of rows after converting the data into long format.
n_long <- nrow(hom_long)
# Remove duplicate records generated by the wide BioMart export.
hom_long <- unique(hom_long)
# Report the original wide row count, long row count and unique relationship count.
msg("Homology: wide rows", n_hml_rows,
    "| long rows", n_long,
    "| unique pairs", nrow(hom_long))
# Ensure that each mouse–homologue gene ID pair occurs only once.
stopifnot(!anyDuplicated(
  hom_long[, .(mouse_ensembl_gene_id, homolog_ensembl_gene_id)]
))
# Convert percentage identity values from character to numeric.
hom_long[, pct_identity := as.numeric(pct_identity)]

# Create a species lookup table with numeric IDs and Ensembl gene prefixes.
sp_tbl <- data.table(
  species_id = seq_len(length(SPECIES) + 1),
  species_name = c("Mouse", SPECIES),
  ensembl_gene_prefix = c(
    "ENSMUSG",
    sapply(SPECIES, function(s)
      sub("[0-9]+$", "",
          hom_long[species == s, homolog_ensembl_gene_id][1]))
  )
)

# Extract unique mouse Ensembl gene IDs and their associated gene names.
# Only mouse genes represented in the homology table are retained.
mouse_ens <- unique(
  hml[mouse_ensembl_gene_id %in% hom_long$mouse_ensembl_gene_id,
      .(ensembl_gene_id = mouse_ensembl_gene_id,
        gene_name = mouse_gene_name)]
)

# For each mouse Ensembl gene ID, select the first non-empty gene name.
# If no non-empty name exists, assign NA.
mouse_ens <- mouse_ens[
  , .(gene_name = {
    x <- unique(gene_name[gene_name != ""])
    if (length(x)) x[1] else NA_character_
  }),
  by = ensembl_gene_id
]

# Extract unique homologous Ensembl gene IDs and their target species.
homolog_ens <- unique(
  hom_long[, .(ensembl_gene_id = homolog_ensembl_gene_id, species)]
)

# For each homologous gene, select the first non-empty gene name.
# Merge these names back using the Ensembl gene ID.
homolog_ens <- merge(
  homolog_ens,
  hom_long[
    , .(gene_name = {
      x <- unique(homolog_gene_name[homolog_gene_name != ""])
      if (length(x)) x[1] else NA_character_
    }),
    by = .(ensembl_gene_id = homolog_ensembl_gene_id)
  ],
  by = "ensembl_gene_id"
)

# Combine mouse genes and homologous genes into one Ensembl gene table.
# Mouse genes receive species_id = 1.
# Homologous genes receive the corresponding species ID from sp_tbl.
ensembl_gene <- rbind(
  mouse_ens[, .(ensembl_gene_id, species_id = 1L, gene_name)],
  homolog_ens[
    , .(ensembl_gene_id,
        species_id = sp_tbl$species_id[
          match(species, sp_tbl$species_name)
        ],
        gene_name)
  ]
)

# Ensure that each Ensembl gene ID occurs only once in the final table.
stopifnot(!anyDuplicated(ensembl_gene$ensembl_gene_id))
# link mouse Ensembl genes to MGI genes through the official symbol (exact, case-sensitive)
gene_ensembl_xref <- merge(mouse_ens[!is.na(gene_name)], gene[, .(mgi_gene_id, gene_symbol)],
                           by.x = "gene_name", by.y = "gene_symbol")[, .(mgi_gene_id, ensembl_gene_id, match_method = "exact_symbol")]
ht_tbl <- data.table(homology_type_id = seq_along(unique(hom_long$homology_type)),
                            homology_type = sort(unique(hom_long$homology_type)))
tx_tbl <- data.table(taxon_id = seq_along(unique(hom_long$last_common_ancestor)),
                    taxon_name = sort(unique(hom_long$last_common_ancestor)))
homology <- hom_long[, .(mouse_ensembl_gene_id, homolog_ensembl_gene_id,
                         homology_type_id = ht_tbl$homology_type_id[match(homology_type, ht_tbl$homology_type)],
                         lca_taxon_id = tx_tbl$taxon_id[match(last_common_ancestor, tx_tbl$taxon_name)],
                         pct_identity)]
msg("Bespoke homology: mouse genes with >=1 homolog:", nrow(mouse_ens),
    "| linked to MGI genes:", nrow(gene_ensembl_xref),
    "| mouse genes without any homolog excluded:", uniqueN(hml$mouse_ensembl_gene_id) - nrow(mouse_ens))


# Combine all individual cleaning-log entries into one table.
cleaning_log <- rbindlist(log_env$rows)

# Retain log entries for valid analysis IDs.
# Note: '| TRUE' makes this condition always TRUE, so no rows are actually filtered.
cleaning_log <- cleaning_log[analysis_id %in% analysis$analysis_id | TRUE]

# Assign a sequential unique ID to every cleaning-log entry.
cleaning_log[, log_id := .I]

# Add an explanatory note to duplicate source-file records.
# These files are excluded from the main analysis table.
dup_files[, note := "duplicate source file of an existing analysis; excluded from analysis table"]

# Collect all normalised and reference tables into one named list.
# Each list element represents a separate output table.
out <- list(
  procedure = procedure[, .(procedure_id, procedure_name, procedure_code, description)],
  parameter = parameter,
  parameter_version = param_version,
  parameter_group = parameter_group,
  parameter_group_member = parameter_group_member,
  gene = gene,
  allele = allele[, .(allele_id, allele_symbol, mgi_gene_id, allele_valid)],
  strain = strain,
  life_stage = life_stage,
  lacz_tissue = lacz_tissue,
  lacz_call = lacz_call,
  analysis = analysis,
  analysis_lacz = lacz,
  duplicate_source_file = dup_files,
  human_gene = human_gene,
  disease = disease,
  human_gene_disease = human_gene_disease,
  mouse_human_ortholog = mouse_human_ortholog,
  disease_model_genotype = disease_model_genotype,
  species = sp_tbl,
  ensembl_gene = ensembl_gene,
  gene_ensembl_xref = gene_ensembl_xref,
  homology_type = ht_tbl,
  taxon = tx_tbl,
  homology = homology,
  data_cleaning_log = cleaning_log[
    , .(log_id, analysis_id, field_name, raw_value, clean_value, rule)
  ]
)

# Export every table in the output list as an individual TSV file.
# NA values are represented as \N and quotation marks are disabled.
for (n in names(out)) {
  fwrite(out[[n]],
         file.path(OUT_DIR, paste0(n, ".tsv")),
         sep = "\t", na = "\\N", quote = FALSE)
}

# Save all output tables together as one RDS file.
# This preserves their R data.table structures and data types.
saveRDS(out, file.path(OUT_DIR, "clean_tables.rds"))

# Generate a summary table containing the name and row count of each output table.
summ <- data.table(table = names(out), rows = sapply(out, nrow))
fwrite(summ, file.path(OUT_DIR, "table_row_counts.tsv"), sep = "\t")
print(summ)

# Count cleaning operations by field and cleaning rule.
# Sort the results from the most frequent operation to the least frequent.
fwrite(
  cleaning_log[, .N, by = .(field_name, rule)][order(-N)],
  file.path(OUT_DIR, "cleaning_summary.tsv"), sep = "\t"
)

# Print the total number of cleaning-log entries for each field.
print(cleaning_log[, .N, by = .(field_name)][order(-N)])

# Report remaining missing values in important analysis-table columns.
msg("Remaining missing values in analysis table:")

# Count NA values separately for each selected column.
print(
  analysis[, lapply(.SD, function(x) sum(is.na(x))),
           .SDcols = c("mgi_gene_id", "allele_id", "strain_id",
                       "life_stage_id", "parameter_stable_id", "p_value")]
)

# Record the R version, platform and loaded package information.
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "R_sessionInfo.txt"))
msg("Done.")