# ============================================================================
# upset_groups.R  —  UpSet plot of feature detection per experimental group
# ============================================================================
#   1. Groups are read from the ATTRIBUTE_group row of the file.
#   2. Samples sharing a group are replicates of each other.
#   3. A feature is PRESENT in a group if it is detected in >= n replicates of
#      THAT group (absolute threshold). Otherwise it is absent in that group.
#   4. Every feature is evaluated in every group -> presence/absence matrix.
#   5. That matrix feeds the UpSet plot (one set per group).
#
# HOW TO USE
#   1. Upload your file.
#   2. Change INPUT_FILE below.
#   3. Run the whole script. Output lands in 'results/'.
#   4. Define the replicates threshold (n) in MIN_REPLICATES 
#   5. Define groups color in SUFFIX_COLORS 
#   6. Set order on the y axis
# ============================================================================


# ══════════════════════════ CONFIGURATION ══════════════════════════════════

INPUT_FILE <- "metaboanalyst_HILIC_POS.csv"   # <<<<<  YOUR FILE NAME

# Number of replicates a feature must be detected in to count as present in a
# group. Absolute threshold: 2 means 2 replicates regardless of how many the
# group has (2 of 3 and 2 of 4 both qualify).
MIN_REPLICATES <- 2

# Row of the file that defines the groups.
GROUP_ROW <- "ATTRIBUTE_group"

# How many intersections to show (largest first).
N_INTERSECTIONS <- 25


# ── COLORS ──────────────────────────────────────────────────────────────────
# The color depends on the group in ATTRIBUTE. 
# This example is from the “medium_resin” extract, so modify it as needed.
SUFFIX_COLORS <- c(
  "GYM_HP20"    = "#6BAED6",
  "GYM_XAD7"    = "#2171B5",
  "GYMS_HP20"   = "#74C476",
  "GYMS_XAD7"   = "#238B45",
  "GYMSW_HP20"  = "#FDAE6B",
  "GYMSW_XAD7"  = "#E6550D"
)

# Set order on the y axis.
#   "color" = group by medium_resin (same color together), then by strain
#   "name"  = alphabetical (groups)
SET_ORDER <- "color"



# ══════════════════════════ OTHER CONFIGURATION ══════════════════════════════════

STRAIN_FIELDS <- 1

# Optional override: colors keyed by FULL group name. Overrides SUFFIX_COLORS.
# Leave NULL if not needed.
#   e.g. GROUP_COLORS <- c("KRD168_GYM_HP20" = "#000000")
GROUP_COLORS <- NULL
COLOR_INT_BARS <- "#3F5661"   # intersection bars 
COLOR_INT_TEXT <- "#3F5661"   # counts printed above those bars
COLOR_DOT_ON   <- "#3F5661"   # active dot in the matrix
COLOR_DOT_OFF  <- "#FFFFFF"   # inactive dot
COLOR_BG       <- "#FFFFFF"   # matrix panel background
STRIPE_ALPHA   <- 0.30        # transparency of the colored background stripes

# Set-size bar labels
SHOW_SET_SIZES  <- TRUE
COLOR_SIZE_TEXT <- "black"

# Group labels (y axis) colored by their group.
COLOR_LABELS <- TRUE
# ComplexUpset draws sets bottom-up. If label colors end up mismatched against
# the stripes, set this to TRUE.
FLIP_LABELS <- FALSE

###############################################################################
###############################################################################
###############################################################################
# ════════════════════════════════════════════════════════════════════════════
#                   Nothing below needs editing.
# ════════════════════════════════════════════════════════════════════════════
##############################################################################
###############################################################################
###############################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
  library(purrr); library(ggplot2); library(ComplexUpset); library(scales)
})
rule <- function() message(strrep("-", 68))


# ── 1. Read ─────────────────────────────────────────────────────────────────
if (!file.exists(INPUT_FILE))
  stop("Cannot find '", INPUT_FILE, "'. Files here: ",
       paste(list.files(pattern = "\\.(csv|txt|tsv)$", ignore.case = TRUE), collapse = ", "))

l1 <- readLines(INPUT_FILE, n = 1L, warn = FALSE)
delim <- names(which.max(c("\t" = str_count(l1, "\t"),
                           ";"  = str_count(l1, ";"),
                           ","  = str_count(l1, ","))))

raw <- read_delim(INPUT_FILE, delim = delim,
                  col_types = cols(.default = col_character()),
                  show_col_types = FALSE, name_repair = "unique")

id_col      <- names(raw)[1]                 # "Filename"
sample_cols <- names(raw)[-1]
message("File: ", INPUT_FILE, "  |  ", length(sample_cols), " samples")


# ── 2. Split the group row from the feature rows ────────────────────────────
ids  <- as.character(raw[[id_col]])
i_gr <- which(ids == GROUP_ROW)
if (length(i_gr) != 1L)
  stop("Expected exactly one '", GROUP_ROW, "' row, found ", length(i_gr),
       ".\n  Metadata rows present: ",
       paste(ids[str_starts(ids, "ATTRIBUTE_")], collapse = ", "))

group  <- as.character(raw[i_gr, sample_cols])   # one group per sample
groups <- sort(unique(group))

# ── Set order ───────────────────────────────────────────────────────────────
if (identical(SET_ORDER, "color")) {
  suffix_ord <- str_remove(groups, paste0("^([^_]+_){", STRAIN_FIELDS, "}"))
  strain_ord <- str_extract(groups, paste0("^([^_]+_){", STRAIN_FIELDS, "}"))
  # block order follows SUFFIX_COLORS; strain breaks ties
  groups <- groups[order(match(suffix_ord, names(SUFFIX_COLORS)), strain_ord)]
}

# feature rows: everything that is neither metadata nor blank
is_meta <- str_starts(ids, "ATTRIBUTE_") | is.na(ids) | ids == ""
feat_rows    <- raw[!is_meta, , drop = FALSE]


# ── 3. Detection matrix [features x samples] ────────────────────────────────
mat <- feat_rows[sample_cols] |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) |>
  as.matrix()
rownames(mat) <- make.unique(as.character(feat_rows[[id_col]]))

# detected = value present and > 0 (empty or 0 = not detected)
detected <- !is.na(mat) & mat > 0

rule()
message(sprintf("Features: %d  |  samples: %d  |  not detected: %.1f%%",
                nrow(mat), ncol(mat), 100 * mean(!detected)))


# ── 4. Replicates per group ─────────────────────────────────────────────────
n_rep    <- table(group)[groups]
required <- setNames(rep(as.integer(MIN_REPLICATES), length(groups)), groups)

rule()
message("GROUPS (", length(groups), ")  —  present if detected in >= ",
        MIN_REPLICATES, " replicates")
walk(groups, ~ message(sprintf("  %-24s %d replicates", .x, n_rep[[.x]])))
if (any(n_rep < MIN_REPLICATES))
  warning("Groups with fewer than ", MIN_REPLICATES,
          " replicates (they can never be present): ",
          paste(groups[n_rep < MIN_REPLICATES], collapse = ", "))


# ── 5. Presence per group (the core step) ───────────────────────────────────
#   For each feature and each group: in how many replicates was it detected?
#   Present if that count reaches the threshold for THAT group.
presence <- vapply(groups, function(g) {
  cols <- which(group == g)
  rowSums(detected[, cols, drop = FALSE]) >= required[[g]]
}, logical(nrow(mat)))
colnames(presence) <- groups
rownames(presence) <- rownames(mat)

# a feature present in no group contributes nothing to the UpSet
no_group <- rowSums(presence) == 0
presence <- presence[!no_group, , drop = FALSE]

rule()
message("Features present in >= 1 group: ", nrow(presence),
        "  (dropped ", sum(no_group), " that miss the threshold everywhere)")
message("Set sizes: ", min(colSums(presence)), " - ",
        max(colSums(presence)), " features")
message("Present in all ", length(groups), " groups (core): ",
        sum(rowSums(presence) == length(groups)))
message("Exclusive to a single group: ", sum(rowSums(presence) == 1))


# ── 6. Palette: one color per group ─────────────────────────────────────────
#   suffix = group name minus the first STRAIN_FIELDS fields.
suffix <- str_remove(groups, paste0("^([^_]+_){", STRAIN_FIELDS, "}"))
PAL <- setNames(unname(SUFFIX_COLORS[suffix]), groups)

if (!is.null(GROUP_COLORS)) {                      # literal override
  shared <- intersect(names(GROUP_COLORS), groups)
  PAL[shared] <- GROUP_COLORS[shared]
  extra <- setdiff(names(GROUP_COLORS), groups)
  if (length(extra))
    warning("GROUP_COLORS names groups that do not exist (ignored): ",
            paste(extra, collapse = ", "))
}

if (any(is.na(PAL)))
  stop("No color for: ", paste(groups[is.na(PAL)], collapse = ", "),
       "\n  Suffixes not found: ",
       paste(unique(suffix[is.na(PAL)]), collapse = ", "),
       "\n  Suffixes defined in SUFFIX_COLORS: ",
       paste(names(SUFFIX_COLORS), collapse = ", "),
       "\n  Check STRAIN_FIELDS (currently ", STRAIN_FIELDS, ") or add the suffix.")

rule()
message("COLORS")
walk(groups, ~ message(sprintf("  %-24s %s", .x, PAL[[.x]])))


# ── 7. UpSet ────────────────────────────────────────────────────────────────
d <- as.data.frame(presence)
d$feature_id <- rownames(presence)

# table ComplexUpset uses to paint the background stripes: it needs a column
# named `set` holding the set names.
# NOTE: the extra column must NOT be called `group` or `group_name` — those
# names collide with the internal data.frame ComplexUpset merges against.
stripe_data <- data.frame(set = groups, stripe_key = groups, stringsAsFactors = FALSE)

# order in which ComplexUpset draws the sets on the y axis (bottom-up)
label_colors <- if (FLIP_LABELS) rev(unname(PAL)) else unname(PAL)

p <- upset(
  d,
  intersect             = groups,
  n_intersections       = N_INTERSECTIONS,
  sort_intersections_by = "cardinality",   # with 18 sets, sorting by degree is uninformative
  sort_sets             = FALSE,           # honour the order of `groups`
  encode_sets           = FALSE,           # essential: otherwise set names become
  # numbers and the named color vectors
  # no longer match
  name                  = NULL,
  height_ratio          = 0.60,
  width_ratio           = 0.28,
  
  # background stripes: one color per group
  stripes = upset_stripes(
    mapping = aes(color = stripe_key),
    colors  = scales::alpha(PAL, STRIPE_ALPHA),
    data    = stripe_data
  ),
  
  base_annotations = list(
    "Features" = intersection_size(
      text = list(size = 2.4, vjust = -0.2, color = COLOR_INT_TEXT),
      fill = COLOR_INT_BARS
    ) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = COLOR_BG, color = NA),
            axis.line.y = element_line(linewidth = 0.3))
  ),
  
  # set-size bars: colored by group
  set_sizes = upset_set_size(geom = geom_bar(aes(fill = group, x = group),
                                             width = 0.6)) +
    (if (SHOW_SET_SIZES)
      geom_text(aes(label = after_stat(count)), stat = "count",
                hjust = 1, size = 3, colour = COLOR_SIZE_TEXT)
     else NULL) +
    expand_limits(y = max(colSums(presence)) * 1.12) +
    scale_fill_manual(values = PAL, guide = "none") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7)),
  
  matrix = intersection_matrix(geom    = geom_point(size = 1.9),
                               segment = geom_segment(linewidth = 0.35)) +
    scale_color_manual(values = c("TRUE" = COLOR_DOT_ON, "FALSE" = COLOR_DOT_OFF),
                       guide = "none",
                       na.value = COLOR_BG),   # <- the color of the matrix box
  
  # per-component themes: this is what lets us color the matrix y axis
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(
      text             = element_text(size = 8),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = COLOR_BG, color = NA),
      axis.text.y      = element_text(
        size   = 8,
        face   = if (COLOR_LABELS) "bold" else "plain",
        colour = if (COLOR_LABELS) label_colors else "black"
      )
    ),
    "overall_sizes" = theme(
      text             = element_text(size = 8),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = COLOR_BG, color = NA)
    )
  ))
) +
  labs(title = "",
       subtitle = sprintf(
         "n = %d features  |  >= %d replicates  |  top %d intersections",
         nrow(presence), MIN_REPLICATES, N_INTERSECTIONS)) +
  theme(plot.title = element_text(face = "bold", size = 12))


# ── 8. Save ─────────────────────────────────────────────────────────────────
dir.create("results", showWarnings = FALSE)
w <- 6 + 0.35 * length(groups)
h <- 5 + 0.22 * length(groups)
ggsave("results/upset_groups.png", p, width = w, height = h,
       dpi = 150, bg = "white", limitsize = FALSE)
ggsave("results/upset_groups.pdf", p, width = w, height = h, limitsize = FALSE)

# presence table. feature_id embeds ID / annotation (optional) / mz / RT.
out <- as_tibble(presence) |>
  mutate(feature_id = rownames(presence), .before = 1)
emb <- str_match(out$feature_id, "^(\\d+)(?:/(.+?))?/([0-9.]+)mz/([0-9.]+)min$")
if (mean(!is.na(emb[, 1])) > 0.9) {
  out <- out |>
    mutate(row_id = emb[, 2], annotation = emb[, 3],
           mz = as.numeric(emb[, 4]), rt_min = as.numeric(emb[, 5]),
           .after = feature_id)
  message("feature_id split into row_id / annotation / mz / rt_min  (",
          sum(!is.na(emb[, 3])), " annotated)")
}
out$n_groups <- rowSums(presence)
write_csv(out, "results/presence_by_group.csv")

write_csv(tibble(sample = sample_cols, group = group),
          "results/samples_and_groups.csv")

rule()
message("Done. Files in 'results/':")
message("  - upset_groups.png / .pdf")
message("  - presence_by_group.csv    (TRUE/FALSE per group + mz/RT/annotation)")
message("  - samples_and_groups.csv   (which sample fell in which group)")

print(p)
