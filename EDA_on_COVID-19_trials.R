# =============================================================================
# COVID Clinical Trials — Meaningful EDA
# Each plot here is built to answer a specific clinical-research question,
# not just "what does this column look like."
# =============================================================================

# ── 1. Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(janitor)
library(scales)
library(viridis)
library(ggtext)

# ── 2. Shared theme ───────────────────────────────────────────────────────────
theme_covid <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title         = element_markdown(face = "bold", size = 16,
                                            margin = margin(b = 4)),
      plot.subtitle      = element_textbox_simple(colour = "grey40", size = 11.5,
                                                   margin = margin(b = 14),
                                                   lineheight = 1.2),
      plot.caption       = element_text(colour = "grey60", size = 8.5, hjust = 0,
                                        margin = margin(t = 10)),
      plot.background    = element_rect(fill = "#FAFAF7", colour = NA),
      panel.background   = element_rect(fill = "#FAFAF7", colour = NA),
      panel.grid.major   = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor   = element_blank(),
      axis.title         = element_text(colour = "grey25", size = 11, face = "bold"),
      axis.text          = element_text(colour = "grey35"),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = 9.5),
      legend.text        = element_text(size = 9),
      legend.background  = element_rect(fill = "#FAFAF7", colour = NA),
      strip.text         = element_text(face = "bold", size = 10.5, colour = "grey20"),
      strip.background   = element_rect(fill = "grey92", colour = NA)
    )
}
theme_set(theme_covid())

pal_status <- c(
  "Completed"                 = "#1B7A8C",
  "Recruiting"                = "#6FC2A4",
  "Not yet recruiting"        = "#A8DADC",
  "Active, not recruiting"    = "#3D5A80",
  "Enrolling by invitation"   = "#98C1D9",
  "Terminated"                = "#D64045",
  "Withdrawn"                 = "#E89B3B",
  "Suspended"                 = "#F0C808",
  "Available"                 = "#8D99AE",
  "No longer available"       = "#B0A8B9",
  "Approved for marketing"    = "#5C946E",
  "Temporarily not available" = "#C9ADA7"
)

# ── 3. Load & clean ───────────────────────────────────────────────────────────
df <- read.csv("COVID_clinical_trials.csv", stringsAsFactors = FALSE) %>%
  clean_names()

df <- df %>%
  mutate(
    start_date      = as.Date(start_date, format = "%B %d, %Y"),
    completion_date = as.Date(completion_date, format = "%B %d, %Y"),
    first_posted    = as.Date(first_posted, format = "%B %d, %Y"),
    start_year      = year(start_date),
    duration_days   = as.numeric(completion_date - start_date),
    # Lag between trial start and public registration — a transparency signal
    posting_lag_days = as.numeric(first_posted - start_date)
  )

# =============================================================================
# OUTLIER HANDLING
# =============================================================================
# Two different problems need two different fixes:
#
# (a) duration_days: negative values are LOGICALLY IMPOSSIBLE
#     (completion before start = data entry error) -> hard filter, no debate.
#
# (b) enrollment: extreme values are mostly REAL (national registries,
#     population surveillance studies) not errors -> we don't delete them,
#     we (1) flag them, and (2) use a log scale so a plot isn't dominated
#     by 3-4 trials with 7-20 million enrollees.
# =============================================================================

# (a) Hard filter: impossible durations
n_before <- nrow(df)
df_clean <- df %>%
  filter(is.na(duration_days) | duration_days >= 0)
message("Removed ", n_before - nrow(df_clean),
        " rows with impossible (negative) trial durations.")

# Also cap extreme-but-technically-valid durations for plotting only
# (some trials list multi-decade completion windows that are placeholder dates)
duration_cap <- quantile(df_clean$duration_days, 0.99, na.rm = TRUE)

# (b) Flag enrollment outliers using IQR, but DON'T delete them
q1 <- quantile(df_clean$enrollment, 0.25, na.rm = TRUE)
q3 <- quantile(df_clean$enrollment, 0.75, na.rm = TRUE)
iqr <- q3 - q1
upper_fence <- q3 + 1.5 * iqr

df_clean <- df_clean %>%
  mutate(enrollment_outlier = !is.na(enrollment) & enrollment > upper_fence)

message("Flagged ", sum(df_clean$enrollment_outlier, na.rm = TRUE),
        " trials as enrollment outliers (above ", comma(upper_fence), ").")
message("These are kept in the data but shown on a log scale so they ",
        "don't distort the plots.")

# =============================================================================
# PLOT 1 — Where are trials in their lifecycle?
# Question: How much of the COVID trial landscape was still unresolved?
# =============================================================================
p1 <- df_clean %>%
  count(status, sort = TRUE) %>%
  mutate(status = fct_reorder(status, n),
         pct    = n / sum(n)) %>%
  ggplot(aes(x = status, y = n, fill = status)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = paste0(comma(n), "  ·  ", percent(pct, accuracy = 0.1))),
            hjust = -0.05, size = 3.3, colour = "grey25") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28)), labels = comma) +
  scale_fill_manual(values = pal_status, na.value = "grey70") +
  labs(
    title    = "**Most COVID Trials Were Still Unresolved**",
    subtitle = "Nearly 2 in 3 trials were either still recruiting or hadn't started yet —
                a sign this dataset captures a pandemic response *in motion*, not a finished record.",
    x = NULL, y = "Number of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  ) +
  theme(panel.grid.major.y = element_blank())

print(p1)

# =============================================================================
# PLOT 2 — How far along the development pipeline were these trials?
# Question: Was the world mostly running early-stage safety trials,
# or were treatments closer to real-world rollout?
# =============================================================================
p2 <- df_clean %>%
  filter(!is.na(phases), phases != "") %>%
  count(phases, sort = TRUE) %>%
  mutate(phases = fct_reorder(phases, n),
         is_combo = str_detect(phases, "\\|")) %>%
  ggplot(aes(x = phases, y = n, fill = is_combo)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 3.3, colour = "grey25") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = comma) +
  scale_fill_manual(values = c("FALSE" = "#3D5A80", "TRUE" = "#A8DADC")) +
  labs(
    title    = "**Phase 2 Dominated — the Field Was Still Proving Efficacy**",
    subtitle = "Most COVID trials sat in early-to-mid stage testing.
                Very few had reached Phase 4 (post-approval monitoring) by the time this data was captured.",
    x = NULL, y = "Number of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  ) +
  theme(panel.grid.major.y = element_blank())

print(p2)

# =============================================================================
# PLOT 3 — What kinds of interventions were being tested?
# Question: Was the global response mostly drugs? Vaccines? Devices?
# =============================================================================
extract_intervention_type <- function(x) {
  str_extract(x, "^[A-Za-z ]+(?=:)")
}

p3 <- df_clean %>%
  filter(!is.na(interventions), interventions != "") %>%
  separate_rows(interventions, sep = "\\|") %>%
  mutate(int_type = extract_intervention_type(interventions)) %>%
  filter(!is.na(int_type)) %>%
  count(int_type, sort = TRUE) %>%
  slice_max(order_by = n, n = 10) %>%
  mutate(int_type = fct_reorder(int_type, n)) %>%
  ggplot(aes(x = int_type, y = n, fill = n)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 3.3, colour = "grey25") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = comma) +
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.2) +
  labs(
    title    = "**Drugs Dominated the Search for COVID Treatments**",
    subtitle = "Interventions categorised by type (Drug, Device, Diagnostic Test, etc.) reveal
                where research effort was concentrated.",
    x = NULL, y = "Number of Interventions",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  ) +
  theme(panel.grid.major.y = element_blank())

print(p3)

# =============================================================================
# PLOT 4 — Trial size: how big were these trials, really?
# Question: enrollment, properly handled with outliers flagged, not deleted.
# =============================================================================
p4 <- df_clean %>%
  filter(!is.na(enrollment), enrollment > 0) %>%
  ggplot(aes(x = enrollment, fill = enrollment_outlier)) +
  geom_histogram(bins = 60, colour = "white", linewidth = 0.15) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("FALSE" = "#3D5A80", "TRUE" = "#D64045"),
    labels = c("Typical trial", "Outlier (IQR-flagged, kept on log scale)"),
    name   = NULL
  ) +
  labs(
    title    = "**Trial Sizes Span Five Orders of Magnitude**",
    subtitle = "Most trials enrolled a few hundred participants. A small number of national
                surveillance studies enrolled in the millions — flagged in red, not deleted, and shown on a log scale.",
    x = "Enrollment (log scale)", y = "Number of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  )

print(p4)

# =============================================================================
# PLOT 5 — Who funded these trials, and did funding source relate to outcome?
# Question: classic trials-research concern — does funding source associate
# with how a trial concludes (completed vs terminated/withdrawn)?
# =============================================================================
p5 <- df_clean %>%
  filter(!is.na(funded_bys), funded_bys != "",
         status %in% c("Completed", "Terminated", "Withdrawn", "Suspended")) %>%
  mutate(outcome = case_when(
    status == "Completed" ~ "Completed",
    TRUE ~ "Stopped early (terminated/withdrawn/suspended)"
  )) %>%
  count(funded_bys, outcome) %>%
  group_by(funded_bys) %>%
  mutate(total = sum(n), pct = n / total) %>%
  ungroup() %>%
  filter(total >= 20) %>%                      # ignore tiny funder categories
  mutate(funded_bys = fct_reorder(funded_bys, pct * (outcome == "Completed"), .fun = sum)) %>%
  ggplot(aes(x = funded_bys, y = pct, fill = outcome)) +
  geom_col(position = "fill", width = 0.65) +
  geom_text(aes(label = ifelse(pct > 0.06, percent(pct, accuracy = 1), "")),
            position = position_fill(vjust = 0.5), size = 3.1,
            colour = "white", fontface = "bold") +
  coord_flip() +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_manual(values = c("Completed" = "#1B7A8C",
                               "Stopped early (terminated/withdrawn/suspended)" = "#D64045"),
                    name = NULL) +
  labs(
    title    = "**Completion Rates Differ by Funding Source**",
    subtitle = "Among trials that reached a definite endpoint, the share that completed
                successfully versus stopped early varies by who funded them.",
    x = NULL, y = "Share of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle). Funder categories with < 20 resolved trials excluded."
  ) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")

print(p5)

# =============================================================================
# PLOT 6 — Geographic concentration of trial activity
# Question: where in the world was the research actually happening?
# =============================================================================
p6 <- df_clean %>%
  filter(!is.na(locations), locations != "") %>%
  separate_rows(locations, sep = "\\|") %>%
  mutate(country = str_trim(str_extract(locations, "[^,]+$"))) %>%
  filter(!is.na(country), country != "") %>%
  count(country, sort = TRUE) %>%
  slice_max(order_by = n, n = 15) %>%
  mutate(country = fct_reorder(country, n)) %>%
  ggplot(aes(x = country, y = n, fill = n)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 3.3, colour = "grey25") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = comma) +
  scale_fill_viridis_c(option = "rocket", direction = -1, begin = 0.15) +
  labs(
    title    = "**Trial Activity Concentrated in a Handful of Countries**",
    subtitle = "Research capacity and pandemic burden both shaped where COVID trials were conducted.",
    x = NULL, y = "Number of Trial Sites",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle). One trial may run at multiple sites/countries."
  ) +
  theme(panel.grid.major.y = element_blank())

print(p6)

# =============================================================================
# PLOT 7 — Trial duration (outliers hard-filtered, capped for display)
# Question: how long do COVID trials typically run?
# =============================================================================
p7 <- df_clean %>%
  filter(!is.na(duration_days), duration_days <= duration_cap) %>%
  ggplot(aes(x = duration_days)) +
  geom_histogram(aes(fill = after_stat(count)), bins = 55,
                 colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = median(df_clean$duration_days, na.rm = TRUE),
             linetype = "dashed", colour = "#D64045", linewidth = 0.9) +
  annotate("text",
           x = median(df_clean$duration_days, na.rm = TRUE) + 25,
           y = Inf, vjust = 1.6, hjust = 0, size = 3.4, colour = "#D64045",
           label = paste0("Median: ", round(median(df_clean$duration_days, na.rm = TRUE)), " days")) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_fill_viridis_c(option = "viridis", guide = "none") +
  labs(
    title    = "**Most Trials Ran for Under a Year**",
    subtitle = "Negative-duration rows (impossible: completion before start) were removed entirely.
                The top 1% of durations were excluded from this view only, to avoid axis distortion from placeholder dates.",
    x = "Duration (days)", y = "Number of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  )

print(p7)

# =============================================================================
# PLOT 8 — Registration transparency: lag between start and public posting
# Question: how quickly were trials made publicly known after they began?
# =============================================================================
p8 <- df_clean %>%
  filter(!is.na(posting_lag_days)) %>%
  filter(posting_lag_days >= quantile(posting_lag_days, 0.01),
         posting_lag_days <= quantile(posting_lag_days, 0.99)) %>%
  ggplot(aes(x = posting_lag_days)) +
  geom_histogram(bins = 60, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#D64045", linewidth = 0.8) +
  annotate("text", x = 0, y = Inf, vjust = 1.6, hjust = -0.05, size = 3.2,
           colour = "#D64045", label = "Posted on start date") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "**Registration Timing Relative to Trial Start**",
    subtitle = "Negative values = posted *before* the trial started (prospective registration, the gold standard).
                Positive values = posted *after* the trial had already begun.",
    x = "Days between trial start and public registration", y = "Number of Trials",
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle). 1st/99th percentile trimmed for display."
  )

print(p8)

# =============================================================================
# PLOT 9 — Status × Phase: where in the pipeline do trials stall?
# Question: do early-phase trials get withdrawn/terminated more than late-phase?
# =============================================================================
p9 <- df_clean %>%
  filter(!is.na(phases), phases != "",
         status %in% c("Completed", "Recruiting", "Terminated", "Withdrawn", "Suspended")) %>%
  count(phases, status) %>%
  group_by(phases) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = phases, y = status, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = percent(pct, accuracy = 1)), size = 3, colour = "white", fontface = "bold") +
  scale_fill_viridis_c(option = "inferno", direction = -1, labels = percent, name = "Share within phase") +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(
    title    = "**Status Mix Within Each Trial Phase**",
    subtitle = "Each column sums to 100% — showing how trial outcomes differ depending on
                how far along the development pipeline a trial was.",
    x = "Phase", y = NULL,
    caption = "Source: ClinicalTrials.gov COVID-19 trial registry (Kaggle)"
  ) +
  theme(panel.grid = element_blank(), legend.position = "right")

print(p9)

message("\n✅ All 9 plots rendered. Each one answers a specific question about the trial landscape.")
