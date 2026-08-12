# =============================================================================
# Figure 6 - Particle-size distribution (PSD)
# RINENG-D-26-06081 (Revision R1)
#
# Panels:
#   (a) Cumulative percentile curves
#   (b) D10, D50 and D90 versus cooling time
#   (c) Span versus cooling time
#   (d) Differential particle-size distributions (%Chan)
#
# Data source:
# data/external/Granulo_laser_30.xlsx
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("readxl", "ggplot2", "tidyr", "dplyr", "patchwork"))


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------




# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------


if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}


OUTDIR <- file.path("results", "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

XLSX <- file.path("data", "external", "Granulo_laser_30.xlsx")

if (!file.exists(XLSX)) {
  stop("File not found: ", XLSX)
}

SHEET <- "Granulo-S30%"


# -----------------------------------------------------------------------------
# 2. General settings
# -----------------------------------------------------------------------------

time_levels <- c(
  "t0",
  "t12",
  "t24",
  "t48"
)

time_labels <- c(
  "t0"  = "0 h",
  "t12" = "12 h",
  "t24" = "24 h",
  "t48" = "48 h"
)


# Colors used for sampling times
pal_time <- c(
  "t0"  = "#1b9e77",
  "t12" = "#d95f02",
  "t24" = "#7570b3",
  "t48" = "#e7298a"
)


# Colors used for D10, D50 and D90
pal_D <- c(
  "D10" = "#1b9e77",
  "D50" = "#d95f02",
  "D90" = "#7570b3"
)


# Common graphical settings for all four panels
LINE_WIDTH <- 0.55
POINT_SIZE <- 1.5


# -----------------------------------------------------------------------------
# 3. Publication theme
# -----------------------------------------------------------------------------

theme_pub <- function() {
  
  theme_bw(
    base_family = FAM,
    base_size = 8
  ) +
    theme(
      plot.title = element_blank(),
      
      plot.tag = element_text(
        family = FAM,
        face = "bold",
        size = 10,
        colour = "black"
      ),
      
      plot.tag.position = c(0.03, 0.97),
      
      axis.title = element_text(
        family = FAM,
        face = "bold",
        size = 9,
        colour = "black"
      ),
      
      axis.text = element_text(
        family = FAM,
        face = "bold",
        size = 7.5,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.4
      ),
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        family = FAM,
        size = 7.5,
        colour = "black"
      ),
      
      legend.key.size = grid::unit(3.2, "mm"),
      
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.25
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(4, 4, 4, 4)
    )
}


# =============================================================================
# 4. DATA FOR PANELS (a), (b) AND (c)
# =============================================================================

# Column A = percentile
# Columns B:E = t0, t12, t24 and t48

df_pct <- read_excel(
  XLSX,
  sheet = SHEET,
  range = "A1:E11"
)


# Assign explicit column names
names(df_pct) <- c(
  "pct",
  "t0",
  "t12",
  "t24",
  "t48"
)


df_pct_long <- df_pct %>%
  
  pivot_longer(
    cols = all_of(time_levels),
    names_to = "time",
    values_to = "diam_um"
  ) %>%
  
  mutate(
    pct = as.numeric(pct),
    diam_um = as.numeric(diam_um),
    time = factor(
      time,
      levels = time_levels
    )
  ) %>%
  
  filter(
    !is.na(pct),
    !is.na(diam_um)
  )


# =============================================================================
# 5. Calculate D10, D50, D90 and Span
# =============================================================================

df_D <- df_pct_long %>%
  
  filter(
    pct %in% c(
      10,
      50,
      90
    )
  ) %>%
  
  mutate(
    metric = paste0(
      "D",
      pct
    )
  ) %>%
  
  select(
    time,
    metric,
    diam_um
  ) %>%
  
  pivot_wider(
    names_from = metric,
    values_from = diam_um
  ) %>%
  
  mutate(
    Span = (D90 - D10) / D50
  )


df_D_long <- df_D %>%
  
  pivot_longer(
    cols = c(
      D10,
      D50,
      D90
    ),
    names_to = "Dx",
    values_to = "diam_um"
  ) %>%
  
  mutate(
    Dx = factor(
      Dx,
      levels = c(
        "D10",
        "D50",
        "D90"
      )
    )
  )


# Display values in the R console for verification
print(df_D)


# =============================================================================
# 6. DATA FOR PANEL (d)
# =============================================================================
#
# Detailed granulometry block:
#
# B = particle diameter
# C = t0 %Chan
# D = t0 %Pass
# E = t12 %Chan
# F = t12 %Pass
# G = t24 %Chan
# H = t24 %Pass
# I = t48 %Chan
# J = t48 %Pass
#
# Only %Chan is used for panel (d).
# =============================================================================

df_chan_raw <- read_excel(
  XLSX,
  sheet = SHEET,
  range = "B37:J108",
  col_names = FALSE
)


names(df_chan_raw) <- c(
  "size_um",
  "t0_chan",
  "t0_pass",
  "t12_chan",
  "t12_pass",
  "t24_chan",
  "t24_pass",
  "t48_chan",
  "t48_pass"
)


df_chan <- df_chan_raw %>%
  
  select(
    size_um,
    t0_chan,
    t12_chan,
    t24_chan,
    t48_chan
  ) %>%
  
  rename(
    t0  = t0_chan,
    t12 = t12_chan,
    t24 = t24_chan,
    t48 = t48_chan
  ) %>%
  
  pivot_longer(
    cols = all_of(time_levels),
    names_to = "time",
    values_to = "Chan"
  ) %>%
  
  mutate(
    size_um = as.numeric(size_um),
    Chan = as.numeric(Chan),
    time = factor(
      time,
      levels = time_levels
    )
  ) %>%
  
  filter(
    !is.na(size_um),
    !is.na(Chan),
    size_um > 0
  )


# =============================================================================
# 7. PANEL (a)
#    Cumulative percentile curves
# =============================================================================

p_a <- ggplot(
  df_pct_long,
  aes(
    x = diam_um,
    y = pct,
    colour = time,
    group = time
  )
) +
  
  geom_line(
    linewidth = LINE_WIDTH
  ) +
  
  geom_point(
    size = POINT_SIZE,
    shape = 19
  ) +
  
  scale_x_log10(
    breaks = c(
      1,
      10,
      100
    ),
    labels = c(
      "1",
      "10",
      "100"
    )
  ) +
  
  scale_colour_manual(
    values = pal_time,
    labels = time_labels
  ) +
  
  labs(
    tag = "(a)",
    x = "Particle diameter (µm)",
    y = "Percentile (%)"
  ) +
  
  theme_pub() +
  
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.margin = margin(
      0,
      0,
      1,
      0
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = LINE_WIDTH,
        size = POINT_SIZE
      )
    )
  )


# =============================================================================
# 8. PANEL (b)
#    D10, D50 and D90 versus cooling time
# =============================================================================

p_b <- ggplot(
  df_D_long,
  aes(
    x = time,
    y = diam_um,
    colour = Dx,
    group = Dx
  )
) +
  
  geom_line(
    linewidth = LINE_WIDTH
  ) +
  
  geom_point(
    size = POINT_SIZE,
    shape = 19
  ) +
  
  scale_colour_manual(
    values = pal_D,
    labels = expression(
      D[10],
      D[50],
      D[90]
    )
  ) +
  
  scale_x_discrete(
    labels = time_labels
  ) +
  
  labs(
    tag = "(b)",
    x = "Cooling time (h)",
    y = "Diameter (µm)"
  ) +
  
  theme_pub() +
  
  # Explicit formatting guarantees bold axis titles in panel (b)
  theme(
    axis.title.x = element_text(
      family = FAM,
      face = "bold",
      size = 9,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      family = FAM,
      face = "bold",
      size = 9,
      colour = "black"
    ),
    
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.margin = margin(
      0,
      0,
      1,
      0
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = LINE_WIDTH,
        size = POINT_SIZE
      )
    )
  )


# =============================================================================
# 9. PANEL (c)
#    Span versus cooling time
# =============================================================================

p_c <- ggplot(
  df_D,
  aes(
    x = time,
    y = Span,
    group = 1
  )
) +
  
  geom_line(
    linewidth = LINE_WIDTH,
    colour = "steelblue"
  ) +
  
  geom_point(
    size = POINT_SIZE,
    shape = 19,
    colour = "steelblue"
  ) +
  
  scale_x_discrete(
    labels = time_labels
  ) +
  
  labs(
    tag = "(c)",
    x = "Cooling time (h)",
    y = "Span = (D₉₀ − D₁₀) / D₅₀"
  ) +
  
  theme_pub() +
  
  theme(
    legend.position = "none"
  )


# =============================================================================
# 10. PANEL (d)
#     Differential particle-size distribution
# =============================================================================

p_d <- ggplot(
  df_chan,
  aes(
    x = size_um,
    y = Chan,
    colour = time,
    group = time
  )
) +
  
  geom_line(
    linewidth = LINE_WIDTH
  ) +
  
  # Same point size as panels (a), (b) and (c)

  scale_x_log10(
    breaks = c(
      0.02,
      0.1,
      0.5,
      2,
      10,
      50,
      200,
      1000
    ),
    labels = c(
      "0.02",
      "0.1",
      "0.5",
      "2",
      "10",
      "50",
      "200",
      "1000"
    )
  ) +
  
  scale_colour_manual(
    values = pal_time,
    labels = time_labels
  ) +
  
  labs(
    tag = "(d)",
    x = "Particle diameter (µm)",
    y = "Channel fraction (%)"
  ) +
  
  theme_pub() +
  
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.margin = margin(
      0,
      0,
      1,
      0
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = LINE_WIDTH,
        size = POINT_SIZE
      )
    )
  )


# =============================================================================
# 11. Final 2 x 2 assembly
# =============================================================================

p_all <- (
  p_a | p_b
) / (
  p_c | p_d
) +
  
  plot_layout(
    widths = c(
      1,
      1
    ),
    heights = c(
      1,
      1
    )
  )


print(p_all)


# =============================================================================
# 12. Export
# =============================================================================

W <- 14.5
H <- 10.5
UNITS <- "cm"
DPI <- 1200


tif_path <- file.path(
  OUTDIR,
  "Figure_6.tiff"
)

pdf_path <- file.path(
  OUTDIR,
  "Figure_6.pdf"
)


ggsave(
  filename = tif_path,
  plot = p_all,
  width = W,
  height = H,
  units = UNITS,
  dpi = DPI,
  compression = "lzw",
  bg = "white"
)


ggsave(
  filename = pdf_path,
  plot = p_all,
  device = cairo_pdf,
  width = W,
  height = H,
  units = UNITS,
  bg = "white"
)


message(
  "Figure 6 saved in: ",
  OUTDIR
)

message(
  "TIFF: ",
  tif_path
)

message(
  "PDF: ",
  pdf_path
)

# =============================================================================
