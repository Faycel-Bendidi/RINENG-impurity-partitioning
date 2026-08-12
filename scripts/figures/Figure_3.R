# =============================================================================
# Figure 3 - RINENG-D-26-06081 (Revision R1)
# Refroidissement force et continu de l'acide phosphorique PA30 %
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2"))


# Times New Roman
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

# -----------------------------------------------------------------------------
# Chemins et paramètres
# -----------------------------------------------------------------------------

OUTDIR <- file.path("results", "figures")

W <- 9
H <- 6.5
UNITS <- "cm"
DPI <- 1200

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Données
# -----------------------------------------------------------------------------

d <- data.frame(
  temps = c(0, 12, 24, 36, 48, 60, 72),
  T_C   = c(70.0, 41.1, 28.9, 23.8, 21.6, 20.7, 20.0)
)

d$id <- paste0("t", d$temps)

# Étiquettes des points
d$lab <- sprintf("%s (%d ; %.1f)", d$id, d$temps, d$T_C)

# Décalage des étiquettes par rapport aux points
d$dx <- c(1.5, 1.5, 1.5, 1.5, 2.2, 1.5, 2.0)
d$dy <- c(3.6, 3.6, 3.6, 3.6, 7.6, 3.6, -3.8)
d$hj <- c(0, 0, 0, 0, 0, 0, 1.0)
d$vj <- c(0, 0, 0, 0, 0, 0, 1.0)

# Courbe lissée passant par les sept points
f <- splinefun(
  d$temps,
  d$T_C,
  method = "natural"
)

courbe <- data.frame(
  temps = seq(0, 72, by = 0.25)
)

courbe$T_C <- f(courbe$temps)

cat(
  "Ecart maximal entre la courbe et les points mesures :",
  format(max(abs(f(d$temps) - d$T_C)), digits = 3),
  "degC\n"
)

cat(
  "Courbe strictement decroissante :",
  all(diff(courbe$T_C) <= 1e-9),
  "| minimum :",
  format(min(courbe$T_C), nsmall = 2),
  "degC\n"
)

# -----------------------------------------------------------------------------
# 2. Tracé
# -----------------------------------------------------------------------------

p <- ggplot() +
  geom_line(
    data = courbe,
    aes(temps, T_C),
    colour = "#2E75B6",
    linewidth = 1.1
  ) +
  geom_point(
    data = d,
    aes(temps, T_C),
    colour = "red",
    fill = "red",
    shape = 21,
    size = 3.0,
    stroke = 0.6
  ) +
  geom_text(
    data = d,
    aes(
      x = temps + dx,
      y = T_C + dy,
      label = lab,
      hjust = hj,
      vjust = vj
    ),
    size = 6.4,
    fontface = "bold",
    colour = "black",
    family = FAM
  ) +
  scale_x_continuous(
    breaks = seq(0, 72, by = 6),
    limits = c(-3, 76),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(20, 80, by = 20),
    limits = c(14, 82),
    expand = c(0, 0)
  ) +
  labs(
    x = "Cooling time (hour)",
    y = "Temperature (\u00b0C)"
  ) +
  theme_bw(
    base_size = 18,
    base_family = FAM
  ) +
  theme(
    axis.title = element_text(
      face = "bold",
      size = 9
    ),
    axis.text = element_text(
      face = "bold",
      size = 8,
      colour = "black"
    ),
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "grey35",
      linewidth = 0.7
    ),
    plot.margin = margin(8, 12, 8, 8)
  )

print(p)

# -----------------------------------------------------------------------------
# 3. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------

tif <- file.path(OUTDIR, "Figure_3.tiff")
pdf <- file.path(OUTDIR, "Figure_3.pdf")

ggsave(
  tif,
  p,
  width = W,
  height = H,
  units = UNITS,
  dpi = DPI,
  compression = "lzw"
)

ggsave(
  pdf,
  p,
  width = W,
  height = H,
  units = UNITS,
  device = cairo_pdf,
  family = FAM
)

mo <- file.info(tif)$size / 1024^2

message("Ecriture dans : ", OUTDIR)

message(
  sprintf(
    "  Figure_3.tiff : %6.1f Mo%s",
    mo,
    if (mo > 10) {
      "  <-- depasse la limite Elsevier de 10 Mo"
    } else {
      ""
    }
  )
)

message(
  sprintf(
    "  Figure_3.pdf  : %6.1f Mo",
    file.info(pdf)$size / 1024^2
  )
)
