# =============================================================================
#  Figure 7 - RINENG-D-26-06081 (Revision R1)
#
#  Une seule image : 10 elements majeurs (5 x 2) au-dessus des 7 elements
#  traces (4 x 2). Mise en forme strictement identique aux Figures 4 et 5.
#
#  Grandeur tracee : D_i = C_sol,i / C_liq,i (Eq. 3), coefficient de partage
#  apparent solide/liquide.
#
#  Reponse aux commentaires du Reviewer 5 :
#    - 14 : lisibilite (grille de 17 panneaux reorganisee, panneaux agrandis)
#    -  9 : barres d'erreur (n = 2, +/- 1 ecart-type)
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("readxl", "ggplot2", "tidyr", "dplyr", "readr", "patchwork"))


# -----------------------------------------------------------------------------
# Chemins
# -----------------------------------------------------------------------------
OUTDIR <- file.path("results", "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
XLSX    <- file.path("data", "external", "Supplementary_Data.xlsx")

if (!file.exists(XLSX))
  stop("Fichier introuvable : ", XLSX,
       "\n  Verifier le chemin relatif depuis la racine du depot.")

# -----------------------------------------------------------------------------
# 1. Lecture du fichier Excel
#    Feuille S2 : 3 blocs empiles (D_i, f, R) ; chaque variable occupe
#    3 colonnes (Test | Mean | SD) et chaque temps occupe 2 lignes, Mean et SD
#    n'etant remplis que sur la premiere des deux.
#    HDR_ROW = 1 : ligne d'en-tete du bloc des coefficients ("DP2O5", ...).
#    La notation D_i retenue ici reprend celle du fichier source.
# -----------------------------------------------------------------------------
TIMES   <- c(0, 12, 24, 36, 48, 60, 72)
HDR_ROW <- 1

raw <- as.data.frame(read_excel(XLSX, sheet = "Table S2_coef_mean_SD",
                                col_names = FALSE))

res <- list()
for (cc in seq(2, ncol(raw), by = 3)) {
  nm <- raw[HDR_ROW, cc]
  if (is.na(nm) || nm == "") next
  rows <- HDR_ROW + 2 + 2 * (seq_along(TIMES) - 1)
  res[[length(res) + 1]] <- data.frame(
    Element = as.character(nm),
    t_hours = TIMES,
    Value   = suppressWarnings(as.numeric(raw[rows, cc + 1])),
    SD      = suppressWarnings(as.numeric(raw[rows, cc + 2])),
    stringsAsFactors = FALSE)
}
df_long <- bind_rows(res)

# Retirer UNIQUEMENT les especes entierement sous la limite de detection
# (DCo et DPb) : leur coefficient de partage est indefini a tous les temps.
df_long <- df_long %>%
  group_by(Element) %>% filter(!all(is.na(Value))) %>% ungroup()

# ---------------------------------------------------------------------------
# 1bis. Valeurs estimees par kNN (feuille S3)
#       Les six valeurs sous la limite de detection dans le solide sont
#       remplacees par leur estimation kNN, celle-la meme qui a servi a
#       l'analyse multivariee. Elles sont marquees d'une croix et n'ont pas de
#       barre d'erreur, n'etant pas issues de mesures repliquees.
#       La feuille S3 nomme les variables "D_K2O" (avec tiret bas) alors que la
#       feuille S2 utilise "DK2O" : le tiret bas est retire pour la jointure.
# ---------------------------------------------------------------------------
s3 <- as.data.frame(read_excel(XLSX, sheet = "Table S3_kNN imputation and v",
                               col_names = FALSE))

ok  <- grepl("^D_", s3[[1]]) & !is.na(suppressWarnings(as.numeric(s3[[4]])))
imp <- data.frame(
  Element = gsub("_", "", s3[[1]][ok]),
  t_hours = as.numeric(s3[[2]][ok]),
  Est     = as.numeric(s3[[4]][ok]),
  stringsAsFactors = FALSE)

message("Valeurs kNN reprises de la feuille S3 : ", nrow(imp))
print(imp)

df_long <- df_long %>%
  left_join(imp, by = c("Element", "t_hours")) %>%
  mutate(imputed = is.na(Value) & !is.na(Est),
         Value   = ifelse(imputed, Est, Value)) %>%
  select(-Est)

if (any(is.na(df_long$Value)))
  warning("Des valeurs restent manquantes apres reprise des estimations kNN.")

# -----------------------------------------------------------------------------
# 2. Familles et titres des panneaux
# -----------------------------------------------------------------------------
major_elements <- c("DAl2O3", "DCaO", "DF", "DFe2O3", "DK2O",
                    "DMgO", "DNa2O", "DP2O5", "DSiO2", "DSO3")
trace_elements <- c("DCd", "DCr", "DCu", "DMn", "DMo", "DNi", "DZn")

lab_map <- c(
  "DP2O5"  = "bold(P[2]*O[5])",
  "DAl2O3" = "bold(Al[2]*O[3])",
  "DCaO"   = "bold(CaO)",
  "DFe2O3" = "bold(Fe[2]*O[3])",
  "DK2O"   = "bold(K[2]*O)",
  "DMgO"   = "bold(MgO)",
  "DNa2O"  = "bold(Na[2]*O)",
  "DSO3"   = "bold(SO[3])",
  "DSiO2"  = "bold(SiO[2])",
  "DF"     = "bold(F)",
  "DCd"    = "bold(Cd)",
  "DCr"    = "bold(Cr)",
  "DCu"    = "bold(Cu)",
  "DMn"    = "bold(Mn)",
  "DMo"    = "bold(Mo)",
  "DNi"    = "bold(Ni)",
  "DZn"    = "bold(Zn)"
)

# Titre de l'axe des y : D avec i en indice, defini comme le rapport des
# concentrations solide/liquide. "sol" et "liq" sont en exposant, comme "liq"
# dans la Figure 5, pour garder une notation coherente entre les trois figures.
YLAB <- expression(bold(italic(D)[italic(i)] ==
                          italic(C)[italic(i)]^{"sol"} / italic(C)[italic(i)]^{"liq"}))

# Ordre alphabetique sur le symbole de l'element, et non sur le nom Excel :
# le prefixe "D" est retire avant le tri.
sym <- function(v) sub("^D", "", v)

# -----------------------------------------------------------------------------
# 3. Fonction de trace, commune aux deux blocs
# -----------------------------------------------------------------------------
# legend_pos : NULL pour ne pas afficher de legende ; sinon un vecteur c(x, y)
# en coordonnees relatives, qui place la legende dans la case vide de la grille.
make_block <- function(species, couleur, ncol_grid, legend_pos = NULL) {
  
  d   <- df_long %>% filter(Element %in% species)
  ord <- species[order(sym(species))]
  d   <- d %>% mutate(Label = factor(lab_map[Element], levels = unname(lab_map[ord])))
  
  ggplot(d, aes(x = t_hours, y = Value, group = Element)) +
    # Barres d'erreur en noir, tracees en premier pour passer sous les points
    geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
                  width = 2.5, linewidth = 0.5, colour = "black", na.rm = TRUE) +
    # Le trace est continu : les valeurs manquantes ont ete remplacees par les
    # estimations kNN de la feuille S3.
    geom_line(linewidth = 0.8, colour = couleur, na.rm = TRUE) +
    # Cercle plein = valeur mesuree ; croix = valeur estimee par kNN. L'ecart
    # type etant NA pour ces dernieres, geom_errorbar les ignore d'office.
    geom_point(aes(shape = imputed, size = imputed),
               colour = couleur, stroke = 0.6, na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 4),
                       labels = c(`FALSE` = "Experimental value",
                                  `TRUE`  = "kNN-estimated value"),
                       name   = NULL) +
    scale_size_manual(values  = c(`FALSE` = 0.75, `TRUE` = 1.8), guide = "none") +
    guides(shape = if (is.null(legend_pos)) "none" else
      guide_legend(override.aes = list(size = 3, stroke = 0.9,
                                       colour = "black"))) +
    facet_wrap(~ Label, scales = "free", ncol = ncol_grid,
               labeller = label_parsed) +
    scale_x_continuous(breaks = c(0, 20, 40, 60)) +
    labs(x = "Time (hour)", y = YLAB, title = NULL) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      
      axis.title = element_text(
        face = "bold",
        size = 20
      ),
      
      # The y-axis title is removed from each block because
      # it is added only once in the final assembly.
      axis.title.y = element_blank(),
      
      axis.text = element_text(
        size = 14,
        colour = "black"
      ),
      
      # Black border around each facet title
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.7
      ),
      
      strip.text = element_text(
        size = 17,
        margin = margin(t = 1, b = 1)
      ),
      
      # Black border around each individual plotting panel
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.7
      ),
      
      panel.spacing = unit(
        2.5,
        "mm"
      ),
      
      plot.margin = margin(
        2,
        4,
        2,
        2
      ),
      
      panel.grid.major = element_line(
        colour = "grey90"
      ),
      
      panel.grid.minor = element_blank()
    ) -> p
  
  # Legende placee dans la case laissee vide par la grille (bloc traces :
  # 7 panneaux sur 8 emplacements). ggplot2 >= 3.5 exige "inside" plutot
  # qu'un vecteur passe directement a legend.position.
  if (!is.null(legend_pos)) {
    theme_leg <- theme(
      legend.justification = c(0.5, 0.5),
      legend.background    = element_rect(fill = "white", colour = "grey40",
                                          linewidth = 0.4),
      legend.text          = element_text(size = 14, colour = "black"),
      legend.key           = element_blank(),
      legend.key.size      = unit(6, "mm"),
      legend.margin        = margin(5, 8, 5, 5)
    )
    if (utils::packageVersion("ggplot2") >= "3.5.0") {
      p <- p + theme_leg +
        theme(legend.position = "inside", legend.position.inside = legend_pos)
    } else {
      p <- p + theme_leg + theme(legend.position = legend_pos)
    }
  }
  p
}

p_major <- make_block(major_elements, "steelblue", 5)   # 10 panneaux : 5 x 2
# La grille des traces compte 8 emplacements pour 7 panneaux : la legende
# occupe la case vide, en bas a droite.
p_trace <- make_block(trace_elements, "red", 4, legend_pos = c(0.88, 0.22))

# -----------------------------------------------------------------------------
# 4. Assemblage en une seule image.
#    Le titre de l'axe des x n'apparait qu'une fois, sous le bloc inferieur, et
#    celui de l'axe des y une seule fois a gauche, centre sur toute la hauteur.
# -----------------------------------------------------------------------------
blocs <- (p_major + theme(axis.title.x = element_blank())) / p_trace +
  plot_layout(heights = c(1, 1))

ylab_grob <- grid::textGrob(YLAB, rot = 90, gp = grid::gpar(fontsize = 20))

p_all <- wrap_elements(ylab_grob) | blocs
p_all <- p_all + plot_layout(widths = c(1, 32))   # largeur reservee au titre

print(p_all)

# -----------------------------------------------------------------------------
# 5. Export : TIFF 1200 dpi + PDF vectoriel, dans le meme dossier
# -----------------------------------------------------------------------------
W <- 14.5 ; H <- 10.5 ; UNITS <- "in" ; DPI <- 1200

tif <- file.path(OUTDIR, "Figure_7.tiff")
pdf <- file.path(OUTDIR, "Figure_7.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS,
       dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_7.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_7.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   geom_point(size = 0.6)                             -> points plus petits
#   geom_errorbar(linewidth = 0.65, width = 3.5)       -> barres plus marquees
#   strip.text = element_text(size = 17, margin = margin(t = 0, b = 0))
#                                                      -> titres plus serres
#   plot_layout(widths = c(1, 26))                     -> titre de l'axe des y
#                                                         plus eloigne des graphes
#   axis.text = element_text(face = "bold", size = 14, colour = "black")
#                                                      -> graduations en gras
#   scale_size_manual(values = c(`FALSE` = 0.75, `TRUE` = 2.4))
#                                                      -> croix plus grandes
#   legend_pos = c(0.88, 0.22)                         -> deplacer la legende
#                                                         dans la case vide
#
# Si le titre de l'axe des y est trop long a l'affichage, la version courte
#   YLAB <- expression(bold(italic(D)[italic(i)]))
# suffit, la definition figurant dans l'Eq. 3 et la legende.
#
# Taille : 14.5 x 10 pouces = 368 x 254 mm. La largeur utile du journal est de
# 190 mm ; si le TIFF depasse 10 Mo, remplacer la ligne des parametres par
#   W <- 190 ; H <- 131 ; UNITS <- "mm"
# ou televerser le PDF, vectoriel et prefere par Elsevier.
# =============================================================================
