# BINF 6210 – Assignment 1
# Bombus BIN Richness vs Latitude (Exploratory)
# Author: Anita Jafari
# Date: 2025-10-05
# Goal: assess whether Bombus BIN richness varies with latitude (N. Hemisphere)

# load the libraries needed
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
})
library(iNEXT)

# set the seed so that it's reproducible
set.seed(203)
#set the theme
theme_set(theme_light(base_size = 12))


####1 - OVERVIEW----
# In this code I will:
# 1) read a BOLD export that I did for Bombus, 
# 2) key fields are then standardized, 
# 3) QC + filter are done (lat > 0),
# 4) make a quick histogram ( this is done for QC), 
# 5) map the points (I'm using sf library),
# 6) I will band latitudes and compute BIN richness,
# 7) optional iNEXT coverage,
# 8) write methods text and session info.

# Prepare output folders (this will be safe if it's being run again (it doesn't overwrite))
dir.create("figs",    showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

####2 - INPUTS----
# Path to BOLD dataset export.
fp <- "../data/result.tsv"
stopifnot(file.exists(fp))  # fail fast if path is wrong

####3 - IMPORT & STANDARDIZATION----
## WHY am I doing this? BOLD headers vary a bit; I'm harmonizing the essentials.

# Read Bombus TSV dataset from BOLD
df_raw <- read_tsv(fp, show_col_types = FALSE)

names(df_raw)

# Standardize key fields
df <- df_raw %>%
  mutate(
    bin          = bin_uri,          # BIN identifier
    country      = `country/ocean`,  # country name 
    processid    = processid,        # process ID
    species_name = species,          # species
    coord_raw    = coord             # coordinate string
  )

# Parse coordinates
extract_num <- function(x, which = 1L) {
  m <- str_match(x, "\\[?\\s*([+-]?[0-9]*\\.?[0-9]+)\\s*,\\s*([+-]?[0-9]*\\.?[0-9]+)\\s*\\]?")
  suppressWarnings(as.numeric(m[, which + 1L]))  # 1 = lat, 2 = lon
}

df <- df %>%
  mutate(
    lat = extract_num(coord_raw, 1L),
    lon = extract_num(coord_raw, 2L)
  )


# Quick snapshot: effort by country (pre-filter)
df %>%
  count(country, sort = TRUE) %>%
  write_csv("outputs/records_by_country.csv")

####4 - FILTERS & QC----
## WHY am I doing this? to keep usable, geo-referenced Northern Hemisphere records with BINs.

df_clean <- df %>%
  filter(!is.na(lat), !is.na(lon)) %>%                       # need coordinates
  filter(lat >= -90, lat <= 90, lon >= -180, lon <= 180) %>% # valid ranges
  filter(!(lat == 0 & lon == 0)) %>%                         # drop (0,0)
  filter(lat > 0) %>%                                        # N. Hemisphere focus
  filter(!is.na(bin) & bin != "") %>%                        # require BIN
  select(processid, bin, species_name, country, lat, lon)    # slim set for downstream

# checking the df_clean
if (nrow(df_clean) == 0) stop("No records after filtering.")

# Save a copy for reproducibility and side analyses if needed further
write_csv(df_clean, "outputs/Bombus_clean_subset.csv")

####5 - FIGURE 1 (QC HISTOGRAM)----
p_qc <- ggplot(df_clean, aes(x = lat)) +
  geom_histogram(bins = 40) +
  labs(
    title    = "Latitude distribution of Bombus records (QC)",
    subtitle = "Northern Hemisphere records with valid coordinates",
    x = "Latitude (°N)", y = "Record count"
  )
ggsave("figs/fig1_QC_latitude_hist.png", p_qc, width = 8, height = 5, dpi = 300)  # write PNG

####6 - FIGURE 2 (MAP)----
## WHY am I doing this? to get a quick view of spatial extent; This is helpful for spotting obvious artifacts.
od_sf <- st_as_sf(df_clean, coords = c("lon","lat"), crs = 4326, remove = FALSE)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p_map <- ggplot() +
  geom_sf(data = world, linewidth = 0.2, fill = "grey95") +
  geom_sf(data = od_sf, alpha = 0.5, size = 0.6) +
  coord_sf(xlim = c(-170, 50), ylim = c(0, 80)) +
  labs(
    title    = "Bombus barcode records used for analysis",
    subtitle = "Northern Hemisphere subset from BOLD",
    x = "Longitude", y = "Latitude",
    caption = "Points = individual georeferenced specimens"
  )
ggsave("figs/fig2_world_map_records.png", p_map, width = 10, height = 6, dpi = 300)

####7 - MAIN: BIN RICHNESS vs LATITUDE----
## WHY am I doing this? to get a simple banding approach (5°) and show richness pattern with latitude.

lat_band_width <- 5  # 5° bands

lat_min <- floor(min(df_clean$lat, na.rm = TRUE))
lat_max <- ceiling(max(df_clean$lat, na.rm = TRUE))

breaks <- seq(lat_min, lat_max, by = lat_band_width)
if (tail(breaks, 1) < lat_max) breaks <- c(breaks, lat_max)  # ensure last break covers max

# Summaries by latitude band
df_rich <- df_clean %>%
  mutate(
    lat_band = cut(lat, breaks = breaks, include.lowest = TRUE, right = FALSE),
    lat_mid  = (as.numeric(lat_band) - 0.5) * lat_band_width + min(breaks)  # midpoint helper
  ) %>%
  filter(!is.na(lat_band)) %>%
  group_by(lat_band, lat_mid) %>%
  summarise(
    BIN_richness = n_distinct(bin),   # unique BINs per band
    records      = n(),               # record count per band (context)
    .groups = "drop"
  ) %>%
  arrange(lat_mid)

write_csv(df_rich, "outputs/Bombus_richness_by_latband.csv")

# Figure 3: Richness vs latitude (midpoint)
p_rich <- ggplot(df_rich, aes(x = lat_mid, y = BIN_richness)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title    = "Bombus BIN richness across latitude",
    subtitle = paste0("Latitude bands = ", lat_band_width, "°; Northern Hemisphere (BOLD)"),
    x = "Latitude band midpoint (°N)",
    y = "BIN richness (unique BINs per band)"
  )
ggsave("figs/fig3_richness_by_latitude.png", p_rich, width = 8, height = 5.5, dpi = 300)

####8 - OPTIONAL: iNEXT COVERAGE----
## WHY doing this? quick sample coverage snapshot by a few data-rich countries.

candidate_countries <- c("Canada","United States","Mexico","Norway","Sweden","Finland")

top_countries <- df_clean %>%
  count(country, sort = TRUE) %>%
  mutate(is_candidate = country %in% candidate_countries) %>%
  arrange(desc(is_candidate), desc(n)) %>%
  slice_head(n = 6) %>%
  pull(country)

# Incidence (presence/absence) table: country x BIN
incidence_table <- df_clean %>%
  filter(country %in% top_countries) %>%
  count(country, bin) %>%
  mutate(inc = 1L) %>%                                        # treat any presence as 1
  select(country, bin, inc) %>%
  pivot_wider(names_from = bin, values_from = inc, values_fill = 0)

# Convert rows to list of incidence vectors (named by country)
inc_list <- apply(incidence_table[,-1, drop = FALSE], 1, as.integer) %>% as.list()
names(inc_list) <- incidence_table$country

# DataInfo from iNEXT (incidence_freq)
covg <- DataInfo(inc_list, datatype = "incidence_freq")
write_csv(covg, "outputs/coverage_estimates.csv")

####10 - SESSION INFO----
capture.output(sessionInfo(), file = "outputs/sessionInfo.txt")
message("\nAll is done! Figures are saved inisde → ./figs, tables are put in → ./outputs.\n")


#Edition

# Create database to filter out: inside vs outside 40–60°
df_rich_t <- df_rich %>%
  mutate(midlat_bin = ifelse(lat_mid >= 40 & lat_mid <= 60, "40-60", "other"))

# t-test 
t.test(BIN_richness ~ midlat_bin, data = df_rich_t)

gam <- ggplot(df_rich, aes(x = lat_mid, y = BIN_richness)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE) +
  labs(x = "Latitude", y = "BIN richness")

p_rich_interactive <- ggplotly(gam)
p_rich_interactive




#Add a VENN Diagram 
#from df_clean, filter for both latitude band 30-45 and 45-60, which contain the most richness in the earlier plots. 
#select only the unique values 
#each bin dataframe should contain two columns: bin and latitudes

bin_30_45 <- unique(df_clean$bin[df_clean$lat >= 30 & df_clean$lat < 45])

bin_45_60 <- unique(df_clean$bin[df_clean$lat >= 45 & df_clean$lat < 60])


#Prepare list of BINs per latitude band from the precious data frame 
bin_list <- list(
  `45°N` = bin_30_45,
  `60°N` = bin_45_60
)

#Plot Venn diagram using ggVennDiagram package, and adjust the marginal size:
ggVennDiagram(bin_list,
              label_color = "white",
              label_alpha = 0, 
              category.names = c("45°N", "60°N")) +
              scale_fill_gradient(low = viridis(1, option = "E"), 
                                 high = viridis(1, option = "D")) +
  labs(
    title = "BIN Overlap Between Latitude Bands",
    subtitle = "Bombus BINs shared and unique in 30–45°N vs 45–60°N"
  ) +
  theme(legend.position = "",
        plot.title = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(2, 2, 2,2), "cm"), 
        plot.subtitle = element_text(size = 14)
)





