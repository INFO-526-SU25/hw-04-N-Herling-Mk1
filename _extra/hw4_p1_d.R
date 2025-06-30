# Load libraries
library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(tmap)
library(here)

# Cache tigris shapefiles
options(tigris_use_cache = TRUE)

# ────────────────────────────────────────────────
# Hyperparameters
queen_setting <- TRUE
snap_distance <- 20000
z_threshold <- 1.96

# ────────────────────────────────────────────────
# Load and clean data
water_2022_clean <- read.csv(here("data", "water_insecurity_2022.csv")) %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

water_2023_clean <- read.csv(here("data", "water_insecurity_2023.csv")) %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

# Filter to common GEOIDs
common_geoids <- intersect(water_2022_clean$geoid, water_2023_clean$geoid)
water_2022_filtered <- water_2022_clean %>% filter(geoid %in% common_geoids)
water_2023_filtered <- water_2023_clean %>% filter(geoid %in% common_geoids)

# ────────────────────────────────────────────────
# Join to shapefile
counties_sf <- counties(cb = TRUE, year = 2022, class = "sf")

# Filter to continental US counties (exclude AK '02', HI '15', PR '72', etc.)
continental_states <- setdiff(unique(counties_sf$STATEFP), c("02", "15", "72"))
counties_sf_continental <- counties_sf %>%
  filter(STATEFP %in% continental_states)

total_counties_continental <- nrow(counties_sf_continental)

counties_2022 <- counties_sf_continental %>%
  left_join(water_2022_filtered %>%
              select(geoid, percent_lacking_plumbing_2022 = percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

counties_2023 <- counties_sf_continental %>%
  left_join(water_2023_filtered %>%
              select(geoid, percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

# ────────────────────────────────────────────────
# Calculate change and project
spatial_data <- counties_2023 %>%
  left_join(counties_2022 %>% st_drop_geometry() %>%
              select(GEOID, percent_lacking_plumbing_2022),
            by = "GEOID") %>%
  mutate(change = percent_lacking_plumbing - percent_lacking_plumbing_2022) %>%
  filter(!is.na(change))

# Project to Albers Equal Area
spatial_data <- st_transform(spatial_data, crs = 5070)

# ────────────────────────────────────────────────
# Create neighbors
nb <- poly2nb(spatial_data, queen = queen_setting, snap = snap_distance)

# Identify and exclude counties with no neighbors
no_neighbors <- which(card(nb) == 0)
num_no_neighbors <- length(no_neighbors)

if (num_no_neighbors > 0) {
  spatial_data <- spatial_data[-no_neighbors, ]
  nb <- poly2nb(spatial_data, queen = queen_setting, snap = snap_distance)
}

# Total counties analyzed after exclusions
total_analyzed <- nrow(spatial_data)

# ────────────────────────────────────────────────
# Print alert summary
cat(
  "Summary:\n",
  "Total counties in continental US: ", total_counties_continental, "\n",
  "Counties excluded due to missing neighbors: ", num_no_neighbors, "\n",
  "Total counties analyzed: ", total_analyzed, "\n"
)

# ────────────────────────────────────────────────
# Create spatial weights
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Compute Local G* statistic
local_gstar <- localG(spatial_data$change, listw = lw, zero.policy = TRUE)
spatial_data$gstar_z <- as.numeric(local_gstar)

# ────────────────────────────────────────────────
# Classify clusters based on z-score thresholds
spatial_data <- spatial_data %>%
  mutate(cluster_type = case_when(
    gstar_z >=  z_threshold  ~ "Hotspot (Worsening)",
    gstar_z <= -z_threshold  ~ "Coldspot (Improving)",
    TRUE                     ~ "Not Significant"
  ))

# ────────────────────────────────────────────────
# Count clusters for annotation and update labels
cluster_counts <- spatial_data %>%
  count(cluster_type) %>%
  complete(cluster_type = c("Hotspot (Worsening)", "Coldspot (Improving)", "Not Significant"), fill = list(n = 0))

cluster_labels <- cluster_counts %>%
  mutate(label = paste0(cluster_type, " [n=", n, "]")) %>%
  arrange(factor(cluster_type, levels = c("Hotspot (Worsening)", "Coldspot (Improving)", "Not Significant"))) %>%
  pull(label)

# Update cluster_type to factor with counts in labels
spatial_data$cluster_type <- factor(spatial_data$cluster_type,
                                    levels = c("Hotspot (Worsening)", "Coldspot (Improving)", "Not Significant"),
                                    labels = cluster_labels)

# ────────────────────────────────────────────────
# Set tmap mode and create the final plot
tmap_mode("plot")

cluster_colors <- c(
  "#FDE725FF",
  "#21908CFF",
  "#440154FF"
)
names(cluster_colors) <- cluster_labels

# ────────────────────────────────────────────────
# Get state boundaries for overlay
states_sf <- states(cb = TRUE, year = 2022, class = "sf") %>%
  filter(STATEFP %in% continental_states) %>%
  st_transform(crs = st_crs(spatial_data))  # match CRS with spatial_data

# ────────────────────────────────────────────────
# Add overlay to existing plot
g1 <- tm_shape(spatial_data) +
  tm_polygons("cluster_type",
              fill.scale = tm_scale(values = cluster_colors,
                                    legend.title = "Plumbing Access Change Clusters")) +
  tm_shape(states_sf) +  # overlay begins here
  tm_borders(lwd = 1.2, col = "black") +
  tm_title("Spatial Clustering of Plumbing Access Changes (2022–2023)", size = 1.2) +
  tm_credits(
    paste0(
      "Local G* Statistic (z ≥ ±", z_threshold, ")\n",
      "Contiguity = ", ifelse(queen_setting, "Queen", "Rook"),
      " | Snap Distance = ", snap_distance, " meters"
    ),
    position = c("center", "top"),
    size = 0.7,
    lines = 2
  ) +
  tm_layout(
    legend.position = c("left", "bottom"),
    legend.just = "right",
    legend.outside = FALSE,
    legend.text.size = 0.6,
    legend.title.size = 0.75,
    inner.margins = c(0.2, 0.02, 0.12, 0.02)
  )

print(g1)