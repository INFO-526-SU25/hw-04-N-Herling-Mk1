# ─────────────────────────────────────────────────────
# Load Libraries
library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(tmap)
library(here)
library(viridis)

# ─────────────────────────────────────────────────────
# Global Options
options(tigris_use_cache = TRUE)

# ─────────────────────────────────────────────────────
# Hyperparameters
queen_setting <- TRUE
snap_distance <- 30000
z_threshold   <- 1.96

# ─────────────────────────────────────────────────────
# Load and Clean Data
raw_2022 <- read.csv(here("data", "water_insecurity_2022.csv"))
raw_2023 <- read.csv(here("data", "water_insecurity_2023.csv"))

n_missing_2022 <- sum(is.na(raw_2022$percent_lacking_plumbing))
n_missing_2023 <- sum(is.na(raw_2023$percent_lacking_plumbing))

cat("Excluded due to missing data (2022): n =", n_missing_2022, "\n")
cat("Excluded due to missing data (2023): n =", n_missing_2023, "\n")

water_2022 <- raw_2022 %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

water_2023 <- raw_2023 %>%
  filter(!is.na(percent_lacking_plumbing)) %>%
  mutate(geoid = as.character(geoid))

common_geoids <- intersect(water_2022$geoid, water_2023$geoid)

water_2022 <- water_2022 %>% filter(geoid %in% common_geoids)
water_2023 <- water_2023 %>% filter(geoid %in% common_geoids)

# ─────────────────────────────────────────────────────
# Load US Counties (Continental Only)
continental_states <- c(
  "01","04","05","06","08","09","10","11","12","13","16","17","18","19","20","21",
  "22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37",
  "38","39","40","41","42","44","45","46","47","48","49","50","51","53","54","55","56"
)

counties_sf <- counties(cb = TRUE, year = 2022, class = "sf") %>%
  filter(STATEFP %in% continental_states) %>%
  mutate(GEOID = as.character(GEOID))

# - simplify spatial geometry
counties_sf <- st_simplify(counties_sf, dTolerance = 1000)

# ─────────────────────────────────────────────────────
# Join Water Data to County Geometry
counties_2022 <- counties_sf %>%
  left_join(water_2022 %>%
              select(geoid, percent_lacking_plumbing_2022 = percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

counties_2023 <- counties_sf %>%
  left_join(water_2023 %>%
              select(geoid, percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

# ─────────────────────────────────────────────────────
# Calculate Change
spatial_data <- counties_2023 %>%
  left_join(
    counties_2022 %>% st_drop_geometry() %>%
      select(GEOID, percent_lacking_plumbing_2022),
    by = "GEOID"
  ) %>%
  mutate(change = percent_lacking_plumbing - percent_lacking_plumbing_2022) %>%
  filter(!is.na(change)) %>%
  st_transform(crs = 5070)

# ─────────────────────────────────────────────────────
# Spatial Weights / Neighborhood Structure
nb <- poly2nb(spatial_data, queen = queen_setting, snap = snap_distance)
no_neighbors <- which(card(nb) == 0)

cat("Excluded due to no spatial neighbors: n =", length(no_neighbors), "\n")

if (length(no_neighbors) > 0) {
  spatial_data <- spatial_data[-no_neighbors, ]
  nb <- poly2nb(spatial_data, queen = queen_setting, snap = snap_distance)
}

lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# ─────────────────────────────────────────────────────
# Compute Local G* Statistic
spatial_data$gstar_z <- as.numeric(localG(spatial_data$change, listw = lw, zero.policy = TRUE))

spatial_data <- spatial_data %>%
  mutate(cluster_type = case_when(
    gstar_z >=  z_threshold  ~ "Hotspot (Worsening)",
    gstar_z <= -z_threshold  ~ "Coldspot (Improving)",
    TRUE                     ~ "Not Significant"
  ))

# ─────────────────────────────────────────────────────
# Identify and Label Excluded Counties
# Counties excluded due to missing data (2022 or 2023)
excluded_missing_geoids <- setdiff(
  counties_sf$GEOID,
  union(water_2022$geoid, water_2023$geoid)
)

excluded_missing <- counties_sf %>%
  filter(GEOID %in% excluded_missing_geoids) %>%
  mutate(cluster_type = "Excluded: Missing Data")

# After removing missing, we filtered for no-neighbor exclusions
# The ones dropped after nb creation are "No Neighbors"

excluded_neighbors <- counties_sf %>%
  filter(GEOID %in% setdiff(counties_sf$GEOID, spatial_data$GEOID),
         !GEOID %in% excluded_missing_geoids) %>%
  mutate(cluster_type = "Excluded: No Neighbors")
# ─────────────────────────────────────────────────────
excluded_missing <- st_transform(excluded_missing, crs = 5070)
excluded_neighbors <- st_transform(excluded_neighbors, crs = 5070)


# Combine Included and Excluded
# Combine Included and Excluded
all_counties <- bind_rows(
  spatial_data %>% select(GEOID, cluster_type, geometry),
  excluded_missing %>% select(GEOID, cluster_type, geometry),
  excluded_neighbors %>% select(GEOID, cluster_type, geometry)
)

# Count Cluster Types
counts <- all_counties %>%
  st_drop_geometry() %>%
  count(cluster_type)

label_vec <- counts %>%
  mutate(label = paste0(cluster_type, " (n=", n, ")")) %>%
  select(cluster_type, label) %>%
  deframe()

# Define cluster_type factor levels explicitly
all_counties$cluster_type <- factor(all_counties$cluster_type, levels = counts$cluster_type)

# Create cluster_label with counts for legend
all_counties$cluster_label <- label_vec[all_counties$cluster_type]

# Filter out any rows where cluster_label is NA (avoid 'Missing' legend)
all_counties <- all_counties %>%
  filter(!is.na(cluster_label))

# Reset cluster_label factor levels to only those present after filtering
all_counties$cluster_label <- factor(all_counties$cluster_label, levels = unique(all_counties$cluster_label))

# Define original colors
cluster_colors_orig <- c(
  "Hotspot (Worsening)"      = "#440154FF",
  "Coldspot (Improving)"     = "#35B779FF",
  "Not Significant"          = "#FDE725FF",
  "Excluded: Missing Data"   = "#CCCFFF",
  "Excluded: No Neighbors"   = "#aaaaaa"
)

# Subset colors and assign names for legend matching filtered data
cluster_colors <- cluster_colors_orig[counts$cluster_type]
names(cluster_colors) <- label_vec[counts$cluster_type]

# Check for any missing colors for cluster_labels
missing_colors <- setdiff(levels(all_counties$cluster_label), names(cluster_colors))
if (length(missing_colors) > 0) {
  stop("Missing colors for these cluster_label levels: ", paste(missing_colors, collapse = ", "))
}

# =====================================
# Count Cluster Types
counts <- all_counties %>%
  st_drop_geometry() %>%
  count(cluster_type)

label_vec <- counts %>%
  mutate(label = paste0(cluster_type, " (n=", n, ")")) %>%
  select(cluster_type, label) %>%
  deframe()

all_counties$cluster_type <- factor(all_counties$cluster_type, levels = counts$cluster_type)
all_counties$cluster_label <- label_vec[all_counties$cluster_type]

# ─────────────────────────────────────────────────────
# Set Colors for Legend
cluster_colors_orig <- c(
  "Hotspot (Worsening)"      = "#440154FF",
  "Coldspot (Improving)"     = "#35B779FF",
  "Not Significant"          = "#FDE725FF",
  "Excluded: Missing Data"   = "#CCCFFF",  # gray or white
  "Excluded: No Neighbors"   = "#aaaaaa"   # slightly different gray for distinction
)


cluster_colors <- cluster_colors_orig[counts$cluster_type]
names(cluster_colors) <- label_vec[counts$cluster_type]

# ─────────────────────────────────────────────────────
# Print Final Summary
cat("\n──────────────────────────────\n")
cat("Summary of Exclusions:\n")
cat("- Missing data (2022): n =", n_missing_2022, "\n")
cat("- Missing data (2023): n =", n_missing_2023, "\n")
cat("- Excluded counties due to missing data: n =", nrow(excluded_missing), "\n")
cat("- No spatial neighbors: n =", length(no_neighbors), "\n")
cat("- Total excluded counties: n =", nrow(excluded_counties), "\n")
cat("Total counties in map: n =", nrow(all_counties), "\n")
cat("──────────────────────────────\n")


# ─────────────────────────────────────────────────────
# Exclude any "Missing..." labels from the color palette and legend
exclude_label <- "Missing..."

# If 'Missing...' actually appears in your factor levels:
if (exclude_label %in% levels(all_counties$cluster_label)) {

  keep_labels <- levels(all_counties$cluster_label)[levels(all_counties$cluster_label) != exclude_label]

  cluster_colors <- cluster_colors[names(cluster_colors) %in% keep_labels]

  all_counties$cluster_label <- factor(all_counties$cluster_label, levels = keep_labels)
}


# Plot the Map
tmap_mode("plot")
tmap_options(component.autoscale = FALSE)


# Prepare your factor and color palette
all_counties$cluster_label <- droplevels(factor(all_counties$cluster_label))
cluster_colors <- cluster_colors[names(cluster_colors) %in% levels(all_counties$cluster_label)]

# Build the map with NO default legend
g1 <- tm_shape(all_counties) +
  tm_polygons(
    col = "cluster_label",
    palette = cluster_colors,
    alpha = 0.8,
    border.col = "black",
    lwd = 0.4,
    title = NULL,
    showNA = FALSE
  ) +
  tm_layout(
    legend.show = FALSE,
    inner.margins = c(0.2, 0.02, 0.12, 0.02)
  ) +
  tm_add_legend(
    type = "fill",
    labels = names(cluster_colors),
    col = unname(cluster_colors),
    title = "Plumbing Access Change Clusters",
    size = 0.7,
    position = "left"  # fixed here
  ) +
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
  )

plot(g1)