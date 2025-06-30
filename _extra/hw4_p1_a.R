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

counties_2022 <- counties_sf %>%
  left_join(water_2022_filtered %>%
              select(geoid, percent_lacking_plumbing_2022 = percent_lacking_plumbing),
            by = c("GEOID" = "geoid"))

counties_2023 <- counties_sf %>%
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

# Drop areas with no neighbors
no_neighbors <- which(card(nb) == 0)
cat("Number of observations with no neighbors (excluded):", length(no_neighbors), "\n")

if (length(no_neighbors) > 0) {
  spatial_data <- spatial_data[-no_neighbors, ]
  nb <- poly2nb(spatial_data, queen = queen_setting, snap = snap_distance)
}

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
# Count clusters for annotation
hotspot_count <- sum(spatial_data$cluster_type == "Hotspot (Worsening)")
coldspot_count <- sum(spatial_data$cluster_type == "Coldspot (Improving)")

# ────────────────────────────────────────────────
# Define function for rounded rectangle polygon (simple rectangle here, radius param currently unused)
rounded_rect <- function(x, y, width, height, radius = 10){
  coords <- matrix(c(
    x,           y,
    x + width,   y,
    x + width,   y + height,
    x,           y + height,
    x,           y
  ), ncol = 2, byrow = TRUE)
  st_polygon(list(coords))
}

# ────────────────────────────────────────────────
# Set tmap mode and create the final plot
tmap_mode("plot")  # ensure static plotting mode

# Calculate bounding box of spatial_data for positioning annotation
bbox <- st_bbox(spatial_data)

# Base position and size for annotation box & text (adjust as needed to move box)
# Adjusted position and size for annotation box & text
x_pos <- bbox["xmin"] + (bbox["xmax"] - bbox["xmin"]) * 0.50  # moved left from 0.80 to 0.65 (smaller fraction)
y_pos <- bbox["ymin"] + (bbox["ymax"] - bbox["ymin"]) * -0.13  # slightly lower vertical start to keep bottom margin
box_width <- (bbox["xmax"] - bbox["xmin"]) * 0.2              # keep same width
box_height <- (bbox["ymax"] - bbox["ymin"]) * 0.18            # increased height (was 0.1)

# Create the annotation box polygon
annot_box <- st_sfc(rounded_rect(x_pos, y_pos, box_width, box_height), crs = st_crs(spatial_data))

# Adjust text point position with some padding inside the bigger box
annot_text_point <- st_sfc(
  st_point(c(x_pos + box_width * 0.45, y_pos + box_height * 0.50)),  # - move the text rel. to the box.
  crs = st_crs(spatial_data)
)


# Build sf objects for box and text
annot_box_sf <- st_sf(
  type = "box",
  label = NA_character_,
  geometry = annot_box
)

annot_text_sf <- st_sf(
  type = "text",
  label = paste0("Hotspots: ", hotspot_count, "\nColdspots: ", coldspot_count),
  geometry = annot_text_point
)

# Combine for easy plotting
annot_sf <- rbind(annot_box_sf, annot_text_sf)

# ────────────────────────────────────────────────
# Compose tmap
g1 <- tm_shape(spatial_data) +
  tm_polygons("cluster_type",
              fill.scale = tm_scale(values = c(
                "Hotspot (Worsening)"   = "#DC143C",
                "Coldspot (Improving)"  = "#1ad3d1",
                "Not Significant"       = "grey95"),
                legend.title = "Plumbing Access Change Clusters")) +
  tm_title("Spatial Clustering of Plumbing Access Changes (2022–2023)", size = 1.2) +

  # Add annotation box (polygon)
  tm_shape(annot_sf %>% filter(type == "box")) +
  tm_polygons(col = "#ffffff", alpha = 0.8, border.col = "#333333", lwd = 1.2, legend.show = FALSE) +

  # Add annotation text
  tm_shape(annot_sf %>% filter(type == "text")) +
  tm_text("label", size = 0.7, just = c("left", "bottom"), col = "black") +

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
