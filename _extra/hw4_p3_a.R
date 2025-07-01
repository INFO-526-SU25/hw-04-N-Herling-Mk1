# ----------------------------------------
# Load necessary libraries
# ----------------------------------------
library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)
library(stringr)  # for string trimming

# ----------------------------------------
# Set tigris options and version
# ----------------------------------------
options(tigris_use_cache = TRUE)
tigris_version <- as.character(packageVersion("tigris"))

# ----------------------------------------
# Step 1: Load and prepare Arizona counties shapefile
# ----------------------------------------
az_counties <- counties(state = "AZ", year = 2021, progress_bar = FALSE) %>%
  st_as_sf() %>%
  mutate(
    NAME = gsub(" County$", "", NAME),
    NAME = str_trim(NAME)  # Trim whitespace for safety
  )

# ----------------------------------------
# Step 2: Load population data with proper skipping and renaming
# ----------------------------------------
pop_data_raw <- suppressMessages(suppressWarnings(
  read_excel(
    "data/co-est2023-pop-04.xlsx",
    skip = 4,
    n_max = 15,
    col_names = FALSE
  )
))

# Assign column names manually based on inspected data:
# Column 1: county (with leading dot and suffix)
# Columns 2-6: some codes or estimates; we want 2021 and 2022 estimates at columns 5 and 6 after skipping
colnames(pop_data_raw) <- c(
  "county_raw", "code1", "code2", "est_2020", "est_2021", "est_2022"
)

# Clean and compute population change
pop_data <- pop_data_raw %>%
  mutate(
    NAME = county_raw %>%
      str_remove_all("^\\.") %>%                    # Remove leading dot
      str_remove_all(" County, Arizona$") %>%       # Remove suffix
      str_trim(),                                   # Trim spaces
    total_pop_change_21_22 = est_2022 - est_2021
  ) %>%
  select(NAME, total_pop_change_21_22)

# ----------------------------------------
# Step 3: Check for mismatches in county names
# ----------------------------------------
cat("Counties in shapefile but NOT in pop_data:\n")
print(setdiff(az_counties$NAME, pop_data$NAME))

cat("\nCounties in pop_data but NOT in shapefile:\n")
print(setdiff(pop_data$NAME, az_counties$NAME))

# ----------------------------------------
# Step 4: Merge population data with spatial data
# ----------------------------------------
az_counties_merged <- az_counties %>%
  left_join(pop_data, by = "NAME")

# ----------------------------------------
# Step 5: Print counties and population changes to terminal
# ----------------------------------------
for(i in seq_len(nrow(az_counties_merged))) {
  cat(
    sprintf(
      "county name: '%-15s', population change: '%s'\n",
      az_counties_merged$NAME[i],
      ifelse(is.na(az_counties_merged$total_pop_change_21_22[i]), "NA", az_counties_merged$total_pop_change_21_22[i])
    )
  )
}

# ----------------------------------------
# Step 6: Plot population change map with diverging RdBu palette and legend
# ----------------------------------------
library(RColorBrewer)
rdbu_colors <- brewer.pal(11, "RdBu")

g3 <- ggplot(az_counties_merged) +
  geom_sf(aes(fill = total_pop_change_21_22), color = "black") +
  scale_fill_gradient2(
    low = rdbu_colors[1],    # Blue (negative)
    mid = rdbu_colors[6],    # White (zero)
    high = rdbu_colors[11],  # Red (positive)
    midpoint = 0,
    name = "Population Change\n(2021–2022)",
    labels = scales::comma
  ) +
  labs(
    title = "Population Change in Arizona Counties (2021–2022)",
    x = "Longitude",
    y = "Latitude",
    caption = paste0("Sources: U.S. Census Bureau; shapefile from {tigris} R package (v", tigris_version, ")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.caption = element_text(size = 8, hjust = 1, color = "gray40"),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  ) +
  coord_sf()

print(g3)
