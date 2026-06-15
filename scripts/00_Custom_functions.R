# Custom functions to save images and df #

# Save images ####
# Note: Make sure that, for the custom function `save_image`,
# the plot you wish to save is displayed as the active plot in the “Plots” pane.
save_image <- function(
    name,
    plot = last_plot(),
    orientation = c("landscape", "portrait"),
    width = 16,
    height = 9,
    units = "in",
    dpi = 300,
    ...
) {
  orientation <- match.arg(orientation)
  
  # swap width/height if portrait
  if (orientation == "portrait") {
    tmp <- width
    width <- height
    height <- tmp
  }
  
  # extract extension
  ext <- sub(".*\\.", "", name)
  base <- sub("\\.[^.]*$", "", name)
  
  filename <- paste0("output/figures/", base, "_", Sys.Date(), ".", ext)
  
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    ...
  )
}


# Save df ####
save_df_csv <- function(df, name, na_as_blank = FALSE) {
  
  # extract extension
  ext <- sub(".*\\.", "", name)
  
  # extract base name
  base <- sub("\\.[^.]*$", "", name)
  
  # build filename
  filename <- paste0("output/data_frames/", base, "_", Sys.Date(), ".", ext)
  
  # Decide if NAs should be printed as blanks or kept as is
  na_value <- if (na_as_blank) "" else "NA"
  
  write.csv(df, file = filename, row.names = FALSE, na = na_value)
}
