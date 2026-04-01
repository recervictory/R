# plotting/save_plot.R
# Utility for saving ggplot/base R plot objects to disk in common image formats.

library(Cairo)
library(futile.logger)

# Save a plot object to plot_dir using the specified format (png, pdf, jpeg, tiff).
# File name is constructed from count, project_name, batch_name, file_name, and plot_type.
# Uses CairoPNG when available for higher-quality PNG output.
save_plot <- function(plot_object, plot_dir,
                      count = "00",
                      file_name = "image",
                      project_name = "",
                      format = "png",
                      batch_name = "",
                      plot_type = "scatter",
                      width = 800, height = 800,
                      dpi = 72, ...) {
  flog.info("😎 Function Name: save_plot")

  # Ensure the output directory exists
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Resolve file extension from format
  file_extension <- switch(format,
    png  = "png",
    pdf  = "pdf",
    jpeg = "jpg",
    jpg  = "jpg",
    tiff = "tiff",
    "png"
  )

  # Build the full file path
  file_name <- paste(count, project_name, batch_name, file_name, plot_type, sep = "_")
  file_name <- paste0(file_name, ".", file_extension)
  file_path <- file.path(plot_dir, file_name)
  file_path <- normalizePath(file_path, mustWork = FALSE)

  # Write the plot
  tryCatch({
    if (format == "png") {
      if (capabilities("cairo")) {
        Cairo::CairoPNG(file_path, width = width, height = height, res = dpi)
      } else {
        png(file_path, width = width, height = height, res = dpi)
      }
      print(plot_object)
      dev.off()
    } else if (format %in% c("jpeg", "jpg")) {
      jpeg(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    } else if (format == "tiff") {
      tiff(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    } else if (format == "pdf") {
      width_in  <- width  / dpi
      height_in <- height / dpi
      pdf(file_path, width = width_in, height = height_in)
      print(plot_object)
      dev.off()
    } else {
      png(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    }

    flog.info("%s file saved to %s", toupper(format), file_path)
  }, error = function(e) {
    flog.error("Error saving plot: %s", e$message)
    stop(e)
  })
}
