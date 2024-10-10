# sketch from 
# https://stackoverflow.com/questions/78607499/pass-parameters-to-quarto-render
# Load the Quarto library

render_one <- function(dir) {
  require(quarto)
  quarto::quarto_render(input = 'ensemble-pipeline.qmd',
                        output_file = paste0(basename(dir), '_report.html'),
                        execute_params = list(datadir = dir),
                        execute_dir = getwd())
}

dirs <- c("data/C57BL6J-638850_6800",
          "data/C57BL6J-638850_7200")

for (dir in dirs) {
  cat(dir,"\n")
  tryCatch({
    render_one(dir)
  }, error = function(e) {
    message(paste("Error in rendering report for", dir, ":", e))
  })
}

render_one()