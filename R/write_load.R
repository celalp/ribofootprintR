write_profiling_data <- function(data, name) {
    list.save(x = data, file = paste(name, "rdata", sep = "."))
}

load_profiling_data <- function(file_name) {
    data <- list.load(file_name)
    data
}
