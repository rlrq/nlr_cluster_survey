library(stringr)

gen_split_extract <- function(sep, col, default_col = 1){
    output_func <- function(x){
        cols <- unlist(str_split(x, sep))
        if (length(cols) < col){
            return(cols[[default_col]])
        } else {
            return(cols[[col]])
        }
    }
    return(output_func)
}

split_extract <- function(str, sep, col){
    return(gen_split_extract(sep, col)(str))
}
