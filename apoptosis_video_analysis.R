library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)


define_cell_type <- function(temp_df, cell) {
  cell_count <- unlist(lapply(deframe(temp_df[,c('cell_id',cell)]), function(x) length(x)))
  cell_ratio <- unname(cell_count / temp_df$frame_count)
  return(cell_ratio)
}


# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory


# Load data
df_test <- read.csv('test.csv')

### Convert all lists to numeric lists
# columns_to_numlist <- c("spot_ids", "frames", "apoptotic_frames", "apoptotic_spots", 
#                         "macrophage_frames", "macrophage_spots", "touching_frames", 
#                         "touching_cells", "neighbor_frames", "neighbor_cells")
columns_to_numlist <- c(4,6:16)

for (i in columns_to_numlist) {
  # get rid of brackets on beginning and end:  str_replace_all(a, '\\[|\\]', '')
  # convert from one big string to substrings:  strsplit(str_replace_all(a, '\\[|\\]', ''), ', ')
  # convert to numeric, remove 'name' and make into a list:   list(unname(sapply(a, as.numeric)))
  # all together:  list(unname(sapply(strsplit(str_replace(a, '\\[|\\]', ''), ', ')[[1]], as.numeric)))
  # need the I(list()) part to make an 'asis' list, only way to store a list in df 
  # df_test[,i] <- I(list(lapply(df_test[,i], function(x) unname(sapply(strsplit(str_replace(x, '\\[|\\]', ''), ', ')[[1]], as.numeric)))))
  temp <- lapply(df_test[,i], function(x) unname(na.omit(sapply(strsplit(str_replace_all(x, '\\[|\\]', ''), ', ')[[1]], as.numeric))))
  df_test[,i] <- enframe(temp)[,2]
}

### cleanup data
# delete cells only present for one frame -- excluding children
df_test$parent_frame_count <- unlist(lapply(deframe(df_test[,c('cell_id','frames')]), function(x) length(x)))
df_test <- df_test %>% filter(parent_frame_count > 1)

# add in the actual frame count to include children
df_test$frame_count <- unlist(lapply(deframe(df_test[,c('cell_id','independent_frames')]), function(x) length(x)))
df_test <- df_test %>% mutate(frame_count = case_when(frame_count == 0 ~ parent_frame_count,
                                                      TRUE ~ frame_count))


### define cells
# macrophages
# I think if they are macrophage+ 2/3 of the time, they are macrophages.
# They likely won't be 100% of the time because of how the fluorescent 
# detection works with trackmate and python
# TODO: need to make a way for them to stop being macrophages after a child
# event if they are no longer fluorescing
df_test$mac_frame_ratio <- define_cell_type(df_test, 'macrophage_frames')
df_test <- df_test %>% mutate(mac_bool = case_when(mac_frame_ratio > 0.66 ~ TRUE,
                                                     .default = FALSE))

total_mac <- df_test %>% filter(mac_bool == TRUE) %>% nrow()
mac_ratio <- total_mac / nrow(df_test) * 100

# define apoptotic
# These can be more brief flashes. In the future I may make it where they need to
# stay lit for a bit, but for now we'll just do a straight binary
df_test$dying_frame_ratio <- define_cell_type(df_test, 'apoptotic_frames')
# dead_count <- pmin(unlist(lapply(deframe(df_test[,c(3,7)]), function(x) length(x))))
df_test <- df_test %>% mutate(dying_bool = case_when(dying_frame_ratio > 0 ~ TRUE,
                                                   .default = FALSE))
total_dying <- df_test %>% filter(dying_bool == TRUE) %>% nrow()
dying_ratio <- total_dying / nrow(df_test) * 100


### when touched
# right now touching is determined by two lists, one of all the cell ids they touch
# and one of the frame it occured on. They are in the same order... so:
# touching_cells: [3,6,1] and touching_frames: [1,1,2]
# means taht in frame 1 it is touching cell 3 and 6 and in frame 2 it touches cell 1
# we can just make another list of booleans of whether the touch was a mac or not
# we want the touch history to be transferred to children as well, but only from
# before they split... this is a fun complicated problem!

# first make the mac_touch bool
# temp <- df_test %>% filter(map_lgl(cell_id, ~ any(unlist(test['touching_cells']) %in% .x)))
# temp<- df_test[is.element(df_test$cell_id,  unlist(test['touching_cells'])),]
for (i in 1:nrow(df_test)) {
  if (length(df_test[i,'touching_cells'][[1]]) != 0) {
    df_test[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(df_test[i,'touching_cells']), function(x) df_test[df_test$cell_id == x,]))['mac_bool'])))[,2]
  }
}






