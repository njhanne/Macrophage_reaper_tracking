library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

library(tibble)
library(data.table)

library(ggplot2)
library(cowplot)


define_cell_type <- function(temp_df, cell) {
  cell_count <- unlist(lapply(deframe(temp_df[,c('cell_id',cell)]), function(x) length(x)))
  cell_ratio <- unname(cell_count / temp_df$frame_count)
  return(cell_ratio)
}


# Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("./data/")

sample_info <- read.csv('image_log_out.csv')

# Load data
setwd("./results/tracks_csv/")

samples <- unique(sample_info$new_filename_timeless)
results_csvs <- list.files(pattern = ".csv$")
combined_csvs <- results_csvs %>% str_subset(pattern = "^.*_combined.csv")
drop_csvs <- combined_csvs %>% str_extract(pattern = "^(.*)_combined.csv$", group=1)

search_str <- str_c(drop_csvs, collapse='|')
uncombined_samples <- str_subset(samples, search_str, negate = TRUE)

samples_to_keep <- sample_info[sample_info$new_filename_timeless %in% uncombined_samples,]$new_filename
samples_to_keep <- str_c(samples_to_keep, '.csv')
samples_to_load <- samples_to_keep[samples_to_keep %in% results_csvs]
samples_to_load <- c(samples_to_load, combined_csvs)

# read them all into a list, these can be very slow
df_all <- lapply(samples_to_load, read.csv)
samples_loaded_names <- str_extract(samples_to_load, '^(.*).csv$', group=1)
df_all <- rbindlist(df_all, idcol=TRUE)
df_all$.id = samples_loaded_names[df_all$.id] # put file name into the df
df_all <- df_all[,!2] # gets rid of first useless added column
df_all <- as.data.frame(df_all) # convert datatable back to df

### Convert all lists to numeric lists
# columns_to_numlist <- c("spot_ids", "frames", "apoptotic_frames", "apoptotic_spots",
                        # "macrophage_frames", "macrophage_spots", "touching_frames",
                        # "touching_cells", "neighbor_frames", "neighbor_cells")
columns_to_numlist <- c(4:5,7:length(df_all))

for (i in columns_to_numlist) {
  # get rid of brackets on beginning and end:  str_replace_all(a, '\\[|\\]', '')
  # convert from one big string to substrings:  strsplit(str_replace_all(a, '\\[|\\]', ''), ', ')
  # convert to numeric, remove 'name' and make into a list:   list(unname(sapply(a, as.numeric)))
  # all together:  list(unname(sapply(strsplit(str_replace(a, '\\[|\\]', ''), ', ')[[1]], as.numeric)))
  # need the I(list()) part to make an 'asis' list, only way to store a list in df 
  # df_test[,i] <- I(list(lapply(df_test[,i], function(x) unname(sapply(strsplit(str_replace(x, '\\[|\\]', ''), ', ')[[1]], as.numeric)))))
  temp <- lapply(df_all[,i], function(x) unname(na.omit(sapply(strsplit(str_replace_all(x, '\\[|\\]', ''), ', ')[[1]], as.numeric))))
  df_all[,i] <- enframe(temp)[,2]
}

### cleanup data
# delete cells only present for one frame -- excluding children
df_all$parent_frame_count <- unlist(lapply(deframe(df_all[,c('cell_id','frames')]), function(x) length(x)))
df_all <- df_all %>% filter(parent_frame_count > 1)

# add in the actual frame count to include children
df_all$frame_count <- unlist(lapply(deframe(df_all[,c('cell_id','independent_frames')]), function(x) length(x)))
df_all <- df_all %>% mutate(frame_count = case_when(frame_count == 0 ~ parent_frame_count,
                                                      TRUE ~ frame_count))


### define cells
# macrophages
# I think if they are macrophage+ 1/2 of the time, they are macrophages.
# They likely won't be 100% of the time because of how the fluorescent 
# detection works with trackmate and python
# event if they are no longer fluorescing
df_all$mac_frame_ratio <- define_cell_type(df_all, 'macrophage_frames')
df_all <- df_all %>% mutate(mac_bool = case_when(mac_frame_ratio > 0.5 ~ TRUE,
                                                     .default = FALSE))

# define apoptotic
# These can be more brief, but still want them to be more than a few frames
df_all$dying_frame_ratio <- define_cell_type(df_all, 'apoptotic_frames')
# dead_count <- pmin(unlist(lapply(deframe(df_test[,c(3,7)]), function(x) length(x))))
df_all <- df_all %>% mutate(dying_bool = case_when(dying_frame_ratio > 0.3 ~ TRUE,
                                                   .default = FALSE))



### when touched
# right now touching is determined by two lists, one of all the cell ids they touch
# and one of the frame it occured on. They are in the same order... so:
# touching_cells: [3,6,1] and touching_frames: [1,1,2]
# means taht in frame 1 it is touching cell 3 and 6 and in frame 2 it touches cell 1
# we can just make another list of booleans of whether the touch was a mac or not

# we want the touch history to be transferred to children as well, but only from
# before they split... this is a fun complicated problem!
# actually it's even more complicated than that. All parents cease to exist after children split so probably should not
# even be considered cells. Also that means that children have technically been touching the entire parent line...

# first make the mac_touch bool
# temp <- df_test %>% filter(map_lgl(cell_id, ~ any(unlist(test['touching_cells']) %in% .x)))
# temp<- df_test[is.element(df_test$cell_id,  unlist(test['touching_cells'])),]
# this is now way too slow, I should change this to a LUT
for (i in 1:nrow(df_all)) {
  if (length(df_all[i,'touching_cells'][[1]]) != 0) {
    df_all[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(df_all[i,'touching_cells']), function(x) df_all[df_all$cell_id == x,]))['mac_bool'])))[,2]
  }
}

# for now this is how we will handle parents/children:
# consider the parent as just when the children are co-existent
# they will inherit all the 'touching' as the parent,
# and then later we won't look at parent cells at all in analysis
# NOTE: The touching doesn't get inherited to the 'touched'. 
# I think this is ok...

# Since the loop goes in row order the first children should always run before
# their own children
for (i in 1:nrow(df_test)) {
  if (!is.na(df_test[i,'parent_id'])) {
    parent_row <- df_test[df_test$cell_id == df_test[i,'parent_id'],]
    df_test[i,'touching_cells'] <- enframe(list(unlist(c(parent_row['touching_cells'][[1]], df_test[i,'touching_cells'][[1]]))))[,2]
    df_test[i,'touching_frames'] <- enframe(list(unlist(c(parent_row['touching_frames'][[1]], df_test[i,'touching_frames'][[1]]))))[,2]
    df_test[i,'neighbor_cells'] <- enframe(list(unlist(c(parent_row['neighbor_cells'][[1]], df_test[i,'neighbor_cells'][[1]]))))[,2]
    df_test[i,'neighbor_frames'] <- enframe(list(unlist(c(parent_row['neighbor_frames'][[1]], df_test[i,'neighbor_frames'][[1]]))))[,2]
    df_test[i,'touching_mac'] <- enframe(list(unlist(c(parent_row['touching_mac'][[1]], df_test[i,'touching_mac'][[1]]))))[,2]
  }
}

# first mac touch
# this is so damn ugly my god
df_test['first_mac_touch'] <- unlist(lapply(df_test['touching_mac'][[1]], function(x) match('TRUE', x[[1]])))
df_test['first_mac_touch'] <- apply(df_test, 1, function(x) {
  ifelse(is.na(x['first_mac_touch']), NA, x['touching_frames'][[1]][ x['first_mac_touch'][[1]][[1]][1] ] )})

# first apoptosis
df_test <- df_test %>% rowwise() %>% mutate(first_apoptosis = ifelse(!is.null(apoptotic_frames[1][[1]]), unlist(apoptotic_frames[1][[1]]),  NA))

# first touch
df_test <- df_test %>% rowwise() %>% mutate(first_touch = ifelse(!is.null(touching_frames[1][[1]]), unlist(touching_frames[1][[1]]),  NA))


# last mac touch before apoptosis

# t between mac touch and apoptosis
df_test <- df_test %>% mutate(reaper_time = case_when((!is.na(first_apoptosis) & !is.na(first_mac_touch)) ~ first_apoptosis - first_mac_touch,
                                                      TRUE ~ NA))
# t between any touch and apoptosis
df_test <- df_test %>% mutate(pre_death_touch_time = case_when((!is.na(first_apoptosis) & !is.na(first_touch)) ~ first_apoptosis - first_touch,
                                                      TRUE ~ NA))


parents <- unique(df_test$parent_id)
df_childless <- df_test[-parents[!is.na(parents)],]
### does touching make them die?
total_dying <- df_childless %>% filter(dying_bool == TRUE) %>% nrow()
dying_ratio <- total_dying / nrow(df_childless) * 100

total_mac <- df_childless %>% filter(mac_bool == TRUE) %>% nrow()
mac_ratio <- total_mac / nrow(df_childless) * 100

total_mac_touch <- df_childless %>% filter(!is.na(first_mac_touch)) %>% nrow()
total_reaped <- df_childless %>% filter(reaper_time >= 0) %>% nrow()

reaper_ratio <- total_reaped / total_dying * 100

# likelihood of any cell dying w/o mac touch
not_reaped_ratio <- df_childless %>% filter(is.na(first_mac_touch) & (dying_bool == TRUE)) %>% nrow() / nrow(df_childless) * 100
# likelihood of any cell dying touched by not mac
touched_not_mac_before_dying <- df_childless %>% filter(pre_death_touch_time >= 0 & (reaper_time < 0 | is.na(reaper_time))) %>% nrow()
total_not_mac_touch <- df_childless %>% filter(!is.na(first_touch) & (reaper_time < 0 | is.na(reaper_time))) %>% nrow()
touched_not_mac_dying_ratio <- touched_not_mac_before_dying / total_not_mac_touch * 100
# likelinhood of cells touched by macrophages dying
reaped_ratio <- total_reaped / total_mac_touch * 100



reaper_histogram <- ggplot() + geom_histogram(data = df_childless %>% filter(reaper_time >= 0), aes(x=-reaper_time/(60/8)), alpha=0.5, fill='red') + 
                 geom_histogram(data = df_childless %>% filter(pre_death_touch_time >= 0 & (reaper_time < 0 | is.na(reaper_time))), aes(x=-pre_death_touch_time/(60/8)), alpha=0.5, fill='black')

# idea https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html


