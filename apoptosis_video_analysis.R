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
samples_loaded_names <- str_extract(samples_to_load, '^(.*).csv$', group=1)
df_all <- readRDS('combined_csvs_processed.Rda')

df_all <- lapply(samples_to_load, read.csv)
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
  # df_all[,i] <- I(list(lapply(df_all[,i], function(x) unname(sapply(strsplit(str_replace(x, '\\[|\\]', ''), ', ')[[1]], as.numeric)))))
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
# dead_count <- pmin(unlist(lapply(deframe(df_all[,c(3,7)]), function(x) length(x))))
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
# temp <- df_all %>% filter(map_lgl(cell_id, ~ any(unlist(test['touching_cells']) %in% .x)))
# temp<- df_all[is.element(df_all$cell_id,  unlist(test['touching_cells'])),]
# this is now way too slow, I should change this to a LUT
df_all[, 'touching_mac'] <- NA
for (sample_id in 1:length(unique(df_all$.id))) {
  print(sample_id)
  temp_df <- df_all[df_all$.id == unique(df_all$.id)[sample_id],]
  temp_LUT <- temp_df[c('cell_id', 'mac_bool')]
  
  for (i in 1:nrow(temp_df)) {
    if (length(temp_df[i,'touching_cells'][[1]]) != 0) {
      # temp_df[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(temp_df[i,'touching_cells']), function(x) temp_df[temp_df$cell_id == x,]))['mac_bool'])))[,2]
      temp_df[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(temp_df[i,'touching_cells']), function(x) temp_LUT[x,][2])))))[2]
    }
  }
  df_all[df_all$.id == unique(df_all$.id)[sample_id],]$touching_mac <- temp_df$touching_mac
}
temp <- apply(df_all, 1, function(x) {ifelse(is.na(x['touching_mac']), list(), list(unlist(unname(x['touching_mac']))))})
temp <- lapply(temp, unlist)
temp <- lapply(temp, unname)
df_all['touching_mac'] <- enframe(temp)[,2]


# for now this is how we will handle parents/children:
# consider the parent as just when the children are co-existent
# they will inherit all the 'touching' as the parent,
# and then later we won't look at parent cells at all in analysis
# NOTE: The touching doesn't get inherited to the 'touched'. 
# I think this is ok...

# Since the loop goes in row order the first children should always run before
# their own children
# this one is also quite slow
# it also totally messes up the dtypes that we setup so nicely above...
for (sample_id in 1:length(unique(df_all$.id))) {
  print(sample_id)
  temp_df <- df_all[df_all$.id == unique(df_all$.id)[sample_id],]
  
  for (i in 1:nrow(temp_df)) {
    if (!is.na(temp_df[i,'parent_id'])) {
      parent_row <- temp_df[temp_df$cell_id == temp_df[i,'parent_id'],]
      temp_df[i,'touching_cells'] <- enframe(list(unlist(c(parent_row['touching_cells'][[1]], temp_df[i,'touching_cells'][[1]]))))[,2]
      temp_df[i,'touching_frames'] <- enframe(list(unlist(c(parent_row['touching_frames'][[1]], temp_df[i,'touching_frames'][[1]]))))[,2]
      temp_df[i,'neighbor_cells'] <- enframe(list(unlist(c(parent_row['neighbor_cells'][[1]], temp_df[i,'neighbor_cells'][[1]]))))[,2]
      temp_df[i,'neighbor_frames'] <- enframe(list(unlist(c(parent_row['neighbor_frames'][[1]], temp_df[i,'neighbor_frames'][[1]]))))[,2]
      temp_df[i,'touching_mac'] <- enframe(list(unlist(c(parent_row['touching_mac'][[1]], temp_df[i,'touching_mac'][[1]]))))[,2]
    }
  }
  df_all[df_all$.id == unique(df_all$.id)[sample_id],] <- temp_df
}

# df_all['touching_mac']<-enframe(lapply(df_all['touching_mac'][[1]], function(x) list(unlist(x))[[1]]))[,2]

temp <- apply(df_all, 1, function(x) {
  ifelse(is.null(unlist(x['touching_mac'])), list(), list(unlist(unname(Map(`[`, unlist(x['touching_frames']), unlist(x['touching_mac']))))))})
temp <- lapply(temp, unlist)
temp <- lapply(temp, unname)
df_all['mac_touch_frames'] <- enframe(temp)[,2]

temp <- apply(df_all, 1, function(x) {
  ifelse(is.null(unlist(x['touching_mac'])), list(), list(unlist(unname(Map(`[`, unlist(x['touching_frames']), !unlist(x['touching_mac']))))))})
temp <- lapply(temp, unlist)
temp <- lapply(temp, unname)
df_all['notmac_touch_frames'] <- enframe(temp)[,2]


# first mac touch
# this is so damn ugly my god
df_all['first_mac_touch'] <- unlist(lapply(df_all['touching_mac'][[1]], function(x) match('TRUE', x)))
df_all['first_mac_touch'] <- apply(df_all, 1, function(x) {
  ifelse(is.na(x['first_mac_touch']), NA, x['touching_frames'][[1]][ x['first_mac_touch'][[1]][[1]][1] ] )})

# how many touches
df_all['cell_touches'] <- unlist(lapply(df_all['touching_mac'][[1]], function(x) length(na.omit(x))))

# how many mac touches # sum only counts 'TRUE'
df_all['mac_touches'] <- unlist(lapply(df_all['touching_mac'][[1]], function(x) sum(na.omit(x))))

# how many not-mac touches
df_all['not_mac_touches'] <- df_all['cell_touches'] - df_all['mac_touches']

# how many spots per track
df_all['track_length'] <- unlist(lapply(df_all['frames'][[1]], function(x) length(x)))

# first apoptosis
df_all <- df_all %>% rowwise() %>% mutate(first_apoptosis = ifelse(!is.null(apoptotic_frames[1][[1]]), unlist(apoptotic_frames[1][[1]]),  NA))

# first touch
df_all <- df_all %>% rowwise() %>% mutate(first_touch = ifelse(!is.null(touching_frames[1][[1]]), unlist(touching_frames[1][[1]]),  NA))

# last mac touch before apoptosis

# t between mac touch and apoptosis
df_all <- df_all %>% mutate(reaper_time = case_when((dying_bool & !is.na(first_mac_touch)) ~ first_apoptosis - first_mac_touch,
                                                    TRUE ~ NA))
# t between any touch and apoptosis
df_all <- df_all %>% mutate(pre_death_touch_time = case_when((dying_bool & !is.na(first_touch)) ~ first_apoptosis - first_touch,
                                                             TRUE ~ NA))

df_all <- df_all %>% mutate(pre_death_touch_time = case_when((dying_bool & !is.na(first_touch)) ~ first_apoptosis - first_touch,
                                                             TRUE ~ NA))

# save the df so we don't have to redo all these slow loops
saveRDS(df_all, file='combined_csvs_processed.Rda') 


### Actual analysis ###

df_childless <- data.frame()
for (sample_id in 1:length(unique(df_all$.id))) {
  print(sample_id)
  temp_df <- df_all[df_all$.id == unique(df_all$.id)[sample_id],]
  
  parents <- unique(temp_df$parent_id)
  temp_df <- temp_df[-parents[!is.na(parents)],]
  df_childless <- rbind(df_childless, temp_df) 
}

# filter out the ones that have bad tracking / too many cells
#TODO move this function up to the part where I load in the csvs to make things 
# a lot faster & more programatic
bad_samples <- c('E2_atdc5_oc_24hr_17', 'E2_atdc5_oc_48hr_17', 'E2_atdc5_oc_17_combined')
df_childless <- rename(df_childless, id_col = .id)
df_childless <- df_childless %>% filter(! id_col %in% bad_samples)

### does touching make them die?
compiled_results_df <- data.frame(sample = unique(df_childless$id_col), total_cells = NA, total_dying = NA, dying_ratio = NA,
                                  total_mac = NA, mac_ratio = NA, total_mac_touch = NA,
                                  total_reaped = NA, reaper_ratio = NA, not_reaped_ratio = NA, 
                                  touched_not_mac_before_dying = NA, total_not_mac_touch = NA,
                                  touched_not_mac_dying_ratio = NA, reaped_ratio = NA, 
                                  total_touches = NA, mac_touches = NA, not_mac_touches = NA)

compiled_tall_results_mac <- data.frame(sample = unique(df_childless$id_col), cell_type = 'osteoclasts')
compiled_tall_results_notmac <- data.frame(sample = unique(df_childless$id_col), cell_type = 'chondrocytes')
compiled_tall_results_dying <- data.frame(sample = unique(df_childless$id_col), cell_type = 'dying')
compiled_tall_results_notdying <- data.frame(sample = unique(df_childless$id_col), cell_type = 'not_dying')
compiled_tall_results_reaped <- data.frame(sample = unique(df_childless$id_col), cell_type = 'reaped')
compiled_tall_results_notreaped <- data.frame(sample = unique(df_childless$id_col), cell_type = 'not_reaped')

for (sample_id in 1:length(unique(df_childless$id_col))) {
  temp_df <- df_childless[df_childless$id_col == unique(df_childless$id_col)[sample_id],]
  
  compiled_results_df[sample_id, 'total_cells'] <- nrow(temp_df)
  
  compiled_results_df[sample_id, 'total_dying'] <- temp_df %>% filter(dying_bool == TRUE) %>% nrow()
  compiled_results_df[sample_id, 'dying_ratio'] <- compiled_results_df[sample_id, 'total_dying'] / compiled_results_df[sample_id, 'total_cells'] * 100
  
  compiled_results_df[sample_id, 'total_mac'] <- temp_df %>% filter(mac_bool == TRUE) %>% nrow()
  compiled_results_df[sample_id, 'mac_ratio'] <- compiled_results_df[sample_id, 'total_mac'] / compiled_results_df[sample_id, 'total_cells'] * 100
  
  compiled_results_df[sample_id, 'total_mac_touch'] <- temp_df %>% filter(!is.na(first_mac_touch)) %>% nrow()
  compiled_results_df[sample_id, 'total_not_mac_touch'] <- temp_df %>% filter(is.na(first_mac_touch)) %>% nrow()
  compiled_results_df[sample_id, 'total_reaped'] <- temp_df %>% filter(reaper_time >= 0) %>% nrow()
  
  compiled_results_df[sample_id, 'reaper_ratio'] <- compiled_results_df[sample_id, 'total_reaped'] / compiled_results_df[sample_id, 'total_dying'] * 100
  
  
  compiled_results_df[sample_id, 'total_touches'] <- mean(temp_df$cell_touches)
  compiled_results_df[sample_id, 'mac_touches'] <- mean(temp_df$mac_touches)
  compiled_results_df[sample_id, 'not_mac_touches'] <- mean(temp_df$not_mac_touches)
  
  # likelihood of any cell dying w/o mac touch
  compiled_results_df[sample_id, 'not_reaped_ratio'] <- temp_df %>% filter(is.na(first_mac_touch) & (dying_bool == TRUE)) %>% nrow() / compiled_results_df[sample_id, 'total_cells'] * 100
  # likelihood of any cell dying touched by not mac
  compiled_results_df[sample_id, 'touched_not_mac_before_dying'] <- temp_df %>% filter(pre_death_touch_time >= 0 & (reaper_time < 0 | is.na(reaper_time))) %>% nrow()
  compiled_results_df[sample_id, 'not_mac_touch_ratio'] <- compiled_results_df[sample_id, 'total_not_mac_touch'] / (compiled_results_df[sample_id, 'total_cells'] - compiled_results_df[sample_id, 'total_mac'])
  compiled_results_df[sample_id, 'touched_not_mac_dying_ratio'] <- compiled_results_df[sample_id, 'touched_not_mac_before_dying'] / compiled_results_df[sample_id, 'total_not_mac_touch'] * 100
  # likelinhood of cells touched by macrophages dying
  compiled_results_df[sample_id, 'reaped_ratio'] <- compiled_results_df[sample_id, 'total_reaped'] / compiled_results_df[sample_id, 'total_mac_touch'] * 100
  

  
  temp_filter_df <- temp_df %>% filter(mac_bool == TRUE)
  compiled_tall_results_mac[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_mac[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_mac[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
  
  temp_filter_df <- temp_df %>% filter(mac_bool != TRUE)
  compiled_tall_results_notmac[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_notmac[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_notmac[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
  
  temp_filter_df <- temp_df %>% filter(dying_bool == TRUE)
  compiled_tall_results_dying[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_dying[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_dying[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
  
  temp_filter_df <- temp_df %>% filter(dying_bool != TRUE)
  compiled_tall_results_notdying[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_notdying[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_notdying[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
  
  temp_filter_df <- temp_df %>% filter(reaper_time >= 0)
  compiled_tall_results_reaped[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_reaped[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_reaped[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
  
  temp_filter_df <- temp_df %>% filter(is.na(first_mac_touch) & (dying_bool == TRUE))
  compiled_tall_results_notreaped[sample_id, 'total_touches'] <- mean(temp_filter_df$cell_touches)
  compiled_tall_results_notreaped[sample_id, 'mac_touches'] <- mean(temp_filter_df$mac_touches)
  compiled_tall_results_notreaped[sample_id, 'not_mac_touches'] <- mean(temp_filter_df$not_mac_touches)
}

# cell type
compiled_results_tall <- rbind(compiled_tall_results_mac, compiled_tall_results_notmac)
compiled_results_tall <- inner_join(compiled_results_tall, compiled_results_df[1:2])

compiled_results_taller <- compiled_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells), names_to = "touch_type", values_to = "touches")


ggplot(data=compiled_results_taller, aes(x = total_cells, y = touches, shape = cell_type, color = touch_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + ylab('touch events per cell')


ggplot(data = compiled_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')


# dying vs not dying
dying_results_tall <- rbind(compiled_tall_results_dying, compiled_tall_results_notdying)
dying_results_tall <- inner_join(dying_results_tall, compiled_results_df[1:2])

dying_results_taller <- dying_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells), names_to = "touch_type", values_to = "touches")

ggplot(data=dying_results_taller, aes(x = total_cells, y = touches, shape = cell_type, color = touch_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + ylab('touch events per cell')

ggplot(data = dying_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean", position='fill') + ylab('relative touch events per cell')


# reaped vs not reaped
reaped_results_tall <- rbind(compiled_tall_results_reaped, compiled_tall_results_notreaped)
reaped_results_tall <- inner_join(reaped_results_tall, compiled_results_df[1:2])

reaped_results_taller <- reaped_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells), names_to = "touch_type", values_to = "touches")

ggplot(data=reaped_results_taller, aes(x = total_cells, y = touches, shape = cell_type, color = touch_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + ylab('touch events per cell')

ggplot(data = reaped_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')


reaped_results_tall2 <- rbind(compiled_tall_results_reaped, compiled_tall_results_notdying)
reaped_results_tall2 <- inner_join(reaped_results_tall2, compiled_results_df[1:2])
reaped_results_taller2 <- reaped_results_tall2 %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells), names_to = "touch_type", values_to = "touches")

ggplot(data = reaped_results_taller2, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')
ggplot(data = reaped_results_taller2, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean", position='fill') + ylab('touch events per cell')


# reaper ratios
reaper_results_tall <- compiled_results_df[c(1,10,14)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'ratio')
reaper_results_tall <- reaper_results_tall %>% mutate(cell_type = case_when(cell_type == 'reaped_ratio' ~ "reaper ratio",
                                                                                cell_type == 'not_reaped_ratio' ~ 'not reaped ratio'))

ggplot(data=reaper_results_tall) + geom_bar(aes(cell_type, ratio), stat = "summary", fun.y = "mean")
t.test(ratio ~ cell_type, data=reaper_results_tall, paired = TRUE, alternative = "two.sided")



reaper_results_tall <- compiled_results_df[c(1,13,14)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'ratio')
reaper_results_tall <- reaper_results_tall %>% mutate(cell_type = case_when(cell_type == 'reaped_ratio' ~ "reaper ratio",
                                                                            cell_type == 'touched_not_mac_dying_ratio' ~ 'not reaped ratio'))

ggplot(data=reaper_results_tall) + geom_bar(aes(cell_type, ratio), stat = "summary", fun.y = "mean")
t.test(ratio ~ cell_type, data=reaper_results_tall, paired = TRUE, alternative = "two.sided")



touch_results_tall <- compiled_results_df[c(1,15,16)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'avg_touches')
touch_results_tall <- touch_results_tall %>% mutate(cell_type = case_when(cell_type == 'mac_touch_ratio' ~ "macrophage",
                                                                            cell_type == 'not_mac_touch_ratio' ~ 'not_macrophage'))

ggplot(data=touch_results_tall) + geom_bar(aes(cell_type, avg_touches), stat = "summary", fun.y = "mean")
t.test(avg_touches ~ cell_type, data=touch_results_tall, paired = TRUE, alternative = "two.sided")


# histograms
test2 <- df_childless[df_childless$id_col == unique(df_childless$id_col)[5],]
test2 <- test2 %>% filter(!is.na(first_apoptosis)) %>% select(first_apoptosis, mac_touch_frames, notmac_touch_frames)

df_childless$mac_t_diff <- apply(df_childless, 1, function(x) {x['mac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_t_diff <- apply(df_childless, 1, function(x) {x['notmac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_touch_count <- apply(df_childless, 1, function(x) {length(x['notmac_touch_frames'][[1]])})
df_childless$mac_touch_count <- apply(df_childless, 1, function(x) {length(x['mac_touch_frames'][[1]])})

tall_test <- data.frame()
mac_list <- list()
notmac_list <- list()
for (r in 1:nrow(df_childless)) {
  if (!is.null(df_childless[r,'mac_touch_frames'][[1]][[1]]) & length(df_childless[r,'mac_touch_frames'][[1]][[1]]) != 0) {
    temp_mac <- as.data.frame(append(enframe(unlist(df_childless[r,'mac_t_diff'][[1]]))[,2],
                                      c(df_childless[r,'mac_touch_count'], df_childless[r,'mac_touch_count']+df_childless[r, 'notmac_touch_count'])))
    mac_list <- append(mac_list, list(temp_mac))
  }
  if (!is.null(df_childless[r,'notmac_touch_frames'][[1]][[1]]) & length(df_childless[r,'notmac_touch_frames'][[1]][[1]]) != 0) {
    temp_notmac <- as.data.frame(append(enframe(unlist(df_childless[r,'notmac_t_diff'][[1]]))[,2],
                                     c(df_childless[r,'notmac_touch_count'], df_childless[r,'mac_touch_count']+df_childless[r, 'notmac_touch_count'])))
    notmac_list <- append(notmac_list, list(temp_notmac))
  }
}
mac_touch_test <- rbindlist(mac_list)
notmac_touch_test <- rbindlist(notmac_list)


ggplot() + geom_histogram(data=mac_touch_test, aes(x=value/(60/8), y=..density.., weight = 1/mac_touch_count), binwidth = 0.5, alpha=0.5, fill='red') + 
           geom_density(data=mac_touch_test, aes(x=value/(60/8), y=..density.., weight = 1/mac_touch_count), color='red') +
           geom_histogram(data=notmac_touch_test, aes(x=value/(60/8), y=..density.., weight = 1/notmac_touch_count), binwidth = 0.5, alpha=0.5, fill='black') + 
           geom_density(data=notmac_touch_test, aes(x=value/(60/8), y=..density.., weight = 1/notmac_touch_count))


reaper_histogram <- ggplot() + geom_histogram(data = df_childless %>% filter(reaper_time >= 0), aes(x=-reaper_time/(60/8)), alpha=0.5, fill='red') + 
  geom_histogram(data = df_childless %>% filter(pre_death_touch_time >= 0 & (reaper_time < 0 | is.na(reaper_time))), aes(x=-pre_death_touch_time/(60/8)), alpha=0.5, fill='black')

# idea https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
