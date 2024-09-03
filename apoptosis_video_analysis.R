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


process_batch_data <- function(batch_data) {
  ### Convert all lists to numeric lists
  # columns_to_numlist <- c("spot_ids", "frames", "apoptotic_frames", "apoptotic_spots",
  # "macrophage_frames", "macrophage_spots", "touching_frames",
  # "touching_cells", "neighbor_frames", "neighbor_cells")
  columns_to_numlist <- c(4:6,8:length(batch_data))
  
  for (i in columns_to_numlist) {
    # get rid of brackets on beginning and end:  str_replace_all(a, '\\[|\\]', '')
    # convert from one big string to substrings:  strsplit(str_replace_all(a, '\\[|\\]', ''), ', ')
    # convert to numeric, remove 'name' and make into a list:   list(unname(sapply(a, as.numeric)))
    # all together:  list(unname(sapply(strsplit(str_replace(a, '\\[|\\]', ''), ', ')[[1]], as.numeric)))
    # need the I(list()) part to make an 'asis' list, only way to store a list in df 
    # batch_data[,i] <- I(list(lapply(batch_data[,i], function(x) unname(sapply(strsplit(str_replace(x, '\\[|\\]', ''), ', ')[[1]], as.numeric)))))
    temp <- lapply(batch_data[,i], function(x) unname(na.omit(sapply(strsplit(str_replace_all(x, '\\[|\\]', ''), ', ')[[1]], as.numeric))))
    batch_data[,i] <- enframe(temp)[,2]
  }
  
  ### cleanup data
  # delete cells only present for one frame -- excluding children
  batch_data$parent_frame_count <- unlist(lapply(deframe(batch_data[,c('cell_id','frames')]), function(x) length(x)))
  batch_data <- batch_data %>% filter(parent_frame_count > 1)
  
  # add in the actual frame count to include children
  batch_data$frame_count <- unlist(lapply(deframe(batch_data[,c('cell_id','independent_frames')]), function(x) length(x)))
  batch_data <- batch_data %>% mutate(frame_count = case_when(frame_count == 0 ~ parent_frame_count,
                                                      TRUE ~ frame_count))
  
  
  ### define cells
  # macrophages
  # I think if they are macrophage+ 1/2 of the time, they are macrophages.
  # They likely won't be 100% of the time because of how the fluorescent 
  # detection works with trackmate and python
  # event if they are no longer fluorescing
  batch_data$mac_frame_ratio <- define_cell_type(batch_data, 'macrophage_frames')
  batch_data <- batch_data %>% mutate(mac_bool = case_when(mac_frame_ratio > 0.5 ~ TRUE,
                                                   .default = FALSE))
  
  # define apoptotic
  # These can be more brief, but still want them to be more than a few frames
  batch_data$dying_frame_ratio <- define_cell_type(batch_data, 'apoptotic_frames')
  # dead_count <- pmin(unlist(lapply(deframe(batch_data[,c(3,7)]), function(x) length(x))))
  batch_data <- batch_data %>% mutate(dying_bool = case_when(dying_frame_ratio > 0.3 ~ TRUE,
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
  # temp <- batch_data %>% filter(map_lgl(cell_id, ~ any(unlist(test['touching_cells']) %in% .x)))
  # temp<- batch_data[is.element(batch_data$cell_id,  unlist(test['touching_cells'])),]
  # this is now way too slow, I should change this to a LUT
  batch_data[, 'touching_mac'] <- NA
  for (sample_id in 1:length(unique(batch_data$.id))) {
    print(sample_id)
    temp_df <- batch_data[batch_data$.id == unique(batch_data$.id)[sample_id],]
    temp_LUT <- temp_df[c('cell_id', 'mac_bool')]
    
    for (i in 1:nrow(temp_df)) {
      if (length(temp_df[i,'touching_cells'][[1]]) != 0) {
        # temp_df[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(temp_df[i,'touching_cells']), function(x) temp_df[temp_df$cell_id == x,]))['mac_bool'])))[,2]
        temp_df[i,'touching_mac'] <- enframe(list(unname(do.call("rbind",lapply(unlist(temp_df[i,'touching_cells']), function(x) temp_LUT[x,][2])))))[2]
      }
    }
    batch_data[batch_data$.id == unique(batch_data$.id)[sample_id],]$touching_mac <- temp_df$touching_mac
  }
  temp <- apply(batch_data, 1, function(x) {ifelse(is.na(x['touching_mac']), list(), list(unlist(unname(x['touching_mac']))))})
  temp <- lapply(temp, unlist)
  temp <- lapply(temp, unname)
  batch_data['touching_mac'] <- enframe(temp)[,2]
  
  
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
  for (sample_id in 1:length(unique(batch_data$.id))) {
    print(sample_id)
    temp_df <- batch_data[batch_data$.id == unique(batch_data$.id)[sample_id],]
    
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
    batch_data[batch_data$.id == unique(batch_data$.id)[sample_id],] <- temp_df
  }
  
  # batch_data['touching_mac']<-enframe(lapply(batch_data['touching_mac'][[1]], function(x) list(unlist(x))[[1]]))[,2]
  
  temp <- apply(batch_data, 1, function(x) {
    ifelse(is.null(unlist(x['touching_mac'])), list(), list(unlist(unname(Map(`[`, unlist(x['touching_frames']), unlist(x['touching_mac']))))))})
  temp <- lapply(temp, unlist)
  temp <- lapply(temp, unname)
  batch_data['mac_touch_frames'] <- enframe(temp)[,2]
  
  temp <- apply(batch_data, 1, function(x) {
    ifelse(is.null(unlist(x['touching_mac'])), list(), list(unlist(unname(Map(`[`, unlist(x['touching_frames']), !unlist(x['touching_mac']))))))})
  temp <- lapply(temp, unlist)
  temp <- lapply(temp, unname)
  batch_data['notmac_touch_frames'] <- enframe(temp)[,2]
  
  
  # first mac touch
  # this is so damn ugly my god
  batch_data['first_mac_touch'] <- unlist(lapply(batch_data['touching_mac'][[1]], function(x) match('TRUE', x)))
  batch_data['first_mac_touch'] <- apply(batch_data, 1, function(x) {
    ifelse(is.na(x['first_mac_touch']), NA, x['touching_frames'][[1]][ x['first_mac_touch'][[1]][[1]][1] ] )})
  
  # how many touches
  batch_data['cell_touches'] <- unlist(lapply(batch_data['touching_mac'][[1]], function(x) length(na.omit(x))))
  
  # how many mac touches # sum only counts 'TRUE'
  batch_data['mac_touches'] <- unlist(lapply(batch_data['touching_mac'][[1]], function(x) sum(na.omit(x))))
  
  # how many not-mac touches
  batch_data['not_mac_touches'] <- batch_data['cell_touches'] - batch_data['mac_touches']
  
  # how many spots per track
  batch_data['track_length'] <- unlist(lapply(batch_data['frames'][[1]], function(x) length(x)))
  
  # first apoptosis
  batch_data <- batch_data %>% rowwise() %>% mutate(first_apoptosis = ifelse(!is.null(apoptotic_frames[1][[1]]), unlist(apoptotic_frames[1][[1]]),  NA))
  
  # first touch
  batch_data <- batch_data %>% rowwise() %>% mutate(first_touch = ifelse(!is.null(touching_frames[1][[1]]), unlist(touching_frames[1][[1]]),  NA))
  
  # last mac touch before apoptosis
  
  # t between mac touch and apoptosis
  batch_data <- batch_data %>% mutate(reaper_time = case_when((dying_bool & !is.na(first_mac_touch)) ~ first_apoptosis - first_mac_touch,
                                                      TRUE ~ NA))
  # t between any touch and apoptosis
  batch_data <- batch_data %>% mutate(pre_death_touch_time = case_when((dying_bool & !is.na(first_touch)) ~ first_apoptosis - first_touch,
                                                               TRUE ~ NA))
  
  # batch_data <- batch_data %>% mutate(pre_death_touch_time = case_when((dying_bool & !is.na(first_touch)) ~ first_apoptosis - first_touch,
  #                                                              TRUE ~ NA))
  return(batch_data)
}


#### MAIN ####
#### 0.1 Directory ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd() #check our working directory
setwd("./data/")

sample_info <- read.csv('image_log_out.csv')

#### 0.2 Load data ####
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

#### 0.2.1 Load previous batches ####
# here we need some code to handle different batches, we don't want to rerun all 
# of them every time, the code is quite slow
setwd('../.')
completed_batches <- list.files(pattern = ".Rda$")
batch_nums_to_load <- c('1', '2', '3', '4')
completed_batch_nums <- str_extract(completed_batches, '(?<=batch_)(.*)(?=.Rda)', group=1)
batches_to_load <- list()
for (batch_num in batch_nums_to_load) {
  match_batch <- str_c('batch_', batch_num, '.Rda')
  if (batch_num %in% completed_batch_nums) {
    batches_to_load <- append(batches_to_load, match_batch)
  }
}

# load in all the batches as a list, then rbind them into a df
# this step will be slow
df_all_list <- lapply(batches_to_load, readRDS)
df_all <- rbindlist(df_all_list, fill=TRUE)

# now remove the already analyzed ones from the 'to analyze' list
samples_to_drop <- lapply(unique(df_all[,1]), function(x) str_c(x, '.csv'))[[1]]
samples_to_keep <- !(samples_to_load %in% samples_to_drop)
samples_to_load <- samples_to_load[samples_to_keep] 

#### 0.2.2 Load new batches ####
setwd('./tracks_csv/')

batch_nums_to_analyze <- c('4')
for (batch_num in batch_nums_to_analyze) {
  # get all the filenames, add in the underscore, and remove the first match which is NA
  batch_samples <- str_c(unique(sample_info[sample_info$processing_batch == batch_num,]$new_filename_timeless), '_')[-1] 
  batch_samples_bool <- str_starts(samples_to_load, str_c(batch_samples, collapse='|'))
  batch_samples <- samples_to_load[batch_samples_bool]
  
  # load the csvs and rbind them into a df, this is slow!
  df_batch <- lapply(batch_samples, read.csv)
  df_batch <- rbindlist(df_batch, idcol=TRUE)
  samples_loaded_names <- str_extract(batch_samples, '^(.*).csv$', group=1)
  df_batch$.id = samples_loaded_names[df_batch$.id] # put file name into the df
  df_batch <- df_batch[,!2] # gets rid of first useless added column
  df_batch <- as.data.frame(df_batch) # convert datatable back to df
  
  # process the data into a functional df
  # makes sense to put this in a function and clean up code
  df_batch <- process_batch_data(df_batch)
  
  # save the batch info
  setwd('../.')
  savename <- str_c('batch_',batch_num,'.Rda')
  saveRDS(df_batch, file=savename)
  setwd('./tracks_csv/')
  
  # add the new batch to combined df_all
  df_all <- append(df_all_list, list(as_tibble(df_batch)))
}

# finalize loading all the lists of dfs
df_all <- rbindlist(df_all_list, fill=TRUE)


#### 1.0 Actual analysis ####

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


#TODO rerun some of the samples on python to complete the density measures
#TODO redo the thresholding in the monoculture samples


compiled_confluence_df <- data.frame(sample = unique(df_childless$id_col), cell_area = NA, avg_cell_num = NA)
for (sample_num in 1:length(unique(df_childless$id_col))) {
  sample_id <- unique(df_childless$id_col)[sample_num]
  if (grepl('combined', sample_id)) {
    match <- sample_info[sample_info$new_filename_timeless == sub('_combined', '', sample_id),]
  } else if (grepl('24hr', sample_id) | grepl('48hr', sample_id)) {
    match <- sample_info[sample_info$new_filename == sample_id,]
  }
  else {
    match <- NA
  }
  compiled_confluence_df[sample_num, 'cell_area'] <- mean(match$avg_cell_area)
  compiled_confluence_df[sample_num, "avg_cell_num"] <- mean(match$avg_cell_count)
}

### does touching make them die?
compiled_results_df <- data.frame(sample = unique(df_childless$id_col), total_cells = NA, total_dying = NA, dying_ratio = NA,
                                  total_mac = NA, mac_ratio = NA, touched_mac = NA, 
                                  not_touched_mac = NA, total_reaped = NA, reaper_ratio = NA, 
                                  not_reaped_ratio = NA, touched_not_mac_before_dying = NA, 
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
  compiled_results_df[sample_id, 'dying_ratio'] <- compiled_results_df[sample_id, 'total_dying'] / compiled_results_df[sample_id, 'total_cells']
  
  compiled_results_df[sample_id, 'total_mac'] <- temp_df %>% filter(mac_bool == TRUE) %>% nrow()
  compiled_results_df[sample_id, 'mac_ratio'] <- compiled_results_df[sample_id, 'total_mac'] / compiled_results_df[sample_id, 'total_cells'] 
  
  # number of cells that did / did not touch macrophages
  compiled_results_df[sample_id, 'touched_mac'] <- temp_df %>% filter(!is.na(first_mac_touch)) %>% nrow()
  compiled_results_df[sample_id, 'not_touched_mac'] <- temp_df %>% filter(is.na(first_mac_touch)) %>% nrow()
  compiled_results_df[sample_id, 'total_reaped'] <- temp_df %>% filter(reaper_time >= 0) %>% nrow()
  
  compiled_results_df[sample_id, 'reaper_ratio'] <- 0
  if (compiled_results_df[sample_id, 'total_dying'] !=0){
    compiled_results_df[sample_id, 'reaper_ratio'] <- compiled_results_df[sample_id, 'total_reaped'] / compiled_results_df[sample_id, 'total_dying']
  }
  
  compiled_results_df[sample_id, 'total_touches'] <- mean(temp_df$cell_touches)
  compiled_results_df[sample_id, 'mac_touches'] <- mean(temp_df$mac_touches)
  compiled_results_df[sample_id, 'not_mac_touches'] <- mean(temp_df$not_mac_touches)
  

  # likelihood of any cell dying w/o mac touch
  compiled_results_df[sample_id, 'not_reaped_ratio'] <- temp_df %>% filter(is.na(first_mac_touch) & (dying_bool == TRUE)) %>% nrow() / compiled_results_df[sample_id, 'total_cells'] 
  # likelihood of any cell dying touched by not mac
  compiled_results_df[sample_id, 'touched_not_mac_before_dying'] <- temp_df %>% filter(pre_death_touch_time >= 0 & (reaper_time < 0 | is.na(reaper_time))) %>% nrow()
  compiled_results_df[sample_id, 'not_mac_touch_ratio'] <- compiled_results_df[sample_id, 'not_touched_mac'] / (compiled_results_df[sample_id, 'total_cells'] - compiled_results_df[sample_id, 'total_mac'])
  compiled_results_df[sample_id, 'touched_not_mac_dying_ratio'] <- compiled_results_df[sample_id, 'touched_not_mac_before_dying'] / compiled_results_df[sample_id, 'not_touched_mac']
  
  compiled_results_df[sample_id, 'reaped_ratio'] <- compiled_results_df[sample_id, 'total_reaped'] / compiled_results_df[sample_id, 'touched_mac'] 
  
  #TODO this needs to be moved. this is terrible
  if (grepl('atdc5', compiled_results_df[sample_id, 'sample']) & !grepl('atdc5_oc', compiled_results_df[sample_id, 'sample'])) {
    compiled_results_df[sample_id, 'total_mac'] <- 0
    compiled_results_df[sample_id, 'mac_ratio'] <- 0
    compiled_results_df[sample_id, 'touched_mac'] <- 0
    compiled_results_df[sample_id, 'total_reaped'] <- 0
    compiled_results_df[sample_id, 'reaper_ratio'] <- 0
    compiled_results_df[sample_id, 'not_mac_touches'] <- compiled_results_df[sample_id, 'total_touches'] 
    compiled_results_df[sample_id, 'mac_touches'] <- 0
    compiled_results_df[sample_id, 'reaped_ratio'] <- 0
  }

  
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

# cleanup the compiled data
compiled_results_df <- compiled_results_df %>% mutate(cell_type = case_when(grepl('atdc5_oc', sample) ~ 'ATDC5_coculture',
                                                                            grepl('callus_oc', sample) ~ 'Callus_coculture',
                                                                            grepl('atdc5', sample) ~ 'ATDC5',
                                                                            grepl('callus', sample) ~ 'Callus',
                                                                            grepl('oc', sample) ~ 'Osteoclasts',
                                                                            TRUE ~ 'Other'))


# add data from info_csv into our compiled results
compiled_results_df <- left_join(compiled_results_df, compiled_confluence_df)


### Graph helping ###
pal <- c('#cc3311', '#bbbbbb', '#ee7733', '#0077bb')


# test different confluence calcs
ggplot(data=compiled_results_df, aes(x = total_cells, y = avg_cell_num)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=FALSE)

ggplot(data=compiled_results_df, aes(x = total_cells, y = cell_area)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=FALSE)

ggplot(data=compiled_results_df, aes(x = avg_cell_num, y = cell_area)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=FALSE)


### Probabilities / Bayes
# we can see if the probabilities are 'independent' if they multiply together
# independent if P(mactouch and death) = P(mactouch)*P(death)
compiled_results_df$reaped_independence <- (compiled_results_df$touched_mac / compiled_results_df$total_cells)*(compiled_results_df$total_dying / compiled_results_df$total_cells)
compiled_results_df$reap_prob <- compiled_results_df$total_reaped / compiled_results_df$total_cells

compiled_results_tall <- compiled_results_df %>% select(sample, cell_type, reap_prob, reaped_independence) 
compiled_results_tall <- compiled_results_tall %>% pivot_longer(-c(cell_type, sample), names_to = "prob_type", values_to = "probability")
compiled_results_tall$cell_type <- factor(compiled_results_tall$cell_type)


ggplot(data=compiled_results_tall, aes(cell_type, probability, alpha=prob_type, fill=cell_type)) +
  geom_bar(position='dodge', stat = "summary", fun = "mean") +
  geom_jitter(position = position_jitterdodge(0.1)) +
  ylab('probability') +
  scale_fill_manual(values=pal) +
  scale_alpha_manual(values=c(.66, 1))

for (group_i in 1:length(levels(compiled_results_tall$cell_type))) {
  temp_df <- compiled_results_tall %>% filter(cell_type == levels(compiled_results_tall$cell_type)[group_i])
  print(levels(compiled_results_tall$cell_type)[group_i])
  print(t.test(probability ~ prob_type, data=temp_df, paired = TRUE, alternative = "two.sided"))
}
# the atdc5 coculture appears to be independent...


# P(death given observed mactouch) = x = P(death) * P(observed death after mactouch) / P(mactouch)
# P(mactouch and notdeath) = y
# P(nomactouch and death) = a = P(death) * P(observed death after not touching mac) / P(notmactouch)
# P(nomactouch and notdeath) = b
compiled_results_df$bayes_prob <- (compiled_results_df$total_dying / compiled_results_df$total_cells) * compiled_results_df$reap_prob / (compiled_results_df$touched_mac / compiled_results_df$total_cells)
compiled_results_df$bayes_prob_no <- (compiled_results_df$total_dying / compiled_results_df$total_cells) * compiled_results_df$not_reaped_ratio / (compiled_results_df$not_touched_mac / compiled_results_df$total_cells)

compiled_results_tall <- compiled_results_df %>% select(sample, cell_type, bayes_prob_no, bayes_prob) 
compiled_results_tall <- compiled_results_tall %>% pivot_longer(-c(cell_type, sample), names_to = "prob_type", values_to = "probability")
compiled_results_tall$cell_type <- factor(compiled_results_tall$cell_type)

ggplot(data=compiled_results_tall, aes(cell_type, probability*100, alpha=prob_type, fill=cell_type)) +
  geom_bar(position='dodge', stat = "summary", fun = "mean") +
  geom_jitter(position = position_jitterdodge(0.1)) +
  ylab('probability %') +
  scale_fill_manual(values=pal) +
  scale_alpha_manual(values=c(.66, 1))

for (group_i in 2:length(levels(compiled_results_tall$cell_type))) {
  temp_df <- compiled_results_tall %>% filter(cell_type == levels(compiled_results_tall$cell_type)[group_i])
  print(levels(compiled_results_tall$cell_type)[group_i])
  print(t.test(probability ~ prob_type, data=temp_df, paired = TRUE, alternative = "two.sided"))
}
# bayesian says touching macrophage does increase risk of dying (except in monoculture)



### cell touching analysis
# cell type
compiled_results_tall <- rbind(compiled_tall_results_mac, compiled_tall_results_notmac)
compiled_results_tall <- inner_join(compiled_results_tall, compiled_results_df[c(1,2,19,20)])

compiled_results_taller <- compiled_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, cell_area, total_cells, avg_cell_num), names_to = "touch_type", values_to = "touches")
compiled_results_taller <- compiled_results_taller %>% mutate(cell_type = case_when(cell_type == 'chondrocytes' & grepl('atdc5_oc', sample) ~ 'ATDC5_coculture',
                                                                                    cell_type == 'chondrocytes' & grepl('callus_oc', sample) ~ 'Callus_coculture',
                                                                                    cell_type == 'chondrocytes' & grepl('atdc5', sample) ~ 'ATDC5',
                                                                                    cell_type == 'Osteoclasts' & grepl('callus_oc', sample) ~ 'Osteoclasts_coculture',
                                                                                    cell_type == 'Osteoclasts' & grepl('atdc5_oc', sample) ~ 'Osteoclasts_coculture',
                                                                                    cell_type == 'Osteoclasts' ~ 'Osteoclasts', 
                                                                                    TRUE ~ 'Other'))


ggplot(data=compiled_results_taller, aes(x = total_cells, y = touches, shape = cell_type, color = touch_type, linetype=cell_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=FALSE) + ylab('touch events per cell')


ggplot(data = compiled_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')


# touching for dying vs not dying
dying_results_tall <- rbind(compiled_tall_results_dying, compiled_tall_results_notdying)
dying_results_tall <- inner_join(dying_results_tall, compiled_results_df[c(1,2,19,20)])

dying_results_taller <- dying_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells, cell_area, avg_cell_num), names_to = "touch_type", values_to = "touches")

ggplot(data=dying_results_taller, aes(x = total_cells, y = touches, shape = cell_type, color = touch_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + ylab('touch events per cell')

ggplot(data = dying_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean", position='fill') + ylab('relative touch events per cell')



### Death analysis
# dying ratio by cell
compiled_results_df$cell_type <- factor(compiled_results_df$cell_type, levels=c('Osteoclasts', 'ATDC5', 'ATDC5_coculture', 'Callus_coculture'))

ggplot(data=compiled_results_df, aes(x=cell_type, y = dying_ratio*100, fill = cell_type)) +
  geom_bar(stat='summary', fun='mean') + 
  geom_jitter( width=0.1) +
  ylab('Relative apoptosis (% of cells)') +
  scale_fill_manual(values=pal)

ggplot(data=compiled_results_df, aes(x = total_cells, y = dying_ratio*100, color = cell_type)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE) +
  ylab('Relative apoptosis (% of cells)') +
  xlab('Total number of cells ~ confluency') +
  scale_color_manual(values=pal)

death_by_culture <- aov(dying_ratio~cell_type, data=compiled_results_df)
summary(death_by_culture)
TukeyHSD(death_by_culture)

# reaped vs not reaped
reaped_results_tall <- rbind(compiled_tall_results_reaped, compiled_tall_results_notreaped)
reaped_results_tall <- inner_join(reaped_results_tall, compiled_results_df[c(1,2,19,20)], by='sample')

reaped_results_taller <- reaped_results_tall %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells, cell_area, avg_cell_num), names_to = "touch_type", values_to = "touches")

ggplot(data=reaped_results_taller, aes(x = cell_area, y = touches, shape = cell_type, color = touch_type)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + ylab('touch events per cell')

ggplot(data = reaped_results_taller, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')


reaped_results_tall2 <- rbind(compiled_tall_results_reaped, compiled_tall_results_notdying)
reaped_results_tall2 <- inner_join(reaped_results_tall2, compiled_results_df[c(1,2,19,20)])
reaped_results_taller2 <- reaped_results_tall2 %>% pivot_longer(.,-c(total_touches, cell_type, sample, total_cells, cell_area, avg_cell_num), names_to = "touch_type", values_to = "touches")

ggplot(data = reaped_results_taller2, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean") + ylab('touch events per cell')
ggplot(data = reaped_results_taller2, aes(x = cell_type, y = touches, fill = touch_type))+
  geom_bar(stat = "summary", fun = "mean", position='fill') + ylab('touch events per cell')


# reaper ratios
reaper_results_tall <- compiled_results_df[c(1,11,14)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'ratio')
reaper_results_tall <- reaper_results_tall %>% mutate(cell_type = case_when(cell_type == 'reaped_ratio' ~ "reaper ratio",
                                                                                cell_type == 'not_reaped_ratio' ~ 'not reaped ratio'))

ggplot(data=reaper_results_tall) + geom_bar(aes(cell_type, ratio), stat = "summary", fun.y = "mean")
t.test(ratio ~ cell_type, data=reaper_results_tall, paired = TRUE, alternative = "two.sided")



# reap vs never touched mac and dies ## primary result ##
reaper_results_tall <- compiled_results_df[c(1,13,14)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'ratio')
reaper_results_tall <- reaper_results_tall %>% mutate(reap_type = case_when(cell_type == 'reaped_ratio' ~ "reaped",
                                                                            cell_type == 'touched_not_mac_dying_ratio' ~ 'not reaped'))
reaper_results_tall <- reaper_results_tall %>% mutate(cell_type = case_when(grepl('atdc5_oc', sample) ~ 'ATDC5_coculture',
                                                                            grepl('callus_oc', sample) ~ 'Callus_coculture',
                                                                            grepl('atdc5', sample) ~ 'ATDC5',
                                                                            grepl('callus_oc', sample) ~ 'Osteoclasts_coculture',
                                                                            grepl('atdc5_oc', sample) ~ 'Osteoclasts_coculture',
                                                                            grepl('oc', sample) ~ 'Osteoclasts',
                                                                            TRUE ~ 'Other'))


reaper_results_tall$cell_type <- factor(reaper_results_tall$cell_type, levels=c('Osteoclasts', 'ATDC5', 'ATDC5_coculture', 'Callus_coculture'))

ggplot(data=reaper_results_tall, aes(cell_type, ratio*100, alpha=reap_type, fill=cell_type)) +
  geom_bar(position='dodge', stat = "summary", fun = "mean") +
  geom_jitter(position = position_jitterdodge(0.1)) +
  ylab('% of cells') +
  scale_fill_manual(values=pal) +
  scale_alpha_manual(values=c(.66, 1))

for (group_i in 1:length(levels(reaper_results_tall$cell_type))) {
  temp_df <- reaper_results_tall %>% filter(cell_type == levels(reaper_results_tall$cell_type)[group_i])
  print(levels(reaper_results_tall$cell_type)[group_i])
  print(t.test(ratio ~ reap_type, data=temp_df, paired = TRUE, alternative = "two.sided"))
}

test <- reaper_results_tall %>% group_by(cell_type, reap_type) %>% dplyr::summarize(Mean = mean(ratio, na.rm=TRUE))



#TODO relook at what this is
touch_results_tall <- compiled_results_df[c(1,14,18)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'avg_touches')
touch_results_tall <- touch_results_tall %>% mutate(cell_type = case_when(cell_type == 'reaped_ratio' ~ "macrophage",
                                                                            cell_type == 'not_mac_touch_ratio' ~ 'not_macrophage'))

ggplot(data=touch_results_tall) + geom_bar(aes(cell_type, avg_touches), stat = "summary", fun.y = "mean")
t.test(avg_touches ~ cell_type, data=touch_results_tall, paired = TRUE, alternative = "two.sided")


# probability of reap
df_childless$mac_t_diff <- apply(df_childless, 1, function(x) {x['mac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_t_diff <- apply(df_childless, 1, function(x) {x['notmac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_touch_count <- apply(df_childless, 1, function(x) {length(x['notmac_touch_frames'][[1]])})
df_childless$mac_touch_count <- apply(df_childless, 1, function(x) {length(x['mac_touch_frames'][[1]])})

# if touched by cell type, what is probability of it resulting in dying?
# total touches by cell type that are reap / total touches by cell type
df_childless$mac_reap_count <- apply(df_childless, 1, function(x) {length(x['mac_t_diff'][[1]][!is.na(x['mac_t_diff'][[1]]) & x['mac_t_diff'][[1]] <= 0])})
df_childless$notmac_reap_count <- apply(df_childless, 1, function(x) {length(x['notmac_t_diff'][[1]][!is.na(x['notmac_t_diff'][[1]]) & x['notmac_t_diff'][[1]] <= 0])})

cell_reaper_probs <- df_childless %>% group_by(id_col) %>% summarise(mac_reap_count = sum(mac_reap_count), mac_touch_count = sum(mac_touch_count),
                                                                     notmac_reap_count = sum(notmac_reap_count), notmac_touch_count = sum(notmac_touch_count))


cell_reaper_probs$mac_reaper_ratio <- cell_reaper_probs$mac_reap_count / cell_reaper_probs$mac_touch_count 
cell_reaper_probs$notmac_reaper_ratio <- cell_reaper_probs$notmac_reap_count / cell_reaper_probs$notmac_touch_count
#TODO Once again, this BS needs to be in the preliminary analysis, not here!
for (i in 1:nrow(cell_reaper_probs)) {
  # fix atdc5 mono
  if (grepl('atdc5', cell_reaper_probs[i, 'id_col']) & !grepl('atdc5_oc', cell_reaper_probs[i, 'id_col'])) {
    cell_reaper_probs[i, 'notmac_reap_count'] <-  cell_reaper_probs[i, 'notmac_reap_count'] +  cell_reaper_probs[i, 'mac_reap_count']
    cell_reaper_probs[i, 'mac_reap_count'] <- 0
    cell_reaper_probs[i, 'notmac_touch_count'] <-  cell_reaper_probs[i, 'notmac_touch_count'] +  cell_reaper_probs[i, 'mac_touch_count']
    cell_reaper_probs[i, 'mac_touch_count'] <- 0
    cell_reaper_probs[i, 'notmac_reaper_ratio'] <- cell_reaper_probs[i, 'notmac_reap_count'] / cell_reaper_probs[i, 'notmac_touch_count']
    cell_reaper_probs[i, 'mac_reaper_ratio'] <- 0
  }
  # fix OC mono
  else if (grepl('oc', cell_reaper_probs[i, 'id_col']) & !(grepl('atdc5_oc', cell_reaper_probs[i, 'id_col']) | grepl('callus_oc', cell_reaper_probs[i, 'id_col']))) {
    cell_reaper_probs[i, 'mac_reap_count'] <-  cell_reaper_probs[i, 'notmac_reap_count'] +  cell_reaper_probs[i, 'mac_reap_count']
    cell_reaper_probs[i, 'notmac_reap_count'] <- 0
    cell_reaper_probs[i, 'mac_touch_count'] <-  cell_reaper_probs[i, 'notmac_touch_count'] +  cell_reaper_probs[i, 'mac_touch_count']
    cell_reaper_probs[i, 'notmac_touch_count'] <- 0
    cell_reaper_probs[i, 'mac_reaper_ratio'] <- cell_reaper_probs[i, 'mac_reap_count'] / cell_reaper_probs[i, 'mac_touch_count']
    cell_reaper_probs[i, 'notmac_reaper_ratio'] <- 0
  }
}


cell_reaper_probs_tall <- cell_reaper_probs[c(1,6,7)] %>% pivot_longer(cols = 2:3, names_to = 'reap_type', values_to = 'reaper_prob')

cell_reaper_probs_tall <- cell_reaper_probs_tall %>% mutate(cell_type = case_when(reap_type == 'mac_reaper_ratio' & grepl('atdc5_oc', id_col) ~ 'ATDC5_coculture',
                                                                                          reap_type == 'mac_reaper_ratio' & grepl('atdc5', id_col) ~ 'ATDC5',
                                                                                          reap_type == 'mac_reaper_ratio' & grepl('callus_oc', id_col) ~ 'Callus_coculture',
                                                                                          reap_type == 'mac_reaper_ratio' & grepl('oc', id_col) ~ 'Osteoclasts',
                                                                                          reap_type == 'notmac_reaper_ratio' & grepl('atdc5_oc', id_col) ~ 'ATDC5_coculture',
                                                                                          reap_type == 'notmac_reaper_ratio' & grepl('atdc5', id_col) ~ 'ATDC5',
                                                                                          reap_type == 'notmac_reaper_ratio' & grepl('callus_oc', id_col) ~ 'Callus_coculture',
                                                                                          reap_type == 'notmac_reaper_ratio' & grepl('oc', id_col) ~ 'Osteoclasts',
                                                                                          TRUE ~ 'Other'))


cell_reaper_probs_tall$cell_type <- factor(cell_reaper_probs_tall$cell_type, levels=c('Osteoclasts', 'ATDC5', 'ATDC5_coculture', 'Callus_coculture'))
cell_reaper_probs_tall$reap_type <- factor(cell_reaper_probs_tall$reap_type)
levels(cell_reaper_probs_tall$reap_type) <- c('Not OC Reaper', 'Osteoclast Reaper')

ggplot(data=cell_reaper_probs_tall, aes(cell_type, reaper_prob*100, alpha=reap_type, fill=cell_type)) +
  geom_bar(position='dodge', stat = "summary", fun = "mean") +
  geom_jitter(position = position_jitterdodge(0.1)) +
  ylab('% of cell contacts ') +
  scale_fill_manual(values=pal) +
  scale_alpha_manual(values=c(.66, 1))


for (group_i in 1:length(levels(cell_reaper_probs_tall$cell_type))) {
  temp_df <- cell_reaper_probs_tall %>% filter(cell_type == levels(cell_reaper_probs_tall$cell_type)[group_i])
  print(levels(cell_reaper_probs_tall$cell_type)[group_i])
  print(t.test(reaper_prob ~ reap_type, data=temp_df, paired = TRUE, alternative = "two.sided"))
}


# this result, together with the 'touched not mac before dying' above, show that
# dying cells are equally likely to touch mac or notmac 
# and mac and notmac touch dying cells equally likely, 
# BUT cells that don't touch macrophages are much less likely to die







# histograms
df_childless$mac_t_diff <- apply(df_childless, 1, function(x) {x['mac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_t_diff <- apply(df_childless, 1, function(x) {x['notmac_touch_frames'][[1]] - x['first_apoptosis'][[1]]})
df_childless$notmac_touch_count <- apply(df_childless, 1, function(x) {length(x['notmac_touch_frames'][[1]])})
df_childless$mac_touch_count <- apply(df_childless, 1, function(x) {length(x['mac_touch_frames'][[1]])})

# test2 <- df_childless[df_childless$id_col == unique(df_childless$id_col)[5],]
# test2 <- test2 %>% filter(!is.na(first_apoptosis)) %>% select(first_apoptosis, mac_touch_frames, notmac_touch_frames)

dying_childless <- df_childless %>% filter(!is.na(first_apoptosis)) %>% select(id_col, mac_touch_frames, notmac_touch_frames, mac_t_diff, notmac_t_diff, mac_touch_count, notmac_touch_count)
dying_childless$med_neg_mac_t_diff <- apply(dying_childless, 1, function(x) { median(x['mac_t_diff'][[1]][x['mac_t_diff'][[1]] < 0]) })
dying_childless$mean_neg_mac_t_diff <- apply(dying_childless, 1, function(x) { mean(x['mac_t_diff'][[1]][x['mac_t_diff'][[1]] < 0]) })
dying_childless$med_neg_notmac_t_diff <- apply(dying_childless, 1, function(x) { median(x['notmac_t_diff'][[1]][x['notmac_t_diff'][[1]] < 0]) })
dying_childless$mean_neg_notmac_t_diff <- apply(dying_childless, 1, function(x) { mean(x['notmac_t_diff'][[1]][x['notmac_t_diff'][[1]] < 0]) })

avg_touch_t_tall <- dying_childless[c(1,10,11)] %>% pivot_longer(cols = 2:3, names_to = 'cell_type', values_to = 'avg_touch_t') # mean, not median

ggplot(avg_touch_t_tall) + geom_boxplot(aes(avg_touch_t*8/60, id_col, color = cell_type))
 ggplot() + geom_boxplot(data=avg_touch_t_tall, aes(avg_touch_t*8/60, cell_type))


tall_test <- data.frame()
mac_list <- list()
notmac_list <- list()
for (r in 1:nrow(dying_childless)) {
  if (!is.null(dying_childless[r,'mac_touch_frames'][[1]][[1]]) & length(dying_childless[r,'mac_touch_frames'][[1]][[1]]) != 0) {
    temp_mac <- as.data.frame(append(enframe(unlist(dying_childless[r,'mac_t_diff'][[1]]))[,2],
                                      c(dying_childless[r,'mac_touch_count'], dying_childless[r,'mac_touch_count']+dying_childless[r, 'notmac_touch_count'], dying_childless[r, 'id_col'])))
    mac_list <- append(mac_list, list(temp_mac))
  }
  if (!is.null(dying_childless[r,'notmac_touch_frames'][[1]][[1]]) & length(dying_childless[r,'notmac_touch_frames'][[1]][[1]]) != 0) {
    temp_notmac <- as.data.frame(append(enframe(unlist(dying_childless[r,'notmac_t_diff'][[1]]))[,2],
                                     c(dying_childless[r,'notmac_touch_count'], dying_childless[r,'mac_touch_count']+dying_childless[r, 'notmac_touch_count'], dying_childless[r, 'id_col'])))
    notmac_list <- append(notmac_list, list(temp_notmac))
  }
}
mac_touch_test <- rbindlist(mac_list)
notmac_touch_test <- rbindlist(notmac_list)

for (sample in 1:length(unique(dying_childless$id_col))) {
  histo_plot <- ggplot() + geom_histogram(data=mac_touch_test[mac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = 1/mac_touch_count), binwidth = 0.5, alpha=0.5, fill='red') + 
             geom_density(data=mac_touch_test[mac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = 1/mac_touch_count), color='red') +
             geom_histogram(data=notmac_touch_test[notmac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = 1/notmac_touch_count), binwidth = 0.5, alpha=0.5, fill='black') + 
             geom_density(data=notmac_touch_test[notmac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = 1/notmac_touch_count))
  pdf(paste0("../figures/histogram_", as.character(unique(dying_childless$id_col)[sample]), "_mac_and_notmac.pdf"), width=8, height=5)
  plot(histo_plot)
  dev.off()
}

# setup the weighting for the combined histogram
mac_grand_sum <- sum(mac_touch_test$mac_touch_count)
mac_sample_sum <- mac_touch_test %>% group_by(id_col) %>% summarise(sample_sum = sum(mac_touch_count))
mac_sample_sum$relative_weight <- mac_sample_sum$sample_sum / mac_grand_sum
mac_touch_test <- left_join(mac_touch_test, mac_sample_sum)
notmac_grand_sum <- sum(notmac_touch_test$notmac_touch_count)
notmac_sample_sum <- notmac_touch_test %>% group_by(id_col) %>% summarise(sample_sum = sum(notmac_touch_count))
notmac_sample_sum$relative_weight <- notmac_sample_sum$sample_sum / notmac_grand_sum
notmac_touch_test <- left_join(notmac_touch_test, notmac_sample_sum)

# plot combined histogram
histo_plot <- ggplot() + geom_histogram(data=mac_touch_test, aes(x=value/(60/8), y=..density.., weight = relative_weight/mac_touch_count), binwidth = 0.5, alpha=0, fill='red') + 
  geom_density(data=mac_touch_test, aes(x=value/(60/8), y=..density.., weight = relative_weight/mac_touch_count), color='red') +
  geom_histogram(data=notmac_touch_test[notmac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = relative_weight/notmac_touch_count), binwidth = 0.5, alpha=0, fill='black') + 
  geom_density(data=notmac_touch_test[notmac_touch_test$id_col == unique(dying_childless$id_col)[sample]], aes(x=value/(60/8), y=..density.., weight = relative_weight/notmac_touch_count))
pdf(paste0("../figures/histogram_E2_atdc5_oc_combined_mac_and_notmac.pdf"), width=8, height=5)
plot(histo_plot)
dev.off()


# idea https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
