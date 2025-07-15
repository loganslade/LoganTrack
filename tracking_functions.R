euc_dist <- function(x1, x2){
  return(sqrt(sum((x1 - x2)^2)))}


old_track_link <- function(btracked, ftracked, jump, m_thresh, StabilityValue){
  
  max.jump <- jump
   
  btrack_end <- btracked %>% group_by(track_id) %>% slice(which.min(time.y)) %>% 
    filter(time.y > 1) %>% arrange(desc(time.y))
  
  print(paste0("Match testing ", nrow(btrack_end), " tracks"))
  
  
  btracked$old_track_id <- NA
    
    #For each track end in the position, search for a track in the forward track data#
    
    for(i in btrack_end$track_id){
      btrack <- btrack_end[btrack_end$track_id == i,]
      next_frame <- btrack$time.y - 1
      
      ftracked_search <- ftracked %>% filter(time.y == next_frame,
                                             position.x == btrack$position.x)
      btracked_search <- btracked %>% filter(time.y == next_frame,
                                             position.x == btrack$position.x)
      
      t_match <- ftracked_search %>% filter(between(x, (btrack$x - max.jump), (btrack$x+max.jump)) &
                                              between(y, (btrack$y - max.jump), (btrack$y+max.jump)))
      if(nrow(t_match) == 1){
        
        #Check match similarity#
        if((btrack[StabilityValue]/t_match[StabilityValue] > 2 | btrack[StabilityValue]/t_match[StabilityValue] < 0.6 ) & btrack$scale_ms < m_thresh) {next} 
        
        match_df <- filter(ftracked, track_id %in% t_match$track_id,
                           time.y <= next_frame)
        
        if(any(match_df$i_id %in% btracked$i_id)) {
          btrack_ids <- filter(btracked, i_id %in% match_df$i_id, 
                               time.y <= next_frame)
          match_df <- match_df %>% filter(!i_id %in% btrack_ids$i_id & match_df$time.y > max(btrack_ids$time.y))
          if(nrow(match_df) > 0){
            match_df$old_track_id = match_df$track_id 
            match_df$track_id = i
            btracked <- rbind(btracked, match_df)}
          else{next}
        }                         
        else  {
          match_df$old_track_id = match_df$track_id 
          match_df$track_id = i
          btracked <- rbind(btracked, match_df)
        }
      }
      
    }
  return(btracked) }


unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

DetectMitosis <- function(df, pos, rise_thresh, fall_thresh, gap, sep){
  old <- Sys.time()
  mex_i_id <- c()
  sample <- df %>% filter(n_position == pos)
  tracks <- unique(sample$track_id)
  
  for(track in tracks){
    track_df <- sample %>% group_by(track_id) %>% filter(track_id == track) %>% arrange(time.y)
    
    
    rapid_rises <- which(track_df$d_mscore > rise_thresh)
    rapid_falls <- which(track_df$d_mscore < fall_thresh)
    
    
    paired_rises <- sapply(rapid_rises, function(rise) {
      # Find the first fall within the specified time gap
      fall_after_rise <- rapid_falls[rapid_falls > rise & rapid_falls <= (rise + gap)]
      if (length(fall_after_rise) > 0) return(rise) else {return(NA)}
    })
    
    
    # Filter out unpaired rises
    paired_rises <- na.omit(paired_rises)
    if(length(paired_rises) == 0){next}
    
    paired_falls <- sapply(paired_rises, function(rise) {
      min(rapid_falls[rapid_falls > rise & rapid_falls <= (rise + gap)])})
    
    
    #Filter out doublets
    
    paired_events <- data.frame(
      Rise_Time = paired_rises,
      Fall_Time = paired_falls
    )
    
    filtered_events <- paired_events
    filtered_events <- filtered_events[order(filtered_events$Rise_Time), ] # Sort by rise time
    
    if(nrow(filtered_events) >1){
      for (i in nrow(filtered_events):2) {
        if (filtered_events$Rise_Time[i] - filtered_events$Fall_Time[i - 1] < sep) {
          filtered_events <- filtered_events[-i, ]  # Remove the closer event
        }
      }
    }
    
    mex_i_id <- c(mex_i_id, track_df$i_id[filtered_events$Fall_Time])
  }
  new <- Sys.time() - old
  print(new)
  return(mex_i_id)
}





back_track_multijitter <- function(df, positions, i_jump, m_jump, s_jump, jitter_correction,
                                   jitter_frames, smod, m_thresh, frame_limit){
  df$time <- as.numeric(df$time)
  tracking <- list()
  cur_frame <- max(df$time)
  nframe <- max(df$time)
  x <- 0
  
  for(p in positions){
    old <- Sys.time()
    sample <- df %>% filter(n_position == p)
    
    if(jitter_correction == T){
      
      jitters <- sample %>% filter(time %in% jitter_frames) %>% distinct(time, jitter_x, jitter_y) %>% 
        mutate(post_time = time - 1)
    }
    
    track_set <- sample %>% dplyr::select(c("x","y", "i_id", "n_image", "time", "scale_ms"))
    
    for(a in 1:1){                  
      nobj_frame <- sample %>% filter(time == cur_frame) %>% nrow(.)
      print(paste0("Tracking ",nobj_frame, " objects in position ",p))
      for(b in 1:nobj_frame){
        set <- track_set[track_set$time == cur_frame,]
        cur_object <- set[b,]
        track_id <- cur_object$i_id
        if(length(tracking) > 0 && cur_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == cur_frame]) {next}
        else{
          for(c in (cur_frame-1):frame_limit){
            
            if(jitter_correction == T && c %in% jitters$post_time){shift_x <- jitters$jitter_x[jitters$post_time == c]
            shift_y <- jitters$jitter_y[jitters$post_time == c]}
            
            else{shift_x <- 0 
            shift_y <- 0}
            
            
            if(c == nframe-1){
              max.jump <- i_jump
              
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              
              search <- track_set[track_set$time == c,]
              
              
              next_object <- filter(search, between(search$x,(cur_object$x-shift_x - max.jump), (cur_object$x-shift_x + max.jump)) &
                                      between(search$y,(cur_object$y-shift_y - max.jump), (cur_object$y-shift_y + max.jump))) 
              
              #Track objects without motion for the T1 - T2 frames#
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x-shift_x)-x)^2)+(((cur_object$y-shift_y)-y)^2))) 
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) #make it select only short distances
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x)  
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break} 
            }
            
            #Track objects with motion for all other frames except the T1-T2#
            
            else{
              if(cur_object$scale_ms > m_thresh) {max.jump <- m_jump} else{max.jump <- s_jump}
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              cur_object$enext_x <- (x_cur-shift_x) + x_traj #add enext-j calc hear
              cur_object$enext_y <- (y_cur-shift_y) + y_traj 
              search <- track_set[track_set$time == c,]
              next_object <- filter(search, between(search$x,((cur_object$x - shift_x+(x_traj*smod))- max.jump), ((cur_object$x - shift_x+(x_traj*smod)) + max.jump)) &
                                      between(search$y,((cur_object$y - shift_y+(y_traj*smod)) - max.jump), ((cur_object$y - shift_y + (y_traj*smod)) + max.jump))) 
              
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break} 
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                next_object <- next_object %>% mutate(distance = sqrt(((cur_object$enext_x-x)^2)+((cur_object$enext_y-y)^2))) #removed the enext calcu from hear 
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) 
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break}
              
            }
          }
        }
      }
    }
    new <- Sys.time() - old
    print(new)
  }
  tracks <- as.data.frame(tracking)
  return(tracks)
}


forward_track_multijitter <- function(df, positions, i_jump, m_jump, s_jump, jitter_correction,
                                      jitter_frames, smod, m_thresh, frame_limit_f){
  df$time <- as.numeric(df$time)
  tracking <- list()
  cur_frame <- 1
  nframe <- max(df$time)
  x <- 0
  
  for(p in positions){
    old <- Sys.time()
    sample <- df %>% filter(n_position == p)
    
    if(jitter_correction == T){
      jitters <- sample %>% filter(time %in% jitter_frames) %>% distinct(time, jitter_x, jitter_y)}
    
    track_set <- sample %>% dplyr::select(c("x","y", "i_id", "n_image", "time", "scale_ms"))
    for(a in 1:1){                  
      nobj_frame <- sample %>% filter(time == cur_frame) %>% nrow(.)
      print(paste0("Tracking ",nobj_frame, " objects in position ",p))
      for(b in 1:nobj_frame){
        set <- track_set[track_set$time == cur_frame,]
        cur_object <- set[b,]
        track_id <- cur_object$i_id
        if(length(tracking) > 0 && cur_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == cur_frame]) {next}
        else{
          for(c in (cur_frame+1):frame_limit_f){ 
            
            if(jitter_correction == T && c %in% jitters$time){shift_x <- jitters$jitter_x[jitters$time == c]
            shift_y <- jitters$jitter_y[jitters$time == c]}
            
            else{shift_x <- 0
            shift_y <- 0}
            
            
            if(c == 2){
              max.jump <- i_jump
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              search <- track_set[track_set$time == c,]
              next_object <- filter(search, between(search$x,(cur_object$x+shift_x - max.jump), (cur_object$x+shift_x + max.jump)) &
                                      between(search$y,(cur_object$y+shift_y  - max.jump), (cur_object$y+shift_y  + max.jump))) 
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x) 
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x+shift_x)-x)^2)+(((cur_object$y+shift_y)-y)^2)))
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) #make it select only short distances
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x) 
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break} 
            }
            else{
              if(cur_object$scale_ms > m_thresh) {max.jump <- m_jump} else{max.jump <- s_jump}
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              cur_object$enext_x <- (shift_x+x_cur) + x_traj
              cur_object$enext_y <- (shift_y+y_cur) + y_traj
              search <- track_set[track_set$time == c,]
              next_object <- filter(search, between(search$x,((cur_object$x+shift_x+(x_traj*smod))- max.jump), ((cur_object$x+shift_x+(x_traj*smod)) + max.jump)) &
                                      between(search$y,((cur_object$y+shift_y+(y_traj*smod)) - max.jump), ((cur_object$y+shift_y+(y_traj*smod)) + max.jump))) 
              
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break} 
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x)
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                next_object <- next_object %>% mutate(distance = sqrt(((cur_object$enext_x-x)^2)+((cur_object$enext_y-y)^2)))
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) 
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x)
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break}
              
            }
          }
        }
      }
    }
    new <- Sys.time() - old
    print(new)
  }
  tracks <- as.data.frame(tracking)
  return(tracks)
}



TrackLinkJitter <- function(positions, btracked, ftracked, jump, m_thresh, 
                            StabilityValue, jitter_correction, jitter_frames){
  
  max.jump <- jump
  
  #Filter by position for paralellization#
  ftracked_filt <- ftracked %>% filter(n_position == positions)
  btracked_filt <- btracked %>% filter(n_position == positions)
  
  if(jitter_correction == T){
    
    jitters <- btracked_filt %>% filter(time.y %in% jitter_frames) %>% distinct(time.y, jitter_x, jitter_y) %>% 
      mutate(post_time = time.y - 1)
  }
  
  btrack_end <- btracked_filt %>% group_by(track_id) %>% slice(which.min(time.y)) %>% 
    filter(time.y > 1) %>% arrange(desc(time.y))
  
  
  btracked_filt$old_track_id <- NA
  
  #For each track end in the position, search for a track in the forward track data#
  
  for(i in btrack_end$track_id){
    btrack <- btrack_end[btrack_end$track_id == i,]
    next_frame <- btrack$time.y - 1
    
    {if(jitter_correction == T && next_frame %in% jitters$post_time){shift_x <- jitters$jitter_x[jitters$post_time == next_frame]
    shift_y <- jitters$jitter_y[jitters$post_time == next_frame]}
      
      else{shift_x <- 0 
      shift_y <- 0}}
    
    ftracked_search <- ftracked_filt %>% filter(time.y == next_frame)
    btracked_search <- btracked_filt %>% filter(time.y == next_frame)
    
    t_match <- ftracked_search %>% filter(between(x, (btrack$x - shift_x - max.jump), (btrack$x - shift_x + max.jump)) &
                                            between(y, (btrack$y - shift_y - max.jump), (btrack$y - shift_y + max.jump)))
    if(nrow(t_match) == 1){
      
      #Check match similarity#
      if((btrack[StabilityValue]/t_match[StabilityValue] > 2 | btrack[StabilityValue]/t_match[StabilityValue] < 0.6 ) & btrack$scale_ms < m_thresh) {next} 
      
      match_df <- filter(ftracked_filt, track_id %in% t_match$track_id,
                         time.y <= next_frame)
      
      if(any(match_df$i_id %in% btracked_filt$i_id)) {
        btrack_ids <- filter(btracked_filt, i_id %in% match_df$i_id, 
                             time.y <= next_frame)
        match_df <- match_df %>% filter(!i_id %in% btrack_ids$i_id & match_df$time.y > max(btrack_ids$time.y))
        if(nrow(match_df) > 0){
          match_df$old_track_id = match_df$track_id 
          match_df$track_id = i
          btracked_filt <- rbind(btracked_filt, match_df)}
        else{next}
      }                         
      else  {
        match_df$old_track_id = match_df$track_id 
        match_df$track_id = i
        btracked_filt <- rbind(btracked_filt, match_df)
      }
    }
  }
  return(btracked_filt)
}


back_track_hybrid <- function(df, positions, i_jump, m_jump, s_jump, jitter_correction,
                              jitter_frames, smod, m_thresh, frame_limit){
  df$time <- as.numeric(df$time)
  tracking <- list()
  cur_frame <- max(df$time)
  nframe <- max(df$time)
  x <- 0
  
  for(p in positions){
    old <- Sys.time()
    sample <- df %>% filter(n_position == p)
    
    if(jitter_correction == T){
      
      jitters <- sample %>% filter(time %in% jitter_frames) %>% distinct(time, jitter_x, jitter_y) %>% 
        mutate(post_time = time - 1)
    }
    
    track_set <- sample %>% dplyr::select(c("x","y", "i_id", "flow_x", "flow_y", "n_image", "time", "scale_ms"))
    
    for(a in 1:1){                  
      nobj_frame <- sample %>% filter(time == cur_frame) %>% nrow(.)
      print(paste0("Tracking ",nobj_frame, " objects in position ",p))
      for(b in 1:nobj_frame){
        set <- track_set[track_set$time == cur_frame,]
        cur_object <- set[b,]
        track_id <- cur_object$i_id
        if(length(tracking) > 0 && cur_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == cur_frame]) {next}
        else{
          for(c in (cur_frame-1):frame_limit){
            
            if(jitter_correction == T && c %in% jitters$post_time){shift_x <- jitters$jitter_x[jitters$post_time == c]
            shift_y <- jitters$jitter_y[jitters$post_time == c]}
            
            else{shift_x <- 0 
            shift_y <- 0}
            
            
            if(c == nframe-1){
              max.jump <- i_jump
              
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              
              search <- track_set[track_set$time == c,]
              
              #Track objects without optical for the T1 - T2 frames#
              
              next_object <- filter(search, between(search$x,((cur_object$x + (cur_object$flow_x*smod))- max.jump), ((cur_object$x + (cur_object$flow_x*smod)) + max.jump)) &
                                      between(search$y,((cur_object$y + (cur_object$flow_y*smod)) - max.jump), ((cur_object$y + (cur_object$flow_y*smod)) + max.jump))) 
              
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x+(cur_object$flow_x*smod))-x)^2)+(((cur_object$y+cur_object$flow_y*smod)-y)^2))) 
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) #make it select only short distances
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x)  
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break} 
            }
            
            #Track objects with motion for all other frames except the T1-T2, expect if M#
            
            else{
              max.jump <- s_jump
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              cur_object$enext_x <- (x_cur-shift_x) + x_traj #add enext-j calc hear
              cur_object$enext_y <- (y_cur-shift_y) + y_traj 
              search <- track_set[track_set$time == c,]
              
              if(cur_object$scale_ms > m_thresh){
                
                next_object <- filter(search, between(search$x,((cur_object$x + (cur_object$flow_x*smod))- max.jump), ((cur_object$x + (cur_object$flow_x*smod)) + max.jump)) &
                                        between(search$y,((cur_object$y + (cur_object$flow_y*smod)) - max.jump), ((cur_object$y + (cur_object$flow_y*smod)) + max.jump))) 
              }
              
              else{
                next_object <- filter(search, between(search$x,((cur_object$x - shift_x+(x_traj*smod))- max.jump), ((cur_object$x - shift_x+(x_traj*smod)) + max.jump)) &
                                        between(search$y,((cur_object$y - shift_y+(y_traj*smod)) - max.jump), ((cur_object$y - shift_y + (y_traj*smod)) + max.jump))) 
              }
              
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break} 
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                
                if(cur_object$scale_ms > m_thresh){
                  next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x+(cur_object$flow_x*smod))-x)^2)+(((cur_object$y+cur_object$flow_y*smod)-y)^2)))} 
                
                else{
                  next_object <- next_object %>% mutate(distance = sqrt(((cur_object$enext_x-x)^2)+((cur_object$enext_y-y)^2)))} #removed the enext calcu from hear 
                
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) 
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur-shift_x) 
                y_traj <- y_next-(y_cur-shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break}
              
            }
          }
        }
      }
    }
    new <- Sys.time() - old
    print(new)
  }
  tracks <- as.data.frame(tracking)
  return(tracks)
}

forward_track_hybrid <- function(df, positions, i_jump, m_jump, s_jump, jitter_correction,
                                 jitter_frames, smod, m_thresh, frame_limit_f){
  df$time <- as.numeric(df$time)
  tracking <- list()
  cur_frame <- 1
  nframe <- max(df$time)
  x <- 0
  
  for(p in positions){
    old <- Sys.time()
    sample <- df %>% filter(n_position == p)
    
    sample <- sample %>% mutate(flonext_x = x + flow_x,
                                flonext_y = y + flow_y) 
    
    if(jitter_correction == T){
      jitters <- sample %>% filter(time %in% jitter_frames) %>% distinct(time, jitter_x, jitter_y)}
    
    track_set <- sample %>% dplyr::select(c("x","y", "flonext_x", "flonext_y", "i_id", "n_image", "time", "scale_ms"))
    for(a in 1:1){                  
      nobj_frame <- sample %>% filter(time == cur_frame) %>% nrow(.)
      print(paste0("Tracking ",nobj_frame, " objects in position ",p))
      for(b in 1:nobj_frame){
        set <- track_set[track_set$time == cur_frame,]
        cur_object <- set[b,]
        track_id <- cur_object$i_id
        if(length(tracking) > 0 && cur_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == cur_frame]) {next}
        else{
          for(c in (cur_frame+1):frame_limit_f){ 
            
            if(jitter_correction == T && c %in% jitters$time){shift_x <- jitters$jitter_x[jitters$time == c]
            shift_y <- jitters$jitter_y[jitters$time == c]}
            
            else{shift_x <- 0
            shift_y <- 0}
            
            
            if(c == 2){
              max.jump <- i_jump
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              search <- track_set[track_set$time == c,]
              
              if(cur_object$scale_ms > m_thresh){ 
                next_object <- filter(search, between(search$flonext_x, (cur_object$x - max.jump), (cur_object$x + max.jump)) &
                                        between(search$flonext_y,(cur_object$y - max.jump), (cur_object$y + max.jump))) }
              
              
              else{
                next_object <- filter(search, between(search$x,(cur_object$x+shift_x - max.jump), (cur_object$x+shift_x + max.jump)) &
                                        between(search$y,(cur_object$y+shift_y  - max.jump), (cur_object$y+shift_y  + max.jump))) }
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x) 
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                
                if(cur_object$scale_ms > m_thresh){
                  next_object <- next_object %>% mutate(distance = sqrt(((cur_object$x-flonext_x)^2)+((cur_object$y-flonext_y)^2)))}
                
                else{
                  next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x+shift_x)-x)^2)+(((cur_object$y+shift_y)-y)^2)))}
                
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) #make it select only short distances
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x) 
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break} 
            }
            else{
              max.jump <- s_jump
              x_cur <- cur_object$x
              y_cur <- cur_object$y 
              cur_object$enext_x <- (shift_x+x_cur) + x_traj
              cur_object$enext_y <- (shift_y+y_cur) + y_traj
              search <- track_set[track_set$time == c,]
              
              if(cur_object$scale_ms > m_thresh){ 
                next_object <- filter(search, between(search$flonext_x, (cur_object$x - max.jump), (cur_object$x + max.jump)) &
                                        between(search$flonext_y,(cur_object$y - max.jump), (cur_object$y + max.jump))) }
              
              
              else{
                next_object <- filter(search, between(search$x,(cur_object$x+shift_x - max.jump), (cur_object$x+shift_x + max.jump)) &
                                        between(search$y,(cur_object$y+shift_y  - max.jump), (cur_object$y+shift_y  + max.jump))) }
              
              if(nrow(next_object) == 1){
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break} 
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x)
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else if(nrow(next_object) >1){
                
                if(cur_object$scale_ms > m_thresh){
                  next_object <- next_object %>% mutate(distance = sqrt(((cur_object$x-flonext_x)^2)+((cur_object$y-flonext_y)^2)))}
                
                else{
                  next_object <- next_object %>% mutate(distance = sqrt((((cur_object$x+shift_x)-x)^2)+(((cur_object$y+shift_y)-y)^2)))}
                
                d_min <- min(next_object$distance)
                d_max <- max(next_object$distance)
                d_ratio <- d_max/d_min
                if(d_ratio < 1.5) {break}
                next_object <- next_object %>% slice(which.min(distance)) 
                if(length(tracking) > 0 && next_object$i_id %in% tracking$i_id[tracking$position == p & tracking$time == c]) {break}
                x <- x+1
                tracking$position[x] <- p
                tracking$i_id[x] <-   next_object$i_id
                tracking$time[x] <-   next_object$time
                tracking$parent[x] <- cur_object$i_id
                tracking$track_id[x] <- track_id
                x_next <- next_object$x
                y_next <- next_object$y
                x_traj <- x_next-(x_cur+shift_x)
                y_traj <- y_next-(y_cur+shift_y)
                tracking$x_traj[x] <- x_traj
                tracking$y_traj[x] <- y_traj
                tracking$enext_x[x] <- next_object$x+x_traj
                tracking$enext_y[x] <- next_object$y+y_traj
                cur_object <- next_object}
              else {break}
              
            }
          }
        }
      }
    }
    new <- Sys.time() - old
    print(new)
  }
  tracks <- as.data.frame(tracking)
  return(tracks)
}

