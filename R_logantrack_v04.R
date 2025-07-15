{library(ggplot2)
library(dplyr)
library(tibble)
library(forcats)
library(ggridges)
library(readr)
library(plotly)
library(MASS)
library(ggpointdensity)
library(ggbeeswarm)
library(magick)
library(viridis)
library(shiny)
library(cluster)    
library(factoextra)
library(DBI)
library(stringr)
library(reshape2)}

#####Import functions#####
source("https://raw.githubusercontent.com/loganslade/LoganTrack/refs/heads/main/tracking_functions.R")
source("https://raw.githubusercontent.com/loganslade/ImgCytometR/refs/heads/main/plotting_functions.R")

#Make color maps#
mycolorv <- viridis(n=20)
mycolorm <- magma(n=20)
mycolorr <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

#####Import and filter data from cell profiler results######

nuclei <- read_csv(file.choose())

col <- colnames(nuclei)
col

#####Config######
channels <- c("zero", "one", "two", "three")
zero <- "pcna"
one<- "cdt1"
two <- "dhb"
three <- "foci"

cyto <- T
cyto_var <- c("cdt1", "dhb", "foci")

#####Re-Name#####

nuclei_2 <- nuclei %>% dplyr::rename(x = "centroid-1", y = "centroid-0", length = "axis_major_length",
                                       

                                       
) %>%
  
  {if("zero" %in% channels)  dplyr::rename(.,
                                           !!paste0("mean_",zero) := "intensity_mean-0",
                                           !!paste0("int_",zero) := "int_C1")
                                          else . 
  } %>% 
  
  {if("one" %in% channels)  dplyr::rename(.,
                                          !!paste0("mean_",one) := "intensity_mean-1",
                                          !!paste0("int_",one) := "int_C2") 
                                          else . 
  } %>% 
  
  {if("two" %in% channels) dplyr::rename(.,
                                         !!paste0("mean_",two) := "intensity_mean-2",
                                         !!paste0("int_",two) := "int_C3") else .
  } %>% 
  
  {if("three" %in% channels) dplyr::rename(.,
                                         !!paste0("mean_",three) := "intensity_mean-3") else .
  }


col2 <- colnames(nuclei_2)
coltest <- data.frame(col,col2)
coltest

nuclei_2 <- nuclei_2 %>% mutate(i_id = paste(n_image, label, sep="_"))
#####Cyto joining####
#Import cyto data and join with nucleus data#
if(cyto == T){
   
  cyto <- read_csv(file.choose()) 
  
  cyto.2 <- cyto %>% 
    
    
    {if("zero" %in% channels)  dplyr::rename(.,
                                             !!paste0("mean_",zero) := "intensity_mean-0")
      else . 
    } %>% 
    
    {if("one" %in% channels)  dplyr::rename(.,
                                            !!paste0("mean_",one) := "intensity_mean-1") 
      else . 
    } %>% 
    
    {if("two" %in% channels) dplyr::rename(.,
                                           !!paste0("mean_",two) := "intensity_mean-2") else .
    } %>% 
    
    {if("three" %in% channels) dplyr::rename(.,
                                             !!paste0("mean_",three) := "intensity_mean-3") else .
    }
  cyto.2 <- cyto.2 %>% mutate(i_id = paste(n_image, label, sep="_")) 
  
  cytonuc <- inner_join(nuclei_2, cyto.2, by = c("i_id" = "i_id"), suffix = c(x="",y="_cyto"))
  
  
for(channel in cyto_var) {
          
      cytonuc <-  mutate(cytonuc, "{channel}_cn" := .data[[paste0("mean_",channel,"_cyto")]]/.data[[paste0("mean_",channel)]])
                          }
  
                          
    }

#####Config 2#####
#Extract Information from images#
#This only works if the time is 3 characters and the position is 2 characters#
cells <- cytonuc %>% mutate(time =  frame,
                            n_position = gsub("XY","",str_extract(image, "(XY[0-9][0-9])")))

#match positions with treatments, if it is a multi treatment experiment, and order the level of the treatments#
import_info <- T
if(import_info == T){
plate <- read_csv(file.choose())
plate$position <- sprintf("%02d", as.numeric(plate$position)) 
cells <- inner_join(cells, plate, by = c("n_position" = "position"))

varNames <- c("treat_1")
varOrder <- c("order_1")

cells <- cells %>% FactorOrder(varNames, varOrder)

for(i in 1:length(varNames)) print(levels(cells[[varNames[i]]])) #Check order is correct

cells$scenetreat <- interaction(cells$position, cells$time)

if(length(varNames) > 1) {cells <- cells %>% combinetreat(varNames)
cells$treat <- as.factor(cells$treat)} else{cells <- cells %>% mutate(treat = .[[varNames[1]]])}
varNames1 <- "treat"

cells <- cells %>% FactorOrder(varNames1, varOrder)
print(levels(cells$treat))
}

rm("nuclei", "nuclei_2", "cyto", "cyto.2", "cytonuc", "coltest")
gc()

######mitosis stats#######
cells <- cells %>%
  mutate(mitosis_score = (perimeter*area*length/mean_pcna), #change name of mean_mturq to mean_pcna depending on nuclear marker
         scale_ms = mitosis_score/(min(mitosis_score)+1)) %>% 
             corany("scale_ms", "position", "quant")



#####Tracking####
library(foreach)
library(doParallel)

#Recommend short test track of either 1 position or 1/4 of the frames#
#To confirm that the jumps and mitosis threshold are set correctly#


#Tracking settings#
ijump <- 40 #Maximum search area for first into second frame, 80 works for 20x 1x1 binning
#mjump <- 40 #Maximum search area for mitotic cells, 80 works for 20x 1x1 binning
sjump <- 23 #Maximum search area for moving cells, 50 works for 20x 1x1 binning
link_jump <- 15 #Set to 30 for 1x1 binned images
SpeedModifier <- 0.9 #Between 0-1. This sets how much the search area is impacted by cell motion 
jitter_correction <- T
if(jitter_correction == T) {jitter_frames <- c(127, 152, 153)} #First frame when plate has moved, Typically after drug addition#
frame_interval <- 10
parallel <- T #Set F if you are only tracking one position

test <- F

mthresh <- quantile(cells$scale_ms, 0.98) #This is the threshold for estimating mitosis in tracking
                                          #0.98 is the percent of cells below the threshold, 
                                          #this may not be universal though, and might need to be changed

positions <- unique(cells$n_position) #Set number of positions to track#

{if(test == T){FrameLimit <- max(as.numeric(cells$time))*0.95}
         else{FrameLimit <- 0}}
maxframe <- max(as.numeric(cells$time)) 

{if(parallel == T){
  numCores <- detectCores()
  cl <- makeCluster((numCores - 4), outfile = "")
  registerDoParallel(cl)
  old <- Sys.time()
  print(old)
  
  btracked <-  foreach (i=positions, .combine=rbind, .packages='dplyr') %dopar% {back_track_hybrid(df = cells, positions = i,
                                                                             i_jump = ijump, m_jump = mjump, s_jump = sjump, smod = SpeedModifier,
                                                                             m_thresh = mthresh, frame_limit = FrameLimit, jitter_correction = jitter_correction,  jitter_frames = jitter_frames)} #Backtrack frame limit for full movie is 0
  
  ftracked <-  foreach (i=positions, .combine=rbind, .packages='dplyr') %dopar% {forward_track_hybrid(df = cells, positions = i,
                                                                                i_jump = ijump, m_jump = mjump, s_jump = sjump, smod = SpeedModifier,
                                                                                m_thresh = mthresh, frame_limit_f = maxframe, jitter_correction = jitter_correction,  jitter_frames = jitter_frames)} #forwardtrack frame limit for full movie is nframe 
  new <- Sys.time() - old
  print(new)
  stopCluster(cl)
  unregister_dopar()}

else{
old <- Sys.time()
btracked <-  foreach (i=positions, .combine=rbind) %do% {back_track_hybrid(df = cells, positions = i,
                                                                                i_jump = ijump, m_jump = mjump, s_jump = sjump, smod = SpeedModifier,
                                                                                m_thresh = mthresh, frame_limit = FrameLimit, jitter_correction = jitter_correction,  jitter_frames = jitter_frames)} #Backtrack frame limit for full movie is 0
ftracked <-  foreach (i=positions, .combine=rbind) %do% {forward_track_hybrid(df = cells, positions = i,
                                                                                   i_jump = ijump, m_jump = mjump, s_jump = sjump, smod = SpeedModifier,
                                                                                   m_thresh = mthresh, frame_limit_f = maxframe, jitter_correction = jitter_correction,  jitter_frames = jitter_frames)} #forwardtrack frame limit for full movie is nframe
new <- Sys.time() - old
print(new)}}

#Merge track ids with cell data and and count tracks#
btracked <- dplyr::inner_join(cells, btracked, by = c("i_id" = "i_id"))
btracks <- unique(btracked$track_id)
last <- cells %>% filter(i_id %in% btracks) %>% mutate(track_id = i_id,
                                                          time.x = time,
                                                          time.y = as.numeric(time),
                                                          position.x = position)
btracked<- bind_rows(btracked,last)

track.count <- btracked %>% group_by(position.x, time.x) %>% summarise(n())
track.count


length.count <- btracked %>% group_by(track_id, position.x) %>% summarise(n()) %>% filter(`n()`>= maxframe)
bctracks <- length.count$track_id
length(bctracks)


ftracked <- dplyr::inner_join(cells, ftracked, by = c("i_id" = "i_id"))
ftracks <- unique(ftracked$track_id)
first <- cells %>% filter(i_id %in% ftracks) %>% mutate(track_id = i_id,
                                                          time.x = time,
                                                          time.y = as.numeric(time),
                                                          position.x = position)
ftracked <- bind_rows(first,ftracked)

track.count <- ftracked %>% group_by(position.x, time.x) %>% summarise(n())
track.count

length.count <- ftracked %>% group_by(track_id, position.x) %>% summarise(n()) %>% filter(`n()`>=  maxframe)
fctracks <- length.count$track_id
length(fctracks)

####Simple track linking####
#Link track intersections for forward and backward tracking#
#Jump is the search area for linking tracks#
#StabilityValue will usually be the biosensor that will be quantified. i.e. timer or DHB or CDC6
#The StabilityValue is used to test similarity when linking tracks# 

{if(parallel == T){
numCores <- detectCores()
cl <- makeCluster(numCores - 4)
registerDoParallel(cl)
old <- Sys.time()
print(old)

btracked_m <- foreach(i=positions, .combine=rbind,  .packages='dplyr') %dopar% 
                    {TrackLinkJitter(positions = i, btracked, ftracked, jump = link_jump, m_thresh = mthresh, StabilityValue = "int_dhb", jitter_correction = jitter_correction,  jitter_frames = jitter_frames)}
new <- Sys.time() - old
print(new)
stopCluster(cl)
unregister_dopar()}
else{btracked_m <- foreach(i=positions, .combine=rbind,  .packages='dplyr') %dopar% 
    {TrackLinkJitter(positions = i, btracked, ftracked, jump = link_jump, m_thresh = mthresh, StabilityValue = "int_dhb", jitter_correction = jitter_correction,  jitter_frames = jitter_frames)}
}
}

#####Export file with the tracks#####
write.csv(btracked_m, "tracks_merged.csv")
btracked_m <- read_csv(file.choose())
btracked_m <- btracked_m[,-1]

length.count <- btracked_m %>% group_by(track_id, position.x) %>% summarise(n()) %>% filter(`n()`>= maxframe)
mer_ctracks <- length.count$track_id
length(mer_ctracks)

####Calculate mitosis score and align####
#Selects tracks that are longer than a set threshold#
LongTrackThreshold <- 250 #number of frames that defines a long track, I usually pick around 24 hours of frames
length.count <- btracked_m %>% group_by(track_id, position.x) %>% summarise(n()) %>% filter(`n()` > LongTrackThreshold)
mer_longtracks <- length.count$track_id
length(mer_longtracks)


#Calculate the first derivative of the mitosis score#
tracked <- btracked_m %>% group_by(track_id) %>% arrange(time.y) %>% mutate(d_mscore = scale_ms - dplyr::lag(scale_ms, n=1))
tracked$d_mscore[is.na(tracked$d_mscore)] <- 0

#Select only long tracks#
longtracked <- tracked %>% filter(track_id %in% mer_longtracks) 

#Detect mitotic events with the DetectMitosis Function for every track#


riseThresh <- quantile(tracked$d_mscore, 0.995)
fallThresh <- quantile(tracked$d_mscore, 0.0077)  

#Get Ids for the mitotic cells#
#This function can be parallelized if needed#
MitoID <- foreach(i=positions, .combine=c,  .packages='dplyr') %do% 
                  {DetectMitosis(df = longtracked, pos=i, rise_thresh=riseThresh, fall_thresh=fallThresh, gap=4, sep=40)}

MitoID <- na.omit(MitoID)

mtrack_label <- longtracked %>% mutate(m_state = case_when(i_id %in% MitoID ~ "M",
                                                      !i_id %in% MitoID ~ "other"))
#Filter out non mitotic tracks#
mtracked <- mtrack_label %>% group_by(track_id) %>% filter("M" %in% m_state)
nmtracked <- mtrack_label %>% group_by(track_id) %>% filter(!"M" %in% m_state)

a <- mtracked %>% filter() %>% group_by(track_id) %>% summarise()
a %>% summarise(n())

#Filter out tracks with continuously high mitosis score#
mtracks.fil <- mtracked %>% group_by(track_id) %>% summarise(count = sum(scale_ms > mthresh)) %>% filter(count < 30) #filter tracks with prolonged mitosis scores#
mtracks.fil_id <- mtracks.fil$track_id
mtracked <- mtracked %>% filter(track_id %in% mtracks.fil_id)

a <- mtracked %>% filter() %>% group_by(track_id) %>% summarise()
a %>% summarise(n())

fmito <- mtracked %>% filter(m_state == "M") %>% group_by(track_id) %>% 
  slice(which.min(time.y)) %>% dplyr::select(track_id, time.y)

mtrackedf <- dplyr::left_join(mtracked, fmito, by = c("track_id" = "track_id"), suffix = c("","omex"))
mtrackedf <- mtrackedf %>% mutate(t_rel_fm = time.y - time.yomex,
                                  h_rel_fm = (t_rel_fm*frame_interval)/60)

lmito <- mtracked %>% filter(m_state == "M") %>% group_by(track_id) %>% 
  slice(which.max(time.y)) %>% dplyr::select(track_id, time.y)

mtrackedf <- dplyr::left_join(mtrackedf, lmito, by = c("track_id" = "track_id"), suffix = c("","lmex"))
mtrackedf <- mtrackedf %>% mutate(t_rel_lm = time.y - time.ylmex,
                                  h_rel_lm = (t_rel_lm*frame_interval)/60)

####Generation assignment####
track_ids <- unique(mtrackedf$track_id)
gens <- list()
i <- 1

#Assign the generation and increase generation every mitosis#
for(track in track_ids){
  track_df <- mtrackedf %>% filter(track_id == track) %>% arrange(time.y)
  gen <- 0
  timer <- 0
  for(n in 1:nrow(track_df)){
    if(track_df$m_state[n] == "M"){
        gen <- gen+1
        timer <- 0}
      if(gen == 0){gens$time_in_cycle[i] <- track_df$h_rel_fm[n]}
      else{gens$time_in_cycle[i] <- timer}
    gens$i_id[i] <- track_df$i_id[n]
    gens$generation[i] <- gen
    i <- i+1
    timer <- timer+1}
}
gensDF <- as.data.frame(gens)  

mtrackedfGens <- dplyr::inner_join(x = mtrackedf, y = gensDF, by = c("i_id" = "i_id"))
mtrackedfGens <- mtrackedfGens %>% mutate(track_id_gen = paste(track_id, generation, sep="_"))

####Plots to check mitosis alignment####

line <- mtrackedfGens %>% filter() %>% 
  ggplot(aes(h_rel_fm, d_mscore, group = track_id))+
  geom_line(linewidth=0.2, alpha = 0.3)+
  labs(x="Time relative 1st mitosis (h)", y="dMscore/dT", title = "", colour="") + 
  scale_y_continuous(transform = "identity")+
  theme_classic()+
  theme(text = element_text(size=18),
        aspect.ratio=0.6,
        legend.title = element_text(size = 12),
        legend.key = element_rect(color = NA),
        axis.text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  geom_hline(yintercept = c(riseThresh, fallThresh), color = "red",linetype = "dashed", linewidth = 0.5)
ggplotly(line)

####Export aligned tracks####
write.csv(mtrackedfGens, "tracks_aligned_annotated.csv")

####Plots for data viz####
bad_tracks <- mtrackedfGens %>% filter(between(h_rel_fm, -3, -1) & dhb_cn < 0.8 | 
                                       between(h_rel_fm, 1.5, 2) & dhb_cn > 0.9) %>% .$track_id
bad_tracks <- c(unique(bad_tracks), "290_305")
length(bad_tracks)
tracks_filt <- mtrackedfGens %>% filter(!track_id %in% bad_tracks)
tracks_filt <- tracks_filt %>% mutate(cdt1_norm = mean_cdt1 - mean_cdt1_cyto)

line <- tracks_filt %>% filter(track_id %in% sample) %>% 
      ggplot(aes(time.y, cdt1_norm, group = track_id))+
  #geom_smooth(aes(color = treat), se =T, span = 0.1, linewidth = 1, alpha = 0.6, method = "loess")+
  #geom_smooth(aes(color = treat), se =F, linewidth = 1, alpha = 0.6, formula = y ~ s(x, bs = "cs", fx = TRUE, k = 40))+
  #geom_line(linewidth=0.2, alpha = 0.3)+
  #geom_point(size = 5, alpha=1, aes(color = time.y))+
  geom_line(linewidth=0.2, alpha = 0.3)+
  #geom_point(data = . %>% filter(m_state == 127), size = 1, alpha=1, color = "Red")+
  #scale_colour_manual(values = c("black","deepskyblue", "red"))+
  labs(x="Time relative 1st mitosis (h)", y="Normalized CDT1", title = "", colour="") + 
  scale_y_continuous(transform = "identity")+
  theme_classic()+
  theme(text = element_text(size=18),
        aspect.ratio=0.6,
        legend.title = element_text(size = 12),
        legend.key = element_rect(color = NA),
        axis.text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_wrap(vars(NULL))
line
ggplotly(line)

a <- tracks_filt %>% filter(treat == "CDK4/6i -Washout") %>% group_by(track_id) 
sample <- sample(unique(a$track_id), 1)
coeff <- 10

line <- tracks_filt %>% filter(track_id %in% sample) %>% 
  ggplot(aes(time.y, group = track_id))+
  #geom_smooth(aes(color = treat), se =T, span = 0.1, linewidth = 1, alpha = 0.6, method = "loess")+
  #geom_smooth(aes(color = bin), se =F, linewidth = 1, alpha = 0.6, formula = y ~ s(x, bs = "cs", fx = TRUE, k = 40))+
  geom_line(aes(y=mean_foci -25), linewidth=1, alpha = 0.3, color = "blue")+
  #geom_line(aes(y=dhb_cn), linewidth=0.2, alpha = 0.3, color = "blue")+
   geom_line(aes(y=cdt1_norm/coeff), linewidth=1, alpha = 0.3, color = "red")+
  #geom_point(size = 5, alpha=1, aes(color = time.y))+
  #geom_line(aes(color = treat), linewidth=0.2, alpha = 0.3)+
  #geom_point(data = . %>% filter(time.y == 127), size = 1, alpha=1, color = "Red")+
  #scale_colour_manual(values = c("black","deepskyblue", "red"))+
  labs(x="Frame", y="PCNA Foci", title = "", colour="") + 
  scale_y_continuous(transform = "identity",
                     sec.axis = sec_axis(~.*coeff, name="Mean CDT1"))+
  theme_classic()+
  theme(text = element_text(size=18),
        aspect.ratio=0.6,
        legend.title = element_text(size = 12),
        legend.key = element_rect(color = NA),
        axis.text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  geom_vline(xintercept = 152)+
  facet_wrap(vars(track_id))
line

#Heatmaps#
library(pheatmap)

heattrack <- tracks_filt %>% filter(between(h_rel_fm, -4, 12), treat == "NT")
length.count <- heattrack %>% group_by(track_id) %>% summarise(n()) %>% filter(`n()` == max(.$`n()`))
heattrack_id <- length.count$track_id
length(heattrack_id)
#sample <- sample(heattrack_id, 300)

heatmat <- heattrack %>% filter(track_id %in% heattrack_id) %>% dplyr::select(h_rel_fm, cdt1_norm, track_id)
heatmat_w <- dcast(heatmat, track_id ~ as.factor(h_rel_fm), value.var = "cdt1_norm")
row.names(heatmat_w) <- heatmat_w$track_id
heatmat_w <- as.matrix(heatmat_w[,-1])

pheatmap(heatmat_w, cluster_cols = F, scale = "none", col = mycolorm, 
         show_rownames = F, show_colnames =F, border_color="NA", na_col = "white")

loghm <- log10(heatmat_w)

pheatmap(loghm, cluster_cols = F, scale = "none", col = mycolorm, 
         show_rownames = F, show_colnames =T, border_color="NA", na_col = "white")


#####Track images export#####
dir <- gsub("\\\\", "/", readClipboard())
dir

save <- gsub("\\\\", "/", readClipboard())
save
setwd(save)

toidf <- tracks_filt %>% filter() 

toidf <- toidf %>% mutate(path = paste(dir,"/",gsub('.{7}$', '', image),"T",frame,"_C1.png",sep=""))
images_seq <- unique(toidf$path)

image <- image_read(toidf$path[1])
image <- image_annotate(image, paste0(toidf$track_id[1]), size = 10, color = "red", degrees = 0, location = paste0("+",(toidf$x[1]/2),"+",(toidf$y[1]/2)))
crop <- image_crop(image, paste0("90x90+",((toidf$x[1]/2)-45),"+",((toidf$y[1]/2)-45)))
crop <- image_scale(crop, "180x180")
crop
a <- 0
for(f in 1:length(images_seq)){
  image <- image_read(toidf$path[f])
  image <- image_annotate(image, paste0(toidf$track_id[f]), size = 10, color = "red", degrees = 0, location = paste0("+",(toidf$x[f]/2),"+",(toidf$y[f]/2)))
  crop <- image_crop(image, paste0("90x90+",((toidf$x[f]/2)-45),"+",((toidf$y[f]/2)-45)))
  crop <- image_scale(crop, "180x180")
  image_write(crop,paste0("T",toidf$track_id[f],"_T",toidf$n_image[f],".png"), compression="Lossless")
  a <- a+1
  if(a == 50){print(a)
              gc()
              a <- 0}
                                }


for(i in 1:length(images_seq)){
  print(images_seq[i])
  filt <- toidf %>% filter(path == images_seq[i])
  image <- image_read(filt$path[1])

  for (o in 1:nrow(filt)){
    image <- image_annotate(image, paste0(filt$track_id[o]), size = 10, color = "red", degrees = 0, location = paste0("+",(filt$x[o]/2),"+",(filt$y[o]/2))) 
                         }
  image_write(image,paste0("T",filt$time.x[o],"_",filt$position.x[o],".png"), compression="Lossless")
  rm(image, filt)
  gc()
  graphics.off() 
}

#####Last Frame object info#####
last_frame <- mtrackedfGens_binned %>% filter(time.x == maxframe) %>% dplyr::select(n_position,image, track_id, i_id, time.y, dhb_cn, mean_cdc6, 
                                                                             mean_cdc6_cyto, mean_foci, x, y, d_mscore, scale_ms, time.yomex, h_rel_fm, h_rel_lm130, h_rel_lm, bin)

last_frame <- tracked %>% filter(time.x == maxframe) %>% dplyr::select(n_position,image, track_id, i_id, time.y, dhb_cn, mean_cdc6, 
                                                                             mean_cdc6_cyto, mean_foci, x, y, d_mscore, scale_ms)

