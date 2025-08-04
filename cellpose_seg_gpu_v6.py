# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:48:45 2024
Update Thur July 31 25
Reverted segmentation loop to image-by-image to help memory management
Optimized optical flow measurement for speed and memory
@author: lslad
"""
import numpy as np
import time, os, sys
from urllib.parse import urlparse
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import cellpose
from cellpose import utils, io, plot
import skimage
from skimage.io import imsave
import pandas as pd
import re
from pathlib import Path
import cv2
#%%

#LOGAN FIX THE FILENAMING SCHEME SO THAT IT MATCHES BETWEEN IMAGE EXPORTS AND CSV COLUMNS!!!!#
#This works with multitiff exports into different positions#
#Setup#
gpu = True
print("GPU is on?: ", gpu)
if gpu == True: import tensorflow as tf
shading_corr = False  
foci = False
qc_image_export = True
image_overlay = True
last_frame_export = True
segmentation = True #Use False to only export images 
track_channel = "C1"
channels_list = ["C1", "C2"] #should match extensions from export tiffs
calculate_jitter = True 
if calculate_jitter == True: 
        frame_of_jitter = [1, 336]
        frame_post_jitter = [x+1 for x in frame_of_jitter]


if foci == True: foci_channel = "C1"



if shading_corr == True:
    shad_cor_ch = ["C1", "C2"]
    shad_cor_paths = ["D://Logan//ffc_images//3x3_ffc_CFP_fused.tif", "D://Logan//ffc_images//3x3_ffc_mCherry_fused.tif"]

test_mode = True
#%%

#Generate dictionary of shading correction images#
if shading_corr == True:
    sh_dict = dict() 
    for ch, path in zip(shad_cor_ch, shad_cor_paths):  
        sh_dict[f'{ch}'] = io.imread(path)

print(os.getcwd())

#Get the list of files in the directory#
images_dir = os.getcwd()

file_list = [x for x in os.listdir(images_dir) if track_channel in x]

files = [] 
for file_name in file_list: files.append(os.path.join(images_dir, file_name)) 
    
#Set the size of the background rolling square (should be about the diameter of a mitotic cell)#
footprint = np.ones((30, 30))

#Cellpose setup#
from cellpose import models, io

os.environ['CUDA_VISIBLE_DEVICES'] = '0'
cellpose.core.use_gpu(gpu_number=0, use_torch=True)
model = models.Cellpose(gpu=True, model_type='cyto')
channels = [[0,0]]


#Create save directory, if they don't exist
#Create lists of paths for image exports
save_dir = os.path.join(images_dir,'masks')
Path(save_dir).mkdir(parents=True, exist_ok=True)

if qc_image_export == True: 
    img_save_dir = os.path.join(images_dir,'qc_images')
    Path(img_save_dir).mkdir(parents=True, exist_ok=True)
if last_frame_export == True: 
    last_frame_dir = os.path.join(images_dir,'last_frame') 
    Path(last_frame_dir).mkdir(parents=True, exist_ok=True)

#Import a single image and get dimensions and generate the timepoint names#
testimg =  skimage.io.imread(files[0], plugin="tifffile", is_ome=False)
times = [f"{i:03}" for i in range(1, testimg.shape[0]+1)] #If time is less than 3 digits, change 3 to 2

#%%
#Test segmentation parameters on a small number of images before segmenting the whole thing#
footprint = np.ones((36, 36))
diameter = 25
flow_thersh = 0.7
cellprob = 1


if test_mode: 

    timg = skimage.io.imread(files[0], plugin="tifffile", is_ome=False) #Change files[] to change the image you are testing#
    
    img = timg[108] #Change this number to change the image#
    
    
    
    if shading_corr == True: img = img/sh_dict[track_channel]
    #img = skimage.exposure.rescale_intensity(img, in_range=(1000, 3500)) #make sure to verify in imageJ with a few images#
    imgth = skimage.morphology.white_tophat(image = img, footprint=footprint)
    imgsmooth = skimage.filters.gaussian(imgth, sigma=0.8)
    
    
    masks, flows, styles, diams = model.eval(imgsmooth, diameter=diameter, channels=channels, 
                                             flow_threshold= flow_thersh, cellprob_threshold=cellprob) #22,0.6,1.5
    
    
    
    img_ol = np.stack([img] * 3, axis=-1)
    
    outlines = skimage.segmentation.find_boundaries(masks, mode="inner")
    outY, outX = np.nonzero(outlines)
    
    img_ol = img_ol/ np.max(img_ol)
    img_ol[outY, outX] = np.array([1,0,1]) 
    
    img_ol8 = (img_ol*255).astype(np.uint8)
    plt.imshow(img_ol8)

del(timg)
#%%
#Main segmentation loop for stacks# 
    
#Seg loop for single images#

for filename, filelist in zip(files, file_list):
    timg = skimage.io.imread(filename, plugin="tifffile", is_ome=False)
    start_time = time.time()
    for img, frame in zip(timg,times): #Added limiters to stop from processing the blank frames#
        start_time = time.time()
        
        if segmentation == True: 
        
            if shading_corr == True: img = img/sh_dict[track_channel]
            img = skimage.morphology.white_tophat(image = img, footprint=footprint)
            imgsmooth = skimage.filters.gaussian(img, sigma=0.8)
            masks, flows, styles, diams = model.eval(imgsmooth, diameter=diameter, channels=channels, 
                                                     flow_threshold=flow_thersh, cellprob_threshold=cellprob, batch_size = 16) #22,0.6,1.5
        
            # save masks as png
            save_file = f'{filelist[:-4]}_T{frame}_mask.png'
            save_path = os.path.join(save_dir,save_file)
            imsave(save_path,masks.astype("uint16"), check_contrast=False)
        
        #Exports the last frame of each position, used for matching with fixed data later#
        
        if frame == max(times): 
            if last_frame_export == True: 
                save_file = f'{filelist[:-6]}T{frame}_{track_channel}.tiff'
                save_path = os.path.join(last_frame_dir,save_file)
                imsave(save_path,img.astype("uint16"), check_contrast=False)
        
        #Save image with/without segmentation outlines as downsampled png for QC tracks
        
        if qc_image_export == True: 
        
            img_ol = skimage.transform.rescale(img, 0.5, anti_aliasing=True)
            
            if image_overlay == True:
            
                img_ol = np.stack([img_ol] * 3, axis=-1)
                masks_ds = skimage.transform.rescale(masks, 0.5, anti_aliasing=False)
                
                outlines = skimage.segmentation.find_boundaries(masks_ds, mode="inner")
                outY, outX = np.nonzero(outlines)
                
                img_ol = img_ol/ np.max(img_ol)
                img_ol[outY, outX] = np.array([1,0,1]) 
                
                img_ol8 = img_ol*255
                
                save_file = f'{filelist[:-6]}T{frame}_{track_channel}.png'
                save_path = os.path.join(img_save_dir,save_file)
                imsave(save_path,img_ol8.astype("uint8"), check_contrast=False)
            else: 
                img_ol = img_ol/ np.max(img_ol)
                img_ol8 = img_ol*255
                
                save_file = f'{filelist[:-6]}T{frame}_{track_channel}.png'
                save_path = os.path.join(img_save_dir,save_file)
                imsave(save_path,img_ol8.astype("uint8"), check_contrast=False)
        

            
        
        print("Image: ",filename, "Frame: ",frame, " TIME: ", round(time.time() - start_time, 2))

#%%
#Get file names for all images# 

 
base_name = [x[:-9] for x in os.listdir(save_dir)]
positions = [re.search(r'XY\d{2}', i, re.IGNORECASE).group() for i in file_list] #Change 2 to 3 if the n position is larger

clist = dict()
cfile_list = dict()

for position in positions:
    clist[f'{position}_list'] = [x for x in os.listdir(images_dir) if position in x]
    cfile_list[f'{position}'] = [os.path.join(images_dir, file_name) for file_name in clist[f'{position}_list']]


#main measurement loop#
if foci == True: fociprint=np.ones((2, 2))
image=1
nuclei_list=[]
cyto_list=[]


images=[]

for pos in positions:  

    start_time = time.time()
    if calculate_jitter == True: jitter_dict = dict()
      
    #Read stack of times at positions pos#
    
    images = [skimage.io.imread(channel, plugin="tifffile", is_ome=False) for channel in cfile_list[pos]]
    
    flowimg = images[channels_list.index(track_channel)]
    
    flowstack = np.zeros((testimg.shape[0], testimg.shape[1],  testimg.shape[2], 2), dtype = "float16")
    
    flow_time = time.time()
    for i in range(flowimg.shape[0]):   
        if i == 0: pre = flowimg[i] 
        else: 
            post = flowimg[i]
            flowstack[i] = cv2.calcOpticalFlowFarneback(post, pre, pyr_scale=0.6, levels=5, winsize=8, iterations=1, poly_n=5, poly_sigma=1.1, flags=0, flow=None)
            pre = post
    print(pos, "Finished flow in TIME: ", round(time.time() - flow_time, 2))
    del(flowimg, pre, post)
    
    #Get images at time = frame, make dictionary#
    for frame, t_image in zip(times, range(0,(np.shape(images[0])[0]))): 
        image_set = [image[t_image] for image in images]
        img_dict = dict()
        for chan, im in zip(channels_list, image_set):
            img_dict[f'{chan}'] = im
            
        #Get shading correction images and correct#
        if shading_corr == True:
            for sh_ch in shad_cor_ch: 
                uncor_im = img_dict[f'{sh_ch}']
                shading_im = sh_dict[f'{sh_ch}'] 
                cor_im = uncor_im/shading_im
                cor_im = cor_im.astype(dtype="uint16")
                img_dict[f'{sh_ch}'] = cor_im
        
        #make the foci image, add to dict
        if foci == True:
            img_dict["foci"] = skimage.morphology.white_tophat(image = img_dict[foci_channel], footprint=fociprint)
        
        #Stack the dict images in one np array 
        multilist = [value for key,value in img_dict.items()]
        mc = np.stack(multilist,  axis=-1)
        
        #Load the mask image
        mask_file = [x for x in os.listdir(save_dir) if pos in x and f'_T{frame}' in x][0]  #Added underscore to f string to avoid false capture#
        mask_path = os.path.join(save_dir, mask_file)
        mask = io.imread(mask_path) 
        mask = skimage.segmentation.clear_border(mask)  #Remove masks touching edge 
        
        #Get ring for cytosol measurements
        nc_mask = skimage.segmentation.expand_labels(mask, distance = 1) 
        c_mask =skimage.segmentation.expand_labels(mask, distance = 4) 
        c_mask = c_mask -  mask #nc_mask
        
        #Measure nuc and cyto
        nuclei = skimage.measure.regionprops_table(mask, mc, properties=['area', 'axis_major_length', 'centroid', 'eccentricity', 
                                                                            'intensity_max', 'intensity_mean', 'label', 'orientation',
                                                                            'perimeter', 'solidity'])
        
        nuclei['image'] = mask_file[:-9]
        nuclei['n_image']=image
        nuclei['frame']=frame
        nuclei['position']=pos
        
        x_coords = np.round(nuclei["centroid-1"]).astype(int)
        y_coords = np.round(nuclei["centroid-0"]).astype(int)
        
        nuclei['flow_x'] = np.round(flowstack[t_image,y_coords, x_coords, 0], 2) 
        nuclei['flow_y'] = np.round(flowstack[t_image,y_coords, x_coords, 1], 2) 
        
        for i,ch in enumerate(channels_list): 
            nuclei[f'int_{ch}']= nuclei['area']*nuclei[f'intensity_mean-{i}']
            
        #If jitter correction is enabled, capture the pre & post jitter images and calculate offset#
        if calculate_jitter == True: 
            
            if pd.to_numeric(frame) in frame_of_jitter or pd.to_numeric(frame) in frame_post_jitter:   #When images get loaded, capture in dictionary
                jitter_dict[pd.to_numeric(frame)] = img_dict[track_channel]
             

            if pd.to_numeric(frame) in frame_post_jitter: #When frame is post jitter, run correction for frame and prior frame
                
                jitter = skimage.registration.phase_cross_correlation(jitter_dict[pd.to_numeric(frame)], jitter_dict[pd.to_numeric(frame)-1], normalization=None)[0]
                nuclei['jitter_x'] = jitter[1]
                nuclei['jitter_y'] = jitter[0]

            else: 
                nuclei['jitter_x'] = 0
                nuclei['jitter_y'] = 0
        
        cyto = skimage.measure.regionprops_table(c_mask, mc, properties=['intensity_mean', 'label'])
        cyto['image'] = mask_file[:-9]
        cyto['n_image']=image
        cyto['frame']=frame
        
        nuclei_list.append(nuclei)
        cyto_list.append(cyto)
        image=image+1
    print(pos, "Finished in TIME: ", round(time.time() - start_time, 2))


dfnuc = pd.concat([pd.DataFrame(l) for l in nuclei_list],axis=0)
dfcyto = pd.concat([pd.DataFrame(l) for l in cyto_list],axis=0)

dfnuc.to_csv("nuclei.csv", index=False)
dfcyto.to_csv("cyto.csv", index=False)

