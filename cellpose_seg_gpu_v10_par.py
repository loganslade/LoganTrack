# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:47:09 2026

@author: lgnslade
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:48:45 2024
Update Oct 18 25
Reverted segmentation loop to image-by-image to help memory management
Optimized optical flow measurement for speed and memory
Fixed mask overlay image generation so that it looks somewhat normal
Now processes segmentation at sequentailly by position and parallel with other positions
@author: lslad
"""
import numpy as np
import time, os, sys, gc
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
from scipy.ndimage import zoom
from concurrent.futures import ThreadPoolExecutor, as_completed

#%%

#LOGAN FIX THE FILENAMING SCHEME SO THAT IT MATCHES BETWEEN IMAGE EXPORTS AND CSV COLUMNS!!!!#
#This works with multitiff exports into different positions#
#Setup#
gpu = True
print("GPU is on?: ", gpu)
if gpu == True: import tensorflow as tf
shading_corr = False  
foci = True
qc_image_export = True
image_overlay = True
last_frame_export = True
segmentation = True #Use False to only execute measure measurements
if segmentation == False: print("SEGMENTATION IS OFF")
track_channel = "C1"
channels_list = ["C1", "C2", "C3"] #should match extensions from export tiffs
calculate_jitter = True 
if calculate_jitter == True: 
        frame_of_jitter = [1]
        frame_post_jitter = [x+1 for x in frame_of_jitter]
        
tophat_filter = True

ring_size = 4

if foci == True: foci_channel = "C1"

max_workers = 4


if shading_corr == True:
    shad_cor_ch = ["C1", "C2"]
    shad_cor_paths = ["D://Logan//ffc_images//3x3_ffc_CFP_fused.tif", "D://Logan//ffc_images//3x3_ffc_mCherry_fused.tif"]

test_mode = True

# Limit processing to a subset of frames (frame numbers are 1-based and inclusive).
# Use [start, end] for a range, [start, None] to process from start onward, or
# [None, end] to process through end.  Set limit_frame to False to use all frames.
limit_frame = True
frame_limits = [40, 290]
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
times = [f"{i:02}" for i in range(1, testimg.shape[0]+1)] #If time is less than 3 digits, change 3 to 2

if limit_frame:
    if len(frame_limits) != 2:
        raise ValueError("frame_limits must be [start, end], with None allowed on one side.")
    frame_start, frame_end = frame_limits
    if frame_start is not None and (not isinstance(frame_start, int) or frame_start < 1):
        raise ValueError("frame_limits start must be a positive integer or None.")
    if frame_end is not None and (not isinstance(frame_end, int) or frame_end < 1):
        raise ValueError("frame_limits end must be a positive integer or None.")
    if frame_start is not None and frame_end is not None and frame_start > frame_end:
        raise ValueError("frame_limits start cannot be greater than end.")

    first_frame = 1 if frame_start is None else frame_start
    last_frame = len(times) if frame_end is None else frame_end
    if first_frame > len(times) or last_frame > len(times):
        raise ValueError(f"frame_limits must fall within the stack's 1-{len(times)} frames.")
    frame_indices = list(range(first_frame - 1, last_frame))
else:
    frame_indices = list(range(len(times)))

selected_times = [times[i] for i in frame_indices]
print(f"Processing frames {frame_indices[0] + 1}-{frame_indices[-1] + 1} of {len(times)}.")


#%%    
footprint = np.ones((54, 54))
diameter = 27.5
flow_thersh = 0.6
cellprob = 0.9
sigma = 0.6
    
if test_mode: 
    
    timg = skimage.io.imread(files[3], plugin="tifffile", is_ome=False) #Change files[] to change the image you are testing#
    start_time = time.time()
    img = timg[16] #Change this number to change the image#
    
   
    if shading_corr == True: img = img/sh_dict[track_channel]
    if tophat_filter == True:  
        imgth = skimage.morphology.white_tophat(image = img, footprint=footprint)
    else: 
        imgth = img
    imgsmooth = skimage.filters.gaussian(imgth, sigma=sigma)
    
    
    masks, flows, styles, diams = model.eval(imgsmooth, diameter=diameter, channels=channels, 
                                             flow_threshold= flow_thersh, cellprob_threshold=cellprob) #22,0.6,1.5
    
    
    
    img_ol = np.stack([img] * 3, axis=-1)
    
    outlines = skimage.segmentation.find_boundaries(masks, mode="inner")
    outY, outX = np.nonzero(outlines)
    
    img_ol = img_ol/ np.max(img_ol)
    img_ol[outY, outX] = np.array([1,0,1]) 
    
    img_ol8 = (img_ol*255).astype(np.uint8)
    plt.imshow(img_ol8)
    print(" TIME: ", round(time.time() - start_time, 2))

del(timg)
#%%

#Get file names for all images#

positions = sorted({re.search(r'XY\d{3}', i, re.IGNORECASE).group() for i in file_list}) #Change 2 to 3 if the n position is larger

clist = dict()
cfile_list = dict()
track_file_list = dict()

for position in positions:
    clist[f'{position}_list'] = [x for x in os.listdir(images_dir) if position in x]
    cfile_list[f'{position}'] = [os.path.join(images_dir, file_name) for file_name in clist[f'{position}_list']]
    track_file_list[f'{position}'] = [x for x in clist[f'{position}_list'] if track_channel in x]


#main measurement loop#
if foci == True: fociprint=np.ones((2, 2))
nuclei_list=[]
cyto_list=[]


def process_position(pos_idx, pos):
    print("Start measuring:", pos)
    start_time = time.time()
    local_nuclei_list = []
    local_cyto_list = []

    #Read stack of times at positions pos#
    images = [skimage.io.imread(channel) for channel in cfile_list[pos]]

    flowimg = images[channels_list.index(track_channel)]

    flow_time = time.time()
    flow_by_frame = {}
    for t_image in frame_indices:
        if t_image == 0:
            # Frame 1 has no preceding image, so retain the original zero flow.
            flow_by_frame[t_image] = np.zeros((flowimg.shape[1], flowimg.shape[2], 2), dtype="float16")
        else:
            # The first selected frame is compared with its immediately preceding
            # image, without calculating flow for any unselected frames.
            flow_by_frame[t_image] = cv2.calcOpticalFlowFarneback(
                flowimg[t_image], flowimg[t_image - 1], pyr_scale=0.6, levels=5,
                winsize=8, iterations=1, poly_n=5, poly_sigma=1.1, flags=0, flow=None
            )
    print(pos, "Finished flow in TIME: ", round(time.time() - flow_time, 2))
    del(flowimg)

    #Get images at time = frame, make dictionary#
    for t_image in frame_indices:
        frame = times[t_image]
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
            img_dict["foci"] = skimage.morphology.white_tophat(image=img_dict[foci_channel], footprint=fociprint)

        #Stack the dict images in one np array
        multilist = [value for key, value in img_dict.items()]
        mc = np.stack(multilist, axis=-1)

        #Load the mask image
        mask_file = [x for x in os.listdir(save_dir) if pos in x and f'_T{frame}' in x][0]  #Added underscore to f string to avoid false capture#
        mask_path = os.path.join(save_dir, mask_file)
        mask = io.imread(mask_path)
        mask = skimage.segmentation.clear_border(mask)  #Remove masks touching edge

        #Get ring for cytosol measurements
        c_mask = skimage.segmentation.expand_labels(mask, distance=ring_size)
        c_mask = c_mask - mask  #nc_mask

        #Measure nuc and cyto
        nuclei = skimage.measure.regionprops_table(mask, mc, properties=['area', 'axis_major_length', 'centroid',
                                                                         'intensity_max', 'intensity_mean', 'label', 'perimeter'])

        nuclei['image'] = mask_file[:-9]
        nuclei['n_image'] = (pos_idx * len(times)) + t_image + 1
        nuclei['frame'] = frame
        nuclei['position'] = pos

        x_coords = np.round(nuclei["centroid-1"]).astype(int)
        y_coords = np.round(nuclei["centroid-0"]).astype(int)

        nuclei['flow_x'] = np.round(flow_by_frame[t_image][y_coords, x_coords, 0], 2)
        nuclei['flow_y'] = np.round(flow_by_frame[t_image][y_coords, x_coords, 1], 2)

        for i, ch in enumerate(channels_list):
            nuclei[f'int_{ch}'] = nuclei['area']*nuclei[f'intensity_mean-{i}']

        #If jitter correction is enabled, capture the pre & post jitter images and calculate offset#
        if calculate_jitter == True:
            if pd.to_numeric(frame) in frame_of_jitter:   #When images get loaded, capture in dictionary
                pre_jitter_frame = img_dict[track_channel]
                if t_image + 1 < len(times):
                    post_jitter_frame = images[channels_list.index(track_channel)][t_image+1]
                    jitter = skimage.registration.phase_cross_correlation(post_jitter_frame, pre_jitter_frame, normalization=None)[0]
                    nuclei['jitter_x'] = jitter[1]
                    nuclei['jitter_y'] = jitter[0]
                else:
                    nuclei['jitter_x'] = 0
                    nuclei['jitter_y'] = 0

            else:
                nuclei['jitter_x'] = 0
                nuclei['jitter_y'] = 0

        cyto = skimage.measure.regionprops_table(c_mask, mc, properties=['intensity_mean','area','label'])
        cyto['image'] = mask_file[:-9]
        cyto['n_image'] = (pos_idx * len(times)) + t_image + 1
        cyto['frame'] = frame

        local_nuclei_list.append(nuclei)
        local_cyto_list.append(cyto)
    del(images, flow_by_frame)
    print(pos, "Finished in TIME: ", round(time.time() - start_time, 2))
    return pos_idx, local_nuclei_list, local_cyto_list


#Main segmentation loop for stacks#
#After each position is segmented, start measurement for that position while segmentation continues#

measurement_futures = []

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    for pos_idx, position in enumerate(positions):
        seg_start_time = time.time()
        print(f"Starting segmentation for {position}")
        for track_file in track_file_list[position]:
            if segmentation == True:
                filename = os.path.join(images_dir, track_file)
                timg = skimage.io.imread(filename, plugin="tifffile", is_ome=False)
                for t_image in frame_indices: #Added limiters to stop from processing the blank frames#
                    img = timg[t_image]
                    frame = times[t_image]
                    
                    
                    if shading_corr == True: img = img/sh_dict[track_channel]
                    if tophat_filter == True:  
                        imgth = skimage.morphology.white_tophat(image = img, footprint=footprint)
                    else: 
                        imgth = img
                    imgsmooth = skimage.filters.gaussian(imgth, sigma=sigma)
                    masks, flows, styles, diams = model.eval(imgsmooth, diameter=diameter, channels=channels,
                                                             flow_threshold=flow_thersh, cellprob_threshold=cellprob, batch_size = 16) #22,0.6,1.5
        
                    # save masks as png
                    save_file = f'{track_file[:-4]}_T{frame}_mask.png'
                    save_path = os.path.join(save_dir,save_file)
                    imsave(save_path,masks.astype("uint16"), check_contrast=False)
        
                    #Exports the last frame of each position, used for matching with fixed data later#
                    if t_image == frame_indices[-1]:
                        if last_frame_export == True:
                            save_file = f'{track_file[:-6]}T{frame}_{track_channel}.tiff'
                            save_path_lf = os.path.join(last_frame_dir,save_file)
                            imsave(save_path_lf,img.astype("uint16"), check_contrast=False)
        
                    #Save image with/without segmentation outlines as downsampled png for QC tracks
                    if qc_image_export == True:
                        if image_overlay == True:
                            img_ol = np.stack([imgsmooth] * 3, axis=-1)
        
                            outlines = skimage.segmentation.find_boundaries(masks, mode="inner")
                            outY, outX = np.nonzero(outlines)
        
                            img_ol = img_ol/ np.max(img_ol)
                            img_ol[outY, outX] = np.array([1,0,1])
        
                            img_ol8ds = skimage.transform.rescale(img_ol, 0.5, anti_aliasing=True, channel_axis=2)
                            img_ol8ds = (img_ol8ds*255).astype(np.uint8)
        
                            save_file = f'{track_file[:-6]}T{frame}_{track_channel}.png'
                            save_path_qc = os.path.join(img_save_dir,save_file)
                            imsave(save_path_qc,img_ol8ds, check_contrast=False)
                        else:
                            img_ol = img_ol/ np.max(img_ol)
                            img_ol8 = img_ol*255
        
                            save_file = f'{track_file[:-6]}T{frame}_{track_channel}.png'
                            save_path = os.path.join(img_save_dir,save_file)
                            imsave(save_path,img_ol8.astype("uint8"), check_contrast=False)
        
                    
                del(timg)
        print(position, "Finished segmentation in TIME: ", round(time.time() - seg_start_time, 2))
        measurement_futures.append(executor.submit(process_position, pos_idx, position))

    position_results = []
    for future in as_completed(measurement_futures):
        position_results.append(future.result())

for _, pos_nuclei, pos_cyto in sorted(position_results, key=lambda x: x[0]):
    nuclei_list.extend(pos_nuclei)
    cyto_list.extend(pos_cyto)


dfnuc = pd.concat([pd.DataFrame(l) for l in nuclei_list],axis=0)
dfcyto = pd.concat([pd.DataFrame(l) for l in cyto_list],axis=0)

dfnuc.to_csv("nuclei_cp2.csv", index=False)
dfcyto.to_csv("cyto_cp2.csv", index=False)
