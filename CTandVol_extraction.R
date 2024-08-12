library(ANTsR)
library(ANTsRCore)
library(MNITemplate)

#setwd('~/Desktop/ADNI-2-all/BIDS/')
setwd('ADNI2/')

# Cortical Thickness Extraction -------------------------------------------
setwd('ADNI2/Unprocessed Images/') # BIDS input directory

### ***NOTE: If you have already run fMRIPrep, skip to Step 4.

# Step 1: Register T1w image to template and skull-strip
T1_Template = antsImageRead(MNITemplate::getMNIPath(what = "T1", res = "2mm"))
Brain_Template = antsImageRead(MNITemplate::getMNIPath(what = "Brain", res = "2mm"))
Mask_Template = antsImageRead(MNITemplate::getMNIPath(what = "Brain_Mask", res = "2mm"))
plot(Brain_Template)

# note: using BIDS input directory
T1_Img = antsImageRead('sub-ID0/ses-01/anat/sub-ID0_ses-01_T1w.nii.gz')
T1_Reg = antsRegistration(fixed = T1_Template, moving = T1_Img, typeofTransform = "SyN", verbose = T)
Brain_Img_Reg = antsImageClone(T1_Reg$warpedmovout); Brain_Img_Reg[Mask_Template==0] = 0
plot(Brain_Img_Reg)

# Step 2: Get mask of registered image and apply N4 bias field correction
Mask_Img_Reg = getMask(Brain_Img_Reg)
Brain_Img_Reg_BiasCorrected = n4BiasFieldCorrection(Brain_Img_Reg)

# Step 2B. Optional downsample to save time
Brain_Img_Reg_DS = resampleImage(Brain_Img_Reg_BiasCorrected, resampleParams = c(48,48,48), useVoxels = T)

# Step 3: Perform segmentation and save segmentation image results

## 1 - Cerebrospinal fluid (CSF)
## 2 - Gray Matter (GM)
## 3 - White Matter (WM)
SegResults_kmeans = kmeansSegmentation(img = Brain_Img_Reg_BiasCorrected, k = 3, kmask = Mask_Img_Reg, verbose = T)
SegImg = SegResults_kmeans$segmentation
CSFProb = SegResults_kmeans$probabilityimages[[1]]
GMProb = SegResults_kmeans$probabilityimages[[2]]
WMProb = SegResults_kmeans$probabilityimages[[3]]

# Step 4: Run cortical thickness extraction using "kellyKapowski" function in ANTsR
## REQUIRED INPUTS:
## [Note: all of these outputs should be saved to under the anat/ folder after running fMRIPrep]
## [Note: make sure to use the registered outputs of fMRIPrep, which include labels like "...space-MNI152NLin2009cAsym..." in the file names]
## - Tissue Segmentation Image: 0-background, 1-CSF, 2-GM, 3-WM ("...dseg...")
## - GM Probability Map Image ("...GM_probseg...")
## - WM Probability Map Image ("...WM_probseg...")

# change to fMRIPrep outputs directory
setwd('ADNI2/Processed Images/')
ds = T
pds = rep(48,3)

# downsampling optional but currently this code is configured to downsampled images with 48x48x48 resolution 
if (ds == F) {
  t1 = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz')
  mask = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz')
  dseg = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_dseg.nii.gz')
  gm_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz')
  wm_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz')
  csf_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz')
} else {
  t1 = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
  mask = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
  dseg = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_dseg.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
  gm_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
  wm_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
  csf_prob = antsImageRead('sub-ID0_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz') %>% resampleImage(resampleParams = pds, useVoxels = T)
}

### Cortical thickness extraction using kellyKapowski
## Note: r and m are tuning parameters which affect the optimization algorithm (gradient descent).
## Note: Typical ranges of r and m may depend on image resolution. For 48x48x48, I find best performance with r around 1-2 and m around 0.5-1.
## Note: Expected mean/median CT lies somewhere in the range of 1-4, though it could differ slightly.
## Note: If you see mean/median CT much outside of range, it indicates need to change the r/m parameters.
## Note: Unfortunately, there is not a way to set the seed for this function in R. Multiple replicates per image may be advised.
## Note: Lastly, you may find instances where the code returns very large values for the CT. If this values take up a large portion of the image, you may want to adjust r/m parameters and/or re-run.

# CT_Img = kellyKapowski(s = dseg, g = gm_prob, w = wm_prob, r = 2, m = .5, verbose = T)
CT_Img = kellyKapowski(s = dseg, g = gm_prob, w = wm_prob, r = 0.025, m = 1.5, verbose = T)
median(CT_Img[CT_Img>0])
mean(CT_Img[CT_Img>0])
mean(CT_Img[CT_Img>0 & CT_Img < 10])
hist(CT_Img[CT_Img>0])
hist(CT_Img[CT_Img>0 & CT_Img < 10])
plot(CT_Img)


CT_Img2 = kellyKapowski(s = dseg, g = gm_prob, w = wm_prob, r = 2, m = .5, verbose = T)
median(CT_Img2[CT_Img2>0])
mean(CT_Img2[CT_Img2>0])
mean(CT_Img2[CT_Img2>0 & CT_Img2 < 10])
hist(CT_Img2[CT_Img2>0])
hist(CT_Img2[CT_Img2>0 & CT_Img2 < 10])
plot(CT_Img2)

# Voxel-level Tissue Volumes ----------------------------------------------

## We may also want to extract the estimated amounts of each tissue type (including GM and WM)
## within each voxel, following preprocessing.

## To do so, we can either start by running Steps 1-3 as in the CT extraction example above, or run 
## fMRIPrep. This will give us the registered image, segmentation labels, and tissue probabilities.

## **The key inputs now are the segmentation probability images ("GM_probseg", etc.).
## We can simply multiple these voxel-level probabilities by the volume of a given voxel to find
## estimated volumes for GM, WM, and CSF

## Note: I am using fMRIPrep outputs as above
vox_vol = prod(antsGetSpacing(t1)) # in mm
GM_vol = gm_prob*vox_vol
WM_vol = wm_prob*vox_vol
hist(GM_vol[GM_vol > 0])
hist(WM_vol[WM_vol > 0])
plot(GM_vol)
plot(WM_vol)


