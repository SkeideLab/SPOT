#!/bin/bash
set -u -x -e
# Adapted from https://git.fmrib.ox.ac.uk/seanf/asymmetry-analysis/-/blob/master/hcp_surface.sh?ref_type=heads
# This script maps dHCP volumetric fMRI data to the surface
# It is largely a vanilla implementation of the HCP Surface Pipeline
# One of the main changes is that the fMRI is mapped in the fMRI space
# The script expects the data to be in the standard (non-bids) directory structure of the dHCP fMRI pipeline

# --------------------------------------------

# As this is a derivative of the HCP Surface Pipeline, the terms of the HCP license also apply:

# Human Connectome Project Pipelines = THIS SOFTWARE
# Copyright (c) 2011-2021 The Human Connectome Project and The Connectome Coordination Facility
# Redistribution and use in source and binary forms, with or without modification,
# is permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice,
# this list of conditions, and the following disclaimer.

# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions, and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# The names of Washington University in St. Louis, the University of Minnesota,
# Oxford University, the Human Connectome Project, or any contributors
# to this software may not be used to endorse or promote products derived
# from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# --------------------------------------------
# NOTE the second release of dhcp uses dhcp40wk, not extdhcp40wk as template
subid=$1
sesid=$2
path_wbcommand=$3
anat_dir=$4
func_dir=$5
work_dir=$6
WORKBENCHDIR=$(dirname $path_wbcommand)
export WORKBENCHDIR

# inputs
sub_anat="${anat_dir}/sub-${subid}/ses-${sesid}"
sub_func="${func_dir}/sub-${subid}/ses-${sesid}"
func="${sub_func}/func/sub-${subid}_ses-${sesid}_task-rest_desc-preproc_bold.nii.gz"

# warps to templates
struct2func_xfm="${sub_func}/xfm/sub-${subid}_ses-${sesid}_from-T2w_to-bold_mode-image.mat"

#-------------------- outputs -----------------------------------------------#

output_dir="${work_dir}/sub-${subid}/ses-${sesid}"
mkdir -p "${output_dir}/anat" "${output_dir}/work" "${output_dir}/func"
ref_volume="${sub_func}/func/sub-${subid}_ses-${sesid}_task-rest_desc-firstvol_bold.nii.gz"
func_on_surf="${output_dir}/func/sub-${subid}_ses-${sesid}_hemi-{hemi_upper}_mesh-native_bold.func.gii"
#-----------------parameters------------------------------#
NeighborhoodSmoothing="5"
Factor="0.5"
LeftGreyRibbonValue="1"
RightGreyRibbonValue="1"

###############
if [[ -f "${func_on_surf/\{hemi_upper\}/L}" && -f "${func_on_surf/\{hemi_upper\}/R}" ]]; then
  echo "Functional data on surface already exist at ${func_on_surf}. Skipping..."
  exit 0
fi
if [ ! -f "$ref_volume" ]; then
  #create nifti with only first volume to serve as ref_volume
  fslroi "${func}" \
    "$ref_volume" \
    0 1
fi

# ################
# warp native surface from struct space to func space
# adapted from PostFreeSurfer/scripts/FreeSurfer2CaretConvertAndRegisterNonlinear.sh
# ################

for hemi in left right; do
  if [ $hemi = "left" ]; then
    GreyRibbonValue="$LeftGreyRibbonValue"
    hemi_upper="L"
  elif [ $hemi = "right" ]; then
    GreyRibbonValue="$RightGreyRibbonValue"
    hemi_upper="R"
  fi

  # move all surfaces into func space
  for surface in wm pial midthickness inflated veryinflated; do
    # transform native anatomical surfaces: T2w space --> func space
    $path_wbcommand -surface-apply-affine \
      "${sub_anat}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_space-T2w_${surface}.surf.gii" \
      "${struct2func_xfm}" \
      "${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_${surface}.surf.gii" \
      -flirt \
      "${sub_anat}/anat/sub-${subid}_ses-${sesid}_desc-restore_T2w.nii.gz" \
      "$ref_volume"
  done

  # ################
  # ribbon constrained mapping
  # adapted from fMRISurface/scripts/RibbonVolumeToSurfacemapping.sh
  # ################

  $path_wbcommand -create-signed-distance-volume \
    "${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_wm.surf.gii" \
    "${ref_volume}" \
    "${output_dir}/work/${hemi}.white.native.nii.gz"

  $path_wbcommand -create-signed-distance-volume \
    "${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_pial.surf.gii" \
    "${ref_volume}" \
    "${output_dir}/work/${hemi}.pial.native.nii.gz"

  # create ribbon mask between white and pial mask and fill it with GreyRibbonValue
  fslmaths "${output_dir}/work/${hemi}.white.native.nii.gz" -thr 0 -bin -mul 255 "${output_dir}/work/${hemi}.white_thr0.native.nii.gz"
  fslmaths "${output_dir}/work/${hemi}.white_thr0.native.nii.gz" -bin "${output_dir}/work/${hemi}.white_thr0.native.nii.gz"
  fslmaths "${output_dir}/work/${hemi}.pial.native.nii.gz" -uthr 0 -abs -bin -mul 255 "${output_dir}/work/${hemi}.pial_uthr0.native.nii.gz"
  fslmaths "${output_dir}/work/${hemi}.pial_uthr0.native.nii.gz" -bin "${output_dir}/work/${hemi}.pial_uthr0.native.nii.gz"
  fslmaths "${output_dir}/work/${hemi}.pial_uthr0.native.nii.gz" -mas "${output_dir}/work/${hemi}.white_thr0.native.nii.gz" -mul 255 "${output_dir}/work/${hemi}.ribbon.nii.gz"
  fslmaths "${output_dir}/work/${hemi}.ribbon.nii.gz" -bin -mul "$GreyRibbonValue" "${output_dir}/work/${hemi}.ribbon.nii.gz"
  rm "${output_dir}/work/${hemi}.white.native.nii.gz" "${output_dir}/work/${hemi}.white_thr0.native.nii.gz" "${output_dir}/work/${hemi}.pial.native.nii.gz" "${output_dir}/work/${hemi}.pial_uthr0.native.nii.gz"
done

# Combine ribbon masks for both hemispheres
fslmaths ${output_dir}/work/left.ribbon.nii.gz -add ${output_dir}/work/right.ribbon.nii.gz ${output_dir}/work/ribbon_only.nii.gz
rm ${output_dir}/work/left.ribbon.nii.gz ${output_dir}/work/right.ribbon.nii.gz

# compute mean and standard deviation of functional image
fslmaths ${func} -Tmean ${output_dir}/work/mean -odt float
fslmaths ${func} -Tstd ${output_dir}/work/std -odt float
fslmaths ${output_dir}/work/std -div ${output_dir}/work/mean ${output_dir}/work/cov
fslmaths ${output_dir}/work/cov -mas ${output_dir}/work/ribbon_only.nii.gz ${output_dir}/work/cov_ribbon

#TODO what and why is being smoothed here?
fslmaths ${output_dir}/work/cov_ribbon -div $(fslstats ${output_dir}/work/cov_ribbon -M) ${output_dir}/work/cov_ribbon_norm
fslmaths ${output_dir}/work/cov_ribbon_norm -bin -s $NeighborhoodSmoothing ${output_dir}/work/SmoothNorm
fslmaths ${output_dir}/work/cov_ribbon_norm -s $NeighborhoodSmoothing -div ${output_dir}/work/SmoothNorm -dilD ${output_dir}/work/cov_ribbon_norm_s$NeighborhoodSmoothing
fslmaths ${output_dir}/work/cov -div $(fslstats ${output_dir}/work/cov_ribbon -M) -div ${output_dir}/work/cov_ribbon_norm_s$NeighborhoodSmoothing ${output_dir}/work/cov_norm_modulate
fslmaths ${output_dir}/work/cov_norm_modulate -mas ${output_dir}/work/ribbon_only.nii.gz ${output_dir}/work/cov_norm_modulate_ribbon

STD=$(fslstats ${output_dir}/work/cov_norm_modulate_ribbon -S)
echo $STD
MEAN=$(fslstats ${output_dir}/work/cov_norm_modulate_ribbon -M)
echo $MEAN
Lower=$(echo "$MEAN - ($STD * $Factor)" | bc -l)
echo $Lower
Upper=$(echo "$MEAN + ($STD * $Factor)" | bc -l)
echo $Upper

fslmaths ${output_dir}/work/mean -bin ${output_dir}/work/mask
fslmaths ${output_dir}/work/cov_norm_modulate -thr $Upper -bin -sub ${output_dir}/work/mask -mul -1 ${output_dir}/work/goodvoxels

# --------------------------------- PROJECT FUNCTIONAL DATA TO SURF --------------------------------------------------------
for hemi in left right; do

  if [ $hemi = "left" ]; then
    hemi_upper="L"
  elif [ $hemi = "right" ]; then
    hemi_upper="R"
  fi

  # mapping MEAN and COV
  for map in mean cov; do

    # -volume-to-surface-mapping WITH -volume-roi
    # NAME: ${map}_hemi-${hemi}_mesh-native.func.gii
    $path_wbcommand -volume-to-surface-mapping \
      ${output_dir}/work/${map}.nii.gz \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
      ${output_dir}/work/${map}_hemi-${hemi}_mesh-native.func.gii \
      -ribbon-constrained \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_wm.surf.gii \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_pial.surf.gii \
      -volume-roi \
      ${output_dir}/work/goodvoxels.nii.gz

    $path_wbcommand -metric-dilate \
      ${output_dir}/work/${map}_hemi-${hemi}_mesh-native.func.gii \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
      10 \
      ${output_dir}/work/${map}_hemi-${hemi}_mesh-native.func.gii \
      -nearest

    $path_wbcommand -metric-mask \
      ${output_dir}/work/${map}_hemi-${hemi}_mesh-native.func.gii \
      ${sub_anat}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_desc-medialwall_mask.shape.gii \
      ${output_dir}/work/${map}_hemi-${hemi}_mesh-native.func.gii

    # -volume-to-surface-mapping WITHOUT -volume-roi
    # NAME: ${map}_all_hemi-${hemi}_mesh-native.func.gii
    $path_wbcommand -volume-to-surface-mapping \
      ${output_dir}/work/${map}.nii.gz \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
      ${output_dir}/work/${map}_all_hemi-${hemi}_mesh-native.func.gii \
      -ribbon-constrained \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_wm.surf.gii \
      ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_pial.surf.gii

    $path_wbcommand -metric-mask \
      ${output_dir}/work/${map}_all_hemi-${hemi}_mesh-native.func.gii \
      ${sub_anat}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_desc-medialwall_mask.shape.gii \
      ${output_dir}/work/${map}_all_hemi-${hemi}_mesh-native.func.gii
  done

  # mapping GOODVOXELS
  # NAME: goodvoxels_hemi-${hemi}_mesh-native.shape.gii
  $path_wbcommand -volume-to-surface-mapping \
    ${output_dir}/work/goodvoxels.nii.gz \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
    ${output_dir}/work/goodvoxels_hemi-${hemi}_mesh-native.shape.gii \
    -ribbon-constrained \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_wm.surf.gii \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_pial.surf.gii

  $path_wbcommand -metric-mask \
    ${output_dir}/work/goodvoxels_hemi-${hemi}_mesh-native.shape.gii \
    ${sub_anat}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_desc-medialwall_mask.shape.gii \
    ${output_dir}/work/goodvoxels_hemi-${hemi}_mesh-native.shape.gii

  #  ribbon constrained mapping of fMRI volume to native anatomical surface (in func space)
  # NAME: func_hemi-${hemi}_mesh-native.func.gii
  $path_wbcommand -volume-to-surface-mapping \
    ${func} \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    -ribbon-constrained \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_wm.surf.gii \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_pial.surf.gii \
    -volume-roi \
    ${output_dir}/work/goodvoxels.nii.gz

  $path_wbcommand -metric-dilate \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    ${output_dir}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_mesh-native_space-func_midthickness.surf.gii \
    10 \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    -nearest

  $path_wbcommand -metric-mask \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    ${sub_anat}/anat/sub-${subid}_ses-${sesid}_hemi-${hemi_upper}_desc-medialwall_mask.shape.gii \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}"

done

rm -rf "${output_dir}/work"
