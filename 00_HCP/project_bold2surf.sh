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
path_wbcommand=$2
anat_dir=$3
func_dir=$4
output_dir=$5
WORKBENCHDIR=$(dirname $path_wbcommand)
export WORKBENCHDIR

# inputs
sub_anat="${anat_dir}"
sub_func="${func_dir}"
func="${sub_func}/rfMRI_REST_hp0_clean.nii.gz"

#-------------------- outputs -----------------------------------------------#

output_dir_sub="${output_dir}/${subid}"
work_dir="${output_dir_sub}/work"
mkdir -p "${output_dir_sub}/anat" "$work_dir" "${output_dir_sub}/func"
ref_volume="${sub_func}/rfMRI_REST_hp0_clean_mean.nii.gz"

#-----------------parameters------------------------------#
NeighborhoodSmoothing="5"
Factor="0.5"
LeftGreyRibbonValue="1"
RightGreyRibbonValue="1"

###############

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
  func_on_surf="${output_dir_sub}/func/${subid}_hemi-${hemi_upper}_mesh-native_bold.func.gii"
  if [[ -f "${func_on_surf/\{hemi_upper\}/L}" && -f "${func_on_surf/\{hemi_upper\}/R}" ]]; then
  echo "Functional data on surface already exist at ${func_on_surf}. Skipping..."
  exit 0
  fi
  # ################
  # ribbon constrained mapping
  # adapted from fMRISurface/scripts/RibbonVolumeToSurfacemapping.sh
  # ################

  $path_wbcommand -create-signed-distance-volume \
    "${sub_anat}/${subid}.${hemi_upper}.white.native.surf.gii" \
    "${ref_volume}" \
    "$work_dir/${hemi}.white.native.nii.gz"

  $path_wbcommand -create-signed-distance-volume \
    "${sub_anat}/${subid}.${hemi_upper}.pial.native.surf.gii" \
    "${ref_volume}" \
    "${work_dir}/${hemi}.pial.native.nii.gz"

  # create ribbon mask between white and pial mask and fill it with GreyRibbonValue
  fslmaths "${work_dir}/${hemi}.white.native.nii.gz" -thr 0 -bin -mul 255 "${work_dir}/${hemi}.white_thr0.native.nii.gz"
  fslmaths "${work_dir}/${hemi}.white_thr0.native.nii.gz" -bin "${work_dir}/${hemi}.white_thr0.native.nii.gz"
  fslmaths "${work_dir}/${hemi}.pial.native.nii.gz" -uthr 0 -abs -bin -mul 255 "${work_dir}/${hemi}.pial_uthr0.native.nii.gz"
  fslmaths "${work_dir}/${hemi}.pial_uthr0.native.nii.gz" -bin "${work_dir}/${hemi}.pial_uthr0.native.nii.gz"
  fslmaths "${work_dir}/${hemi}.pial_uthr0.native.nii.gz" -mas "${work_dir}/${hemi}.white_thr0.native.nii.gz" -mul 255 "${work_dir}/${hemi}.ribbon.nii.gz"
  fslmaths "${work_dir}/${hemi}.ribbon.nii.gz" -bin -mul "$GreyRibbonValue" "${work_dir}/${hemi}.ribbon.nii.gz"
  rm "${work_dir}/${hemi}.white.native.nii.gz" "${work_dir}/${hemi}.white_thr0.native.nii.gz" "${work_dir}/${hemi}.pial.native.nii.gz" "${work_dir}/${hemi}.pial_uthr0.native.nii.gz"
done

# Combine ribbon masks for both hemispheres
fslmaths ${work_dir}/left.ribbon.nii.gz -add ${work_dir}/right.ribbon.nii.gz ${work_dir}/ribbon_only.nii.gz
rm ${work_dir}/left.ribbon.nii.gz ${work_dir}/right.ribbon.nii.gz

# compute mean and standard deviation of functional image
fslmaths ${func} -Tmean ${work_dir}/mean -odt float
fslmaths ${func} -Tstd ${work_dir}/std -odt float
fslmaths ${work_dir}/std -div ${work_dir}/mean ${work_dir}/cov
fslmaths ${work_dir}/cov -mas ${work_dir}/ribbon_only.nii.gz ${work_dir}/cov_ribbon

#TODO what and why is being smoothed here?
fslmaths ${work_dir}/cov_ribbon -div $(fslstats ${work_dir}/cov_ribbon -M) ${work_dir}/cov_ribbon_norm
fslmaths ${work_dir}/cov_ribbon_norm -bin -s $NeighborhoodSmoothing ${work_dir}/SmoothNorm
fslmaths ${work_dir}/cov_ribbon_norm -s $NeighborhoodSmoothing -div ${work_dir}/SmoothNorm -dilD ${work_dir}/cov_ribbon_norm_s$NeighborhoodSmoothing
fslmaths ${work_dir}/cov -div $(fslstats ${work_dir}/cov_ribbon -M) -div ${work_dir}/cov_ribbon_norm_s$NeighborhoodSmoothing ${work_dir}/cov_norm_modulate
fslmaths ${work_dir}/cov_norm_modulate -mas ${work_dir}/ribbon_only.nii.gz ${work_dir}/cov_norm_modulate_ribbon

STD=$(fslstats ${work_dir}/cov_norm_modulate_ribbon -S)
echo $STD
MEAN=$(fslstats ${work_dir}/cov_norm_modulate_ribbon -M)
echo $MEAN
Lower=$(echo "$MEAN - ($STD * $Factor)" | bc -l)
echo $Lower
Upper=$(echo "$MEAN + ($STD * $Factor)" | bc -l)
echo $Upper

fslmaths ${work_dir}/mean -bin ${work_dir}/mask
fslmaths ${work_dir}/cov_norm_modulate -thr $Upper -bin -sub ${work_dir}/mask -mul -1 ${work_dir}/goodvoxels

# --------------------------------- PROJECT FUNCTIONAL DATA TO SURF --------------------------------------------------------
for hemi in left right; do

  if [ $hemi = "left" ]; then
    hemi_upper="L"
  elif [ $hemi = "right" ]; then
    hemi_upper="R"
  fi
  func_on_surf="${output_dir_sub}/func/${subid}_hemi-${hemi_upper}_mesh-native_bold.func.gii"
  # mapping MEAN and COV
  for map in mean cov; do

    # -volume-to-surface-mapping WITH -volume-roi
    # NAME: ${map}_hemi-${hemi}_mesh-native.func.gii
    $path_wbcommand -volume-to-surface-mapping \
      ${work_dir}/${map}.nii.gz \
      ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
      ${work_dir}/${map}_hemi-${hemi}_mesh-native.func.gii \
      -ribbon-constrained \
      ${sub_anat}/${subid}.${hemi_upper}.white.native.surf.gii \
      ${sub_anat}/${subid}.${hemi_upper}.pial.native.surf.gii \
      -volume-roi \
      ${work_dir}/goodvoxels.nii.gz

    $path_wbcommand -metric-dilate \
      ${work_dir}/${map}_hemi-${hemi}_mesh-native.func.gii \
      ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
      10 \
      ${work_dir}/${map}_hemi-${hemi}_mesh-native.func.gii \
      -nearest

    $path_wbcommand -metric-mask \
      ${work_dir}/${map}_hemi-${hemi}_mesh-native.func.gii \
      ${sub_anat}/${subid}.${hemi_upper}.roi.native.shape.gii \
      ${work_dir}/${map}_hemi-${hemi}_mesh-native.func.gii

    # -volume-to-surface-mapping WITHOUT -volume-roi
    # NAME: ${map}_all_hemi-${hemi}_mesh-native.func.gii
    $path_wbcommand -volume-to-surface-mapping \
      ${work_dir}/${map}.nii.gz \
      ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
      ${work_dir}/${map}_all_hemi-${hemi}_mesh-native.func.gii \
      -ribbon-constrained \
      ${sub_anat}/${subid}.${hemi_upper}.white.native.surf.gii \
      ${sub_anat}/${subid}.${hemi_upper}.pial.native.surf.gii

    $path_wbcommand -metric-mask \
      ${work_dir}/${map}_all_hemi-${hemi}_mesh-native.func.gii \
      ${sub_anat}/${subid}.${hemi_upper}.roi.native.shape.gii \
      ${work_dir}/${map}_all_hemi-${hemi}_mesh-native.func.gii
  done

  # mapping GOODVOXELS
  # NAME: goodvoxels_hemi-${hemi}_mesh-native.shape.gii
  $path_wbcommand -volume-to-surface-mapping \
    ${work_dir}/goodvoxels.nii.gz \
    ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
    ${work_dir}/goodvoxels_hemi-${hemi}_mesh-native.shape.gii \
    -ribbon-constrained \
    ${sub_anat}/${subid}.${hemi_upper}.white.native.surf.gii \
    ${sub_anat}/${subid}.${hemi_upper}.pial.native.surf.gii

  $path_wbcommand -metric-mask \
    ${work_dir}/goodvoxels_hemi-${hemi}_mesh-native.shape.gii \
    ${sub_anat}/${subid}.${hemi_upper}.roi.native.shape.gii \
    ${work_dir}/goodvoxels_hemi-${hemi}_mesh-native.shape.gii

  #  ribbon constrained mapping of fMRI volume to native anatomical surface (in func space)
  # NAME: func_hemi-${hemi}_mesh-native.func.gii
  $path_wbcommand -volume-to-surface-mapping \
    ${func} \
    ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    -ribbon-constrained \
    ${sub_anat}/${subid}.${hemi_upper}.white.native.surf.gii \
    ${sub_anat}/${subid}.${hemi_upper}.pial.native.surf.gii \
    -volume-roi \
    ${work_dir}/goodvoxels.nii.gz

  $path_wbcommand -metric-dilate \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    ${sub_anat}/${subid}.${hemi_upper}.midthickness.native.surf.gii \
    10 \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    -nearest

  $path_wbcommand -metric-mask \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}" \
    ${sub_anat}/${subid}.${hemi_upper}.roi.native.shape.gii \
    "${func_on_surf/\{hemi_upper\}/$hemi_upper}"

done

rm -rf "${work_dir}"
