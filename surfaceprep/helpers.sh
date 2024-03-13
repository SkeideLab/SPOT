#!/bin/bash
set -u -x -e

# TODO add this to flux/surfaceprep/make_fsavg_to_native_warp.sh
cd /data/p_02495/templates/template_fsaverage/fsaverage/surf
FREESURFER mris_expand -thickness lh.white 0.5 lh.midthickness
FREESURFER mris_expand -thickness rh.white 0.5 rh.midthickness

FREESURFER mris_convert \
    /data/p_02495/templates/template_fsaverage/fsaverage/surf/lh.midthickness \
    /data/p_02495/templates/template_fsaverage/fsaverage/surf/lh_midthickness_surf.gii

FREESURFER mris_convert \
    /data/p_02495/templates/template_fsaverage/fsaverage/surf/rh.midthickness \
    /data/p_02495/templates/template_fsaverage/fsaverage/surf/rh_midthickness_surf.gii
