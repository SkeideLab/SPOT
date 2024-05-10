import numpy as np
from nilearn import plotting, surface
import nibabel as nib

# Load your fMRI data (example data, replace with your own data)
# Here, assume data is a 4D numpy array with dimensions (time, x, y, z)

FUNC_PATH = (
    "{root_dir}ccfmodel/sub-{sub}/ses-{ses}/" #_sigma20_noise0_r.gii
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native{simulated}_corr_noise0_rss.gii"
)

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data() == labels_area[1]
        )
    )[0]
    return indices_area

for hemi in ["L", "R"]:
        if hemi == "L":
            hemi_down = "left"
        elif hemi == "R":
            hemi_down = "right"
        ids = {
            "sub": "CC00056XX07",
            "ses": "10700",
            "hemi": hemi,
            "hemi_down": hemi_down,            
            "root_dir": "/data/p_02915/dhcp_derivatives_SPOT/"
            }


        for datasource in ["simulated"]:
            print(f"Modeling {datasource} data on hemisphere {hemi}.")

            # get bold data for V1 and V2
            simulated = "_desc-simulated" if datasource == "simulated" else ""
            surface_data = surface.load_surf_data(FUNC_PATH.format(**ids, simulated=simulated))
            #surface_data=surface_data[surface_data != 0]
            #print(surface_data.shape)
            #surface_data[surface_data < 0.1] = 0
            surface_data=surface_data[surface_data != 0]
            median_correlation = np.mean(surface_data)
            print(surface_data.shape)
            print("Median correlation coefficient:", median_correlation)




