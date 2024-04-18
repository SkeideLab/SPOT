import numpy as np
from nilearn import plotting, surface
import nibabel as nib

# Load your fMRI data (example data, replace with your own data)
# Here, assume data is a 4D numpy array with dimensions (time, x, y, z)
VISPARC_PATH = (
    "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_"
    "desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
)
FUNC_PATH = (
    "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/func/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native{simulated}_bold.func.gii"
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

        # load retinotopy data
        visparc = nib.load(VISPARC_PATH.format(**ids))
        indices_v2 = get_indices_roi((3,4), visparc)

        for datasource in ["real", "simulated"]:
            print(f"Modeling {datasource} data on hemisphere {hemi}.")

            # get bold data for V1 and V2
            simulated = "_desc-simulated" if datasource == "simulated" else ""
            func = surface.load_surf_data(FUNC_PATH.format(**ids, simulated=simulated))
            func_v2 = func[indices_v2, :].astype(np.float64)
            surface_data=np.corrcoef(func_v2)                                                                                                                        
            #surface_data=surface_data[surface_data != 0]
            #print(surface_data.shape)
            #surface_data[surface_data < 0.1] = 0
            #surface_data=surface_data[surface_data != 0]
            median_correlation = np.median(surface_data)
            print(surface_data.shape)
            print("Median correlation coefficient:", median_correlation)




