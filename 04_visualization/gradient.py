import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
import pyvista as pv
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

"""
Color Scheme:
Green: All components are positive.
Blue: R and A components are positive.
Cyan: A and S components are positive.
Magenta: R and S components are positive.
Red: Only the R component is positive.
Orange: Only the A component is positive.
Yellow: Only the S component is positive.
Purple: No components are positive.
"""

os.environ['MALLOC_CHECK_'] = '2'

VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")


def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        visparc.agg_data() == labels_area
    )[0]
    return indices_area


def extract_roi_mesh(brain_mesh, indices):
    # Extract the vertices and faces corresponding to the ROI
    mask = np.zeros(brain_mesh.n_points, dtype=bool)
    mask[indices] = True

    # Extract the ROI mesh
    roi_mesh = brain_mesh.extract_points(mask, include_cells=True)
    return roi_mesh


def determine_color(vector):
    if vector[0] > 0 and vector[1] > 0 and vector[2] > 0:
        return 'green'  # All components positive
    elif vector[0] > 0 and vector[1] > 0:
        return 'blue'   # X and Y components positive
    elif vector[1] > 0 and vector[2] > 0:
        return 'cyan'   # Y and Z components positive
    elif vector[0] > 0 and vector[2] > 0:
        return 'magenta'  # X and Z components positive
    elif vector[0] > 0:
        return 'red'    # Only X component positive
    elif vector[1] > 0:
        return 'orange'  # Only Y component positive
    elif vector[2] > 0:
        return 'yellow'  # Only Z component positive
    else:
        return 'purple'  # No component positive


# Create a plotter
# plotter = pv.Plotter(shape=(6, 4),border=False, off_screen=True)
df = pd.DataFrame()
for param in ["eccentricity"]:
    results_list = []
    print(param)

    for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
        if group == "neonates<37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_less_37.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 2
        elif group == "neonates>37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_less_37.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 3
        elif group == "fetal<29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/"
                "Averaged_younger_fetal_{hemi}_label-{param}_desc-real_roi-v2th00_metric.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 0
        elif group == "fetal>29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/"
                "Averaged_older_fetal_{hemi}_label-{param}_desc-real_roi-v2th00_metric.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 1
        elif group == "12-16y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_young.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 4
        elif group == "18-21y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_old.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 5

        for hemi in ["L", "R"]:
            for region in ["v2v", "v2d"]:
                if region == "v2v":
                    LABELS_V2 = 3
                elif region == "v2d":
                    LABELS_V2 = 4
                group_name = f"{group}_{hemi}"
                parameters = []
                visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
                indices_v2 = get_indices_roi(LABELS_V2, visparc)
            # FORMAT PATHS FOR INPUT AND OUTPUT
                for index, row in subject_info.iterrows():
                    if group == "12-16y" or group == "18-21y":
                        sub_id = subject_info.at[index, "sub_id"]
                        prefix_model = PREFIX_MODEL.format(
                            sub=sub_id,
                            hemi=hemi,
                            param=param,
                        )
                        input_path = prefix_model

                        ccf = surface.load_surf_data(input_path)
                        ccf_v0 = ccf[:, indices_v2].astype(np.float64)
                        parameters.append(ccf_v0)
                    else:
                        sub_id = subject_info.at[index, "sub_id"]
                        sess_id = subject_info.at[index, "sess_id"]
                        sub = sub_id.replace('sub-', '')
                        ses = sess_id.replace('ses-', '')
                        prefix_model = PREFIX_MODEL.format(
                            sub=sub,
                            ses=ses,
                            hemi=hemi,
                            param=param,
                        )

                        input_path = prefix_model

                        ccf = surface.load_surf_data(input_path)

                        ccf_v0 = ccf[:, indices_v2].astype(np.float64)
                        parameters.append(ccf_v0)

                averaged_map_path = AVER_MAP.format(
                    hemi=hemi,
                    param=param,
                )
                averaged_map = surface.load_surf_data(averaged_map_path)

                # Load the fsaverage mesh (GIFTI format)
                # Ensure you have the correct path to your fsaverage GIFTI file

                if hemi == "L":
                    mesh_path = '/data/p_02915/templates/template_fsaverage/fsaverage/surf/lh.inflated'
                    mesh_gii = nib.freesurfer.io.read_geometry(mesh_path)
                    if region == "v2v":
                        col = 1
                    elif region == "v2d":
                        col = 0
                elif hemi == "R":
                    mesh_path = '/data/p_02915/templates/template_fsaverage/fsaverage/surf/rh.inflated'
                    mesh_gii = nib.freesurfer.io.read_geometry(mesh_path)
                    if region == "v2v":
                        col = 2
                    elif region == "v2d":
                        col = 3

                if param == "eccentricity":
                    clim_range = (0, 20)
                elif param == "polarangle":
                    clim_range = (0, 180)

                # Extract vertices and faces
                vertices = mesh_gii[0][indices_v2]
                # Extract the ROI mesh
                roi_mesh = pv.PolyData(vertices)
                # Optionally, you can calculate faces using a Delaunay triangulation for irregular point clouds
                mesh_delaunay = roi_mesh.delaunay_2d()

                # Verify the number of points in roi_mesh
                print(f"Number of points in roi_mesh: {roi_mesh.n_points}")
                # Generate or load scalar data (if needed)
                num_points = roi_mesh.n_points
                stacked_arrays = np.stack(parameters)
                # Compute the mean along the new axis
                average_array = np.mean(stacked_arrays, axis=0)
                # Compute a gradient at each vertex (example: gradient of x-coordinate)
                gradients = average_array.T  # Replace with actual gradient calculation
                # Compute the gradient magnitude (norm)
                gradient_magnitude = np.linalg.norm(gradients, axis=1)
                # Normalize the gradient vectors for arrow direction
                normalized_gradients = gradients / \
                    gradient_magnitude[:, np.newaxis]
                # Create vectors for arrows using gradient directions and magnitudes
                vectors = normalized_gradients * \
                    gradient_magnitude[:, np.newaxis]
                mesh_delaunay.point_data['Scalars'] = averaged_map[indices_v2]
                # Interpolate scalar data to fill the faces of the mesh
                cell_data_mesh = mesh_delaunay.point_data_to_cell_data()
                plotter = pv.Plotter(border=False, off_screen=True)
                # plotter.subplot(row_p, col)

                # Add the filled mesh with interpolated scalars
                plotter.add_mesh(cell_data_mesh, scalars='Scalars', cmap='jet',
                                 clim=clim_range, opacity=0.9, show_scalar_bar=False)

                # Add the gradient vectors as arrows
                # plotter.add_arrows(mesh_delaunay.points, vectors, mag=6, color='red')  # Adjust 'mag' for scaling arrows
                # Add arrows to the plotter with differentiated colors
                for point, vector in zip(mesh_delaunay.points, vectors):
                    color = determine_color(vector)
                    plotter.add_arrows(point, vector, mag=6, color=color)

                # Set background color (optional)
                plotter.set_background('white')

                # Add a title to the plot
                # plot_title = f"{group}_{hemi}_{param}_v2v"
                # plotter.add_text(plot_title, position='upper_edge', font_size=12, color='black')

                # Set the camera position to face the direction of an example arrow
                # For this example, we assume the arrow is pointing from the origin to the point (1, 1, 1)
                arrow_direction = vectors[0]  # Example arrow direction vector
                # Example arrow start point
                arrow_start = np.array(mesh_delaunay.center)

                if hemi == "L" and region == "v2v":
                    arrow_end = arrow_start + (50, -50, -50)
                    plotter.camera_position = [
                        arrow_end, arrow_start, [0, 0, 1]]
                elif hemi == "L" and region == "v2d":
                    arrow_end = arrow_start + (50, -80, 50)
                    plotter.camera_position = [
                        arrow_end, arrow_start, [0, 0, 1]]
                    # plotter.camera.roll=-45
                elif hemi == "R" and region == "v2v":
                    arrow_end = arrow_start + (-50, -50, -50)
                    plotter.camera_position = [
                        arrow_end, arrow_start, [0, 0, 1]]
                elif hemi == "R" and region == "v2d":
                    arrow_end = arrow_start + (-50, -80, 50)
                    plotter.camera_position = [
                        arrow_end, arrow_start, [0, 0, 1]]
                    # plotter.camera.roll=-45

                plotter.camera.zoom(1.5)

                axes = pv.Axes()
                axes_actor = axes.axes_actor
                axes.axes_actor.shaft_type = 0
                axes_actor.x_axis_shaft_properties.color = (1, 0, 0)
                axes_actor.y_axis_shaft_properties.color = (0, 1, 0)
                axes_actor.z_axis_shaft_properties.color = (0, 0, 1)
                axes_actor.x_axis_label = 'R'
                axes_actor.y_axis_label = 'A'
                axes_actor.z_axis_label = 'S'
                plotter.add_orientation_widget(
                    axes_actor,
                    viewport=(0.01, 0.01, 0.2, 0.1),
                )

                # Set the camera position to face the arrow
                # plotter.camera_position = [arrow_end, arrow_start, [0, 0, 1]]

                # if row_p ==0:
                #    plotter.add_text(f"{hemi} {region}", position='upper_edge', font_size=12, color='black')

                # Show the plotter
                screenshot_filename = f"/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/{group}_{param}_{hemi}_{region}_gradient_org.pdf"
                plotter.save_graphic(screenshot_filename)


# Load the screenshot image into Matplotlib
# img = mpimg.imread(screenshot_filename)

# Create a Matplotlib figure and axis
# fig, ax = plt.subplots(figsize=(6, 8))  # Adjust size as needed

# Display the screenshot image without empty spaces
# ax.imshow(img, extent=[0, img.shape[1], 0, img.shape[0]])
# Define text and positions
# text_left = ['a', 'b','c', 'd', 'e', 'f']
# text_top_row1 = ['Left', 'Right']
# text_top_row2 = ['V2d', 'V2v', 'V2v', 'V2d']


# Set font properties
# fontsize = 12  # Adjust as needed

# Add text annotations
# for i, text in enumerate(text_left):
#    fig.text(0.06, 1 - (i * 0.09) - 0.25, text, ha='center', va='center', fontsize=12, weight='bold')

# for i, text in enumerate(text_top_row1):
#    fig.text((i * 0.48) + 0.255, 0.82, text, ha='center', va='center', fontsize=10)

# for i, text in enumerate(text_top_row2):
#     fig.text((i * 0.23) + 0.16, 0.8, text, ha='center', va='center', fontsize=10)


# Remove axis ticks and labels
# ax.axis('off')
# output_filename = f"/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/{param}_gradient_title.png"
# Make the figure tight
# plt.tight_layout()
# Save or show the annotated plot
# plt.savefig(output_filename, bbox_inches='tight')
# plt.show()
