import re
from pathlib import Path
from nilearn._utils.helpers import stringify_path
from nilearn import _utils
import numpy as np
import warnings
from nilearn._utils.niimg import _get_data
from nilearn.image.resampling import (
    get_bounds,
    BoundingBoxError,
    to_matrix_vector,
    _resample_one_img,
)
from nilearn.image import crop_img
from scipy import linalg


def get_mesh(path, mesh_name, hemis=("left", "right")):
    """Get all kinds of matching meshes from input path that
    describes one result file in derivatives folder


    Args:
        path (str): path to results file
        mesh_name (str): 'pial', 'white', 'sulc', 'inflated'


    Returns:
        dict: dictionary with keys 'left','right',
        containing path to found mesh file as string or pathlib.Path
    """
    path = str(path)
    sub = get_sub(path)
    ses = get_ses(path)
    meshes = {}
    for hemi in hemis:
        mesh = Path(path.split("dhcp", maxsplit=1)[0]) / "dhcp_anat_pipeline"
        mesh = list(
            mesh.glob(
                f"sub-{sub}/ses-{ses}/anat/*_hemi-{hemi.capitalize()[0]}*_{mesh_name}.surf.gii"
            )
        )
        assert (
            len(mesh) < 2
        ), f"[ERROR] More than 1 available surface for mesh {mesh_name} and path {path}."
        assert (
            len(mesh) > 0
        ), f"[ERROR] No available surface found for mesh {mesh_name} and path {path}."

        meshes[hemi] = str(mesh[0])

    return meshes


def get_brainmask(path):
    return str(list(path.parent.glob("*brainmask*"))[0])


def get_ses(path):
    return re.search(r"(?<=ses-)[0-9A-Za-z]*", str(path)).group(0)


def get_sub(path):
    return re.search(r"(?<=sub-)[0-9A-Za-z]*", str(path)).group(0)


def resample_img(
    img,
    target_affine=None,
    target_shape=None,
    interpolation="continuous",
    copy=True,
    order="F",
    clip=True,
    fill_value=0,
    force_resample=False,
    transform_affine=None,
):
    """Resample a Niimg-like object.
    Adapted from nilearn.image.resampling

    Parameters
    ----------
    img : Niimg-like object
        See :ref:`extracting_data`.
        Image(s) to resample.

    target_affine : numpy.ndarray, optional
        If specified, the image is resampled corresponding to this new affine.
        target_affine can be a 3x3 or a 4x4 matrix. (See notes)

    target_shape : tuple or list, optional
        If specified, the image will be resized to match this new shape.
        len(target_shape) must be equal to 3.
        If target_shape is specified, a target_affine of shape (4, 4)
        must also be given. (See notes)

    interpolation : str, optional
        Can be 'continuous', 'linear', or 'nearest'. Indicates the resample
        method. Default='continuous'.

    copy : bool, optional
        If True, guarantees that output array has no memory in common with
        input array.
        In all cases, input images are never modified by this function.
        Default=True.

    order : "F" or "C", optional
        Data ordering in output array. This function is slightly faster with
        Fortran ordering. Default='F'.

    clip : bool, optional
        If True (default) all resampled image values above max(img) and
        under min(img) are clipped to min(img) and max(img). Note that
        0 is added as an image value for clipping, and it is the padding
        value when extrapolating out of field of view.
        If False no clip is performed.
        Default=True.

    fill_value : float, optional
        Use a fill value for points outside of input volume. Default=0.

    force_resample : bool, optional
        Intended for testing, this prevents the use of a padding optimization.
        Default=False.

    Returns
    -------
    resampled : nibabel.Nifti1Image
        input image, resampled to have respectively target_shape and
        target_affine as shape and affine.

    See Also
    --------
    nilearn.image.resample_to_img

    Notes
    -----
    **BoundingBoxError**
    If a 4x4 transformation matrix (target_affine) is given and all of the
    transformed data points have a negative voxel index along one of the
    axis, then none of the data will be visible in the transformed image
    and a BoundingBoxError will be raised.

    If a 4x4 transformation matrix (target_affine) is given and no target
    shape is provided, the resulting image will have voxel coordinate
    (0, 0, 0) in the affine offset (4th column of target affine) and will
    extend far enough to contain all the visible data and a margin of one
    voxel.

    **3x3 transformation matrices**
    If a 3x3 transformation matrix is given as target_affine, it will be
    assumed to represent the three coordinate axes of the target space. In
    this case the affine offset (4th column of a 4x4 transformation matrix)
    as well as the target_shape will be inferred by resample_img, such that
    the resulting field of view is the closest possible (with a margin of
    1 voxel) bounding box around the transformed data.

    In certain cases one may want to obtain a transformed image with the
    closest bounding box around the data, which at the same time respects
    a voxel grid defined by a 4x4 affine transformation matrix. In this
    case, one resamples the image using this function given the target
    affine and no target shape. One then uses crop_img on the result.

    **NaNs and infinite values**
    This function handles gracefully NaNs and infinite values in the input
    data, however they make the execution of the function much slower.

    **Handling non-native endian in given Nifti images**
    This function automatically changes the byte-ordering information
    in the image dtype to new byte order. From non-native to native, which
    implies that if the given image has non-native endianness then the output
    data in Nifti image will have native dtype. This is only the case when
    if the given target_affine (transformation matrix) is diagonal and
    homogeneous.

    """
    from nilearn.image import new_img_like  # avoid circular imports

    # Do as many checks as possible before loading data, to avoid potentially
    # costly calls before raising an exception.
    if target_shape is not None and target_affine is None:
        raise ValueError(
            "If target_shape is specified, target_affine should" " be specified too."
        )

    if target_shape is not None and not len(target_shape) == 3:
        raise ValueError(
            "The shape specified should be the shape of "
            "the 3D grid, and thus of length 3. "
            f"{target_shape} was specified."
        )

    if target_shape is not None and target_affine.shape == (3, 3):
        raise ValueError(
            "Given target shape without anchor vector: "
            "Affine shape should be (4, 4) and not (3, 3)"
        )

    allowed_interpolations = ("continuous", "linear", "nearest")
    if interpolation not in allowed_interpolations:
        raise ValueError(
            f"interpolation must be one of {allowed_interpolations}.\n"
            f" Got '{interpolation}' instead."
        )

    if interpolation == "continuous":
        interpolation_order = 3
    elif interpolation == "linear":
        interpolation_order = 1
    elif interpolation == "nearest":
        interpolation_order = 0

    input_img_is_string = False
    img = stringify_path(img)
    if isinstance(img, str):
        # Avoid a useless copy
        input_img_is_string = True

    img = _utils.check_niimg(img)
    shape = img.shape
    affine = img.affine

    # If later on we want to impute sform using qform add this condition
    # see : https://github.com/nilearn/nilearn/issues/3168#issuecomment-1159447771 # noqa:E501
    if hasattr(img, "get_sform"):  # NIfTI images only
        _, sform_code = img.get_sform(coded=True)
        if not sform_code:
            warnings.warn(
                "The provided image has no sform in its header. "
                "Please check the provided file. "
                "Results may not be as expected."
            )

    # noop cases
    if target_affine is None and target_shape is None:
        if copy and not input_img_is_string:
            img = _utils.copy_img(img)
        return img
    if (
        np.shape(target_affine) == np.shape(affine)
        and np.allclose(target_affine, affine)
        and np.array_equal(target_shape, shape)
    ):
        return img
    if target_affine is not None:
        target_affine = np.asarray(target_affine)

    if np.all(np.array(target_shape) == shape[:3]) and np.allclose(
        target_affine, affine
    ):
        if copy and not input_img_is_string:
            img = _utils.copy_img(img)
        return img

    # We now know that some resampling must be done.
    # The value of "copy" is of no importance: output is always a separate
    # array.
    data = _get_data(img)

    # Get a bounding box for the transformed data
    # Embed target_affine in 4x4 shape if necessary
    if target_affine.shape == (3, 3):
        missing_offset = True
        target_affine_tmp = np.eye(4)
        target_affine_tmp[:3, :3] = target_affine
        target_affine = target_affine_tmp
    else:
        missing_offset = False
        target_affine = target_affine.copy()
    if transform_affine is None:
        transform_affine = np.linalg.inv(target_affine).dot(affine)
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = get_bounds(
        data.shape[:3], transform_affine
    )

    # if target_affine is (3, 3), then calculate
    # offset from bounding box and update bounding box
    # to be in the voxel coordinates of the calculated 4x4 affine
    if missing_offset:
        offset = target_affine[:3, :3].dot([xmin, ymin, zmin])
        target_affine[:3, 3] = offset
        (xmin, xmax), (ymin, ymax), (zmin, zmax) = (
            (0, xmax - xmin),
            (0, ymax - ymin),
            (0, zmax - zmin),
        )

    # if target_shape is not given (always the case with 3x3
    # transformation matrix and sometimes the case with 4x4
    # transformation matrix), then set it to contain the bounding
    # box by a margin of 1 voxel
    if target_shape is None:
        target_shape = (
            int(np.ceil(xmax)) + 1,
            int(np.ceil(ymax)) + 1,
            int(np.ceil(zmax)) + 1,
        )

    # Check whether transformed data is actually within the FOV
    # of the target affine
    if xmax < 0 or ymax < 0 or zmax < 0:
        raise BoundingBoxError(
            "The field of view given "
            "by the target affine does "
            "not contain any of the data"
        )

    if np.all(target_affine == affine):
        # Small trick to be more numerically stable
        transform_affine = np.eye(4)
    else:
        if transform_affine is None:
            transform_affine = np.dot(linalg.inv(affine), target_affine)
    A, b = to_matrix_vector(transform_affine)

    data_shape = list(data.shape)
    # Make sure that we have a list here
    if isinstance(target_shape, np.ndarray):
        target_shape = target_shape.tolist()
    target_shape = tuple(target_shape)

    if interpolation == "continuous" and data.dtype.kind == "i":
        # cast unsupported data types to closest support dtype
        aux = data.dtype.name.replace("int", "float")
        aux = aux.replace("ufloat", "float").replace("floatc", "float")
        if aux in ["float8", "float16"]:
            aux = "float32"
        warnings.warn(f"Casting data from {data.dtype.name} to {aux}")
        resampled_data_dtype = np.dtype(aux)
    else:
        resampled_data_dtype = data.dtype

    # Since the release of 0.17, resampling nifti images have some issues
    # when affine is passed as 1D array and if data is of non-native
    # endianness.
    # See issue https://github.com/nilearn/nilearn/issues/1445.
    # If affine is passed as 1D, scipy uses _nd_image.zoom_shift rather
    # than _geometric_transform (2D) where _geometric_transform is able
    # to swap byte order in scipy later than 0.15 for nonnative endianness.

    # We convert to 'native' order to not have any issues either with
    # 'little' or 'big' endian data dtypes (non-native endians).
    if len(A.shape) == 1 and not resampled_data_dtype.isnative:
        resampled_data_dtype = resampled_data_dtype.newbyteorder("N")

    # Code is generic enough to work for both 3D and 4D images
    other_shape = data_shape[3:]
    resampled_data = np.zeros(
        list(target_shape) + other_shape,
        order=order,
        dtype=resampled_data_dtype,
    )

    all_img = (slice(None),) * 3

    # if (A == I OR some combination of permutation(I) and sign-flipped(I)) AND
    # all(b == integers):
    if (
        np.all(np.eye(3) == A)
        and all(bt == np.round(bt) for bt in b)
        and not force_resample
    ):
        # TODO: also check for sign flips
        # TODO: also check for permutations of I

        # ... special case: can be solved with padding alone
        # crop source image and keep N voxels offset before/after volume
        cropped_img, offsets = crop_img(img, pad=False, return_offset=True)

        # TODO: flip axes that are flipped
        # TODO: un-shuffle permuted dimensions

        # offset the original un-cropped image indices by the relative
        # translation, b.
        indices = [
            (int(off.start - dim_b), int(off.stop - dim_b))
            for off, dim_b in zip(offsets[:3], b[:3])
        ]

        # If image are not fully overlapping, place only portion of image.
        slices = []
        for dimsize, index in zip(resampled_data.shape, indices):
            slices.append(slice(np.max((0, index[0])), np.min((dimsize, index[1]))))
        slices = tuple(slices)

        # ensure the source image being placed isn't larger than the dest
        subset_indices = tuple(slice(0, s.stop - s.start) for s in slices)
        resampled_data[slices] = _get_data(cropped_img)[subset_indices]
    else:
        # If A is diagonal, ndimage.affine_transform is clever enough to use a
        # better algorithm.
        if np.all(np.diag(np.diag(A)) == A):
            A = np.diag(A)
        # Iterate over a set of 3D volumes, as the interpolation problem is
        # separable in the extra dimensions. This reduces the
        # computational cost
        for ind in np.ndindex(*other_shape):
            _resample_one_img(
                data[all_img + ind],
                A,
                b,
                target_shape,
                interpolation_order,
                out=resampled_data[all_img + ind],
                copy=not input_img_is_string,
                fill_value=fill_value,
            )

    if clip:
        # force resampled data to have a range contained in the original data
        # preventing ringing artefact
        # We need to add zero as a value considered for clipping, as it
        # appears in padding images.
        vmin = min(np.nanmin(data), 0)
        vmax = max(np.nanmax(data), 0)
        resampled_data.clip(vmin, vmax, out=resampled_data)

    return new_img_like(img, resampled_data, target_affine)


def get_transform(p, from_to):
    assert from_to == "from-bold_to-T2w", "can only transform from bold to t2w atm"
    transform = list((p.parent.parent / "xfm").glob(f"*{from_to}*"))
    assert len(transform) < 2, "multiple transforms found!"
    assert len(transform) > 0, "no transform found!"
    return str(transform[0])
