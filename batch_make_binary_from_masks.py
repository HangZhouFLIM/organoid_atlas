"""
Utility for converting grayscale mask images into binary masks and adding them
as layers in an active napari viewer session.

Run within a napari console with an existing ``viewer`` variable:

>>> from batch_make_binary_from_masks import batch_make_binary_from_masks
>>> batch_make_binary_from_masks()
"""
import os

import imageio.v2 as imageio
import numpy as np
from qtpy.QtWidgets import QFileDialog


def batch_make_binary_from_masks():
    """Load mask images, binarize them, and save the outputs.

    The function prompts the user to pick one or more mask files. Each selected
    file is loaded, converted to a binary mask (non-zero values become 255),
    and added to the current napari viewer. A copy of each binary mask is saved
    alongside the original image with the suffix ``_binaryMask``.
    """
    # Use napari main window as parent for the file dialog
    parent = viewer.window._qt_window  # type: ignore[name-defined]

    # Let user select multiple mask images
    fnames, _ = QFileDialog.getOpenFileNames(
        parent,
        "Select mask images",
        "",
        "Image files (*.tif *.tiff *.png *.jpg *.jpeg)",
    )

    if not fnames:
        print("No files selected.")
        return

    for fname in fnames:
        print(f"Processing: {fname}")
        img = imageio.imread(fname)

        # Create binary mask: all non-zero pixels -> 255 (white)
        binary = (img > 0).astype(np.uint8) * 255

        # Add to napari as a new image layer
        layer_name = os.path.basename(fname) + "_binary"
        viewer.add_image(binary, name=layer_name, colormap="gray")  # type: ignore[name-defined]

        # Build output path in same folder with suffix "_binaryMask"
        folder, base = os.path.split(fname)
        root, ext = os.path.splitext(base)
        out_path = os.path.join(folder, root + "_binaryMask" + ext)

        # Save binary mask
        imageio.imwrite(out_path, binary)
        print(f"Saved binary mask: {out_path}")

    print("All done!")


if __name__ == "__main__":
    batch_make_binary_from_masks()
