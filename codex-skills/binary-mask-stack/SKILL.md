---
name: binary-mask-stack
description: Convert grayscale or labeled mask images into binary masks and ordered TIFF stacks. Use when the user asks to batch-process Cellpose or similar mask PNG/TIF files such as *_cp_masks.png, masks.png, or segmentation masks in project subfolders, save per-image binary masks, and stack them by filename number sequence into multi-page TIFF files in the original folders.
---

# Binary Mask Stack

## Overview

Use this skill to batch-convert mask images into binary masks and save one ordered TIFF stack per source folder. The bundled script applies the same rule used in the original napari helper: every nonzero pixel becomes 255, and every zero pixel remains 0.

## Quick Start

Prefer the bundled script:

```powershell
$PY = "<python.exe>"
$SCRIPT = Join-Path $env:USERPROFILE ".codex\skills\binary-mask-stack\scripts\binarize_stack_masks.py"
& $PY $SCRIPT "Z:\path\to\project"
```

For non-Cellpose filenames, pass a pattern:

```powershell
$PY = "<python.exe>"
$SCRIPT = Join-Path $env:USERPROFILE ".codex\skills\binary-mask-stack\scripts\binarize_stack_masks.py"
& $PY $SCRIPT "Z:\path\to\project" --pattern "*masks.png"
```

## Workflow

1. Locate the project root and identify mask files. Default to recursive search for `*cp_masks.png`; use `--pattern` when the project uses another naming convention.
2. Group files by their parent folder.
3. Sort files within each folder by the numeric tokens in the filename stem, so names like `sample_1_02`, `sample_2_04`, and `sample_10_20` stack in numeric order.
4. Convert each mask to an 8-bit binary image with values 0 and 255.
5. Save per-image binary outputs beside the source masks as `<source_stem>_binaryMask<ext>`, unless `--no-binary-images` is used.
6. Save one TIFF stack in each source folder as `<folder_name>_binaryMask_stack.tif`, unless `--stack-name` is supplied.
7. Verify that each stack has only values `[0, 255]`, expected dimensions, and the expected number of planes.

## Script Notes

- Use `--dry-run` first when processing a new naming scheme.
- Use `--flat` to process only files directly under the root instead of recursing through subfolders.
- Use `--overwrite` when regenerating existing outputs.
- The script skips files whose stem already ends in `_binaryMask` to avoid reprocessing generated masks.
- If source images in the same folder have different XY shapes, stop and report the inconsistent files instead of writing a broken stack.

## Dependencies

Run the script with a Python environment that has `numpy`, `imageio`, and `tifffile`.
