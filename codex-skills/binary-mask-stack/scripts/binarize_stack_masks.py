#!/usr/bin/env python
"""Batch-binarize mask images and save ordered TIFF stacks."""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import imageio.v2 as imageio
    import numpy as np
    import tifffile
except ModuleNotFoundError as exc:
    missing = exc.name or "a required package"
    raise SystemExit(
        f"Missing Python package: {missing}. Install numpy, imageio, and tifffile, "
        "or run this script with an environment that already has them."
    ) from exc


def sequence_key(path: Path) -> tuple[int, tuple[int, ...], str]:
    """Sort by all numeric tokens in the filename, falling back to lowercase name."""
    nums = tuple(int(match) for match in re.findall(r"\d+", path.stem))
    return (0, nums, path.stem.lower()) if nums else (1, (), path.stem.lower())


def find_mask_files(root: Path, pattern: str, recursive: bool) -> list[Path]:
    iterator = root.rglob(pattern) if recursive else root.glob(pattern)
    files = [
        path
        for path in iterator
        if path.is_file() and not path.stem.endswith("_binaryMask")
    ]
    return sorted(files, key=lambda p: (str(p.parent).lower(), sequence_key(p), p.name.lower()))


def make_binary(path: Path) -> "np.ndarray":
    img = imageio.imread(path)
    if img.ndim == 3 and img.shape[-1] in (3, 4):
        # Mask PNGs should normally be 2D. If a reader returns channels, collapse
        # by treating any nonzero channel as foreground.
        img = np.any(img > 0, axis=-1)
    return (img > 0).astype(np.uint8) * 255


def output_binary_path(path: Path, suffix: str) -> Path:
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


def stack_path_for(folder: Path, stack_name: str | None) -> Path:
    name = stack_name if stack_name else f"{folder.name}_binaryMask_stack.tif"
    if not name.lower().endswith((".tif", ".tiff")):
        name = f"{name}.tif"
    return folder / name


def process_folder(
    folder: Path,
    paths: list[Path],
    binary_suffix: str,
    stack_name: str | None,
    write_binary_images: bool,
    overwrite: bool,
    dry_run: bool,
) -> dict[str, object]:
    binaries: list[np.ndarray] = []
    ordered = sorted(paths, key=lambda p: (sequence_key(p), p.name.lower()))

    print(f"Folder: {folder}")
    for path in ordered:
        binary_path = output_binary_path(path, binary_suffix)
        print(f"  {path.name} -> {binary_path.name}")
        if dry_run:
            continue

        binary = make_binary(path)
        binaries.append(binary)

        if write_binary_images:
            if binary_path.exists() and not overwrite:
                print(f"    skip existing binary image: {binary_path.name}")
            else:
                imageio.imwrite(binary_path, binary)

    target_stack = stack_path_for(folder, stack_name)
    if dry_run:
        print(f"  stack -> {target_stack.name}")
        return {"folder": folder, "count": len(ordered), "stack": target_stack, "dry_run": True}

    shapes = defaultdict(list)
    for path, binary in zip(ordered, binaries):
        shapes[binary.shape].append(path.name)
    if len(shapes) != 1:
        details = "; ".join(f"{shape}: {names}" for shape, names in shapes.items())
        raise ValueError(f"Cannot stack {folder}: inconsistent image shapes. {details}")

    stack = np.stack(binaries, axis=0)
    if target_stack.exists() and not overwrite:
        print(f"  skip existing stack: {target_stack.name}")
    else:
        tifffile.imwrite(
            target_stack,
            stack,
            imagej=True,
            photometric="minisblack",
            metadata={"axes": "ZYX"},
        )
        print(f"  stack -> {target_stack.name} shape={stack.shape} dtype={stack.dtype}")

    return {
        "folder": folder,
        "count": len(ordered),
        "stack": target_stack,
        "shape": stack.shape,
        "dtype": str(stack.dtype),
        "values": np.unique(stack).tolist(),
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert mask images to binary masks and save ordered TIFF stacks per folder."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Project root folder to search. Defaults to the current directory.",
    )
    parser.add_argument(
        "--pattern",
        default="*cp_masks.png",
        help='Mask filename glob pattern. Default: "*cp_masks.png".',
    )
    parser.add_argument(
        "--flat",
        action="store_true",
        help="Search only the root folder instead of recursively searching subfolders.",
    )
    parser.add_argument(
        "--binary-suffix",
        default="_binaryMask",
        help='Suffix for per-image binary masks. Default: "_binaryMask".',
    )
    parser.add_argument(
        "--stack-name",
        default=None,
        help="Optional stack filename to use in every folder. Default: <folder>_binaryMask_stack.tif.",
    )
    parser.add_argument(
        "--no-binary-images",
        action="store_true",
        help="Only write TIFF stacks; do not save per-image binary masks.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing binary images and TIFF stacks.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the files and output paths without writing anything.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        raise SystemExit(f"Root folder does not exist: {root}")
    if not root.is_dir():
        raise SystemExit(f"Root path is not a folder: {root}")

    files = find_mask_files(root, args.pattern, recursive=not args.flat)
    if not files:
        raise SystemExit(f"No files matching {args.pattern!r} found under {root}")

    grouped: dict[Path, list[Path]] = defaultdict(list)
    for path in files:
        grouped[path.parent].append(path)

    print(f"Found {len(files)} mask file(s) in {len(grouped)} folder(s).")
    summaries = []
    for folder in sorted(grouped, key=lambda p: str(p).lower()):
        summaries.append(
            process_folder(
                folder=folder,
                paths=grouped[folder],
                binary_suffix=args.binary_suffix,
                stack_name=args.stack_name,
                write_binary_images=not args.no_binary_images,
                overwrite=args.overwrite,
                dry_run=args.dry_run,
            )
        )

    if not args.dry_run:
        print("Verification:")
        for item in summaries:
            print(
                f"  {item['stack']}: planes={item['count']} "
                f"shape={item['shape']} dtype={item['dtype']} values={item['values']}"
            )
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
