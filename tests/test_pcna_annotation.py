import pandas as pd
import numpy as np
import tifffile

from pcna_annotation import (
    append_annotation, export_curated, padded_crop, prepare_manifest, stable_crop_id,
    split_for_experiment,
)


def test_padded_crop_preserves_shape_dtype_and_edge_pixels():
    image = np.arange(16, dtype=np.uint16).reshape(4, 4)
    crop = padded_crop(image, 0, 0, 4)
    assert crop.shape == (4, 4)
    assert crop.dtype == np.uint16
    assert crop[2, 2] == image[0, 0]
    assert crop[0, 0] == 0


def test_split_is_deterministic_and_group_level():
    assert split_for_experiment("2026-08-06_run_a") == split_for_experiment("2026-08-06_run_a")
    assert split_for_experiment("2026-08-06_run_a") in {"train", "validation", "test"}


def test_prepare_reconciles_stack_mask_and_writes_uint16_crop(tmp_path):
    images, masks, output = tmp_path / "images", tmp_path / "masks", tmp_path / "out"
    images.mkdir(); masks.mkdir()
    stack = np.arange(3 * 8 * 8, dtype=np.uint16).reshape(3, 8, 8)
    tifffile.imwrite(images / "sample_XY001_C1.tif", stack)
    mask = np.zeros((8, 8), dtype=np.uint16); mask[3:5, 3:5] = 7
    tifffile.imwrite(masks / "sample_XY001_T002_C1.tiff", mask)
    nuclei = pd.DataFrame([{"position": "XY1", "frame": 2, "label": 7, "centroid-0": 3.5, "centroid-1": 3.5, "area": 4}])
    nuclei_path = tmp_path / "nuclei.csv"; nuclei.to_csv(nuclei_path, index=False)
    manifest = prepare_manifest("day_a", nuclei_path, images, masks, output, crop_size=4)
    assert manifest.iloc[0].slice_index == 1
    assert tifffile.imread(manifest.iloc[0].crop_path).dtype == np.uint16
    assert manifest.iloc[0].crop_id == stable_crop_id("day_a", "XY001", 2, 7)


def test_export_excludes_disagreement_and_keeps_agreed_label(tmp_path):
    manifest = pd.DataFrame([{"crop_id": "a"}, {"crop_id": "b"}])
    manifest_path = tmp_path / "manifest.csv"; manifest.to_csv(manifest_path, index=False)
    log = tmp_path / "annotations.csv"
    append_annotation(log, {"crop_id": "a", "annotator": "one", "timestamp_utc": "2026-01-01T00:00:00Z", "label": "interphase"})
    append_annotation(log, {"crop_id": "a", "annotator": "two", "timestamp_utc": "2026-01-01T00:01:00Z", "label": "active_mitosis"})
    append_annotation(log, {"crop_id": "b", "annotator": "one", "timestamp_utc": "2026-01-01T00:00:00Z", "label": "new_daughter"})
    curated, disagreements = export_curated(manifest_path, log, tmp_path / "curated.csv")
    assert curated.crop_id.tolist() == ["b"]
    assert disagreements.crop_id.tolist() == ["a"]
