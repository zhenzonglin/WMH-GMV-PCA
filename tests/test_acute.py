"""合成影像测试：负荷定义、空间契约、插值与文件完整性。"""
from pathlib import Path
import csv

import nibabel as nib
import numpy as np
import pytest

from wmh_gmv.acute import extract_acute


def save_image(data, affine, path, spatial_unit="mm"):
    image = nib.Nifti1Image(data, affine)
    image.header.set_xyzt_units(spatial_unit)
    nib.save(image, path)


def inputs(tmp_path: Path):
    shape = (4, 4, 4)
    affine = np.diag([2, 2, 2, 1])
    mask = np.zeros(shape, dtype=np.uint8)
    mask[0, 0, 0] = 1
    mask[3, 3, 3] = 1
    labels = np.zeros(shape, dtype=np.uint8)
    labels[:2] = 1
    labels[2:] = 2
    mask_path, atlas_path = tmp_path / "mask.nii.gz", tmp_path / "atlas.nii.gz"
    save_image(mask, affine, mask_path)
    save_image(labels, affine, atlas_path)
    label_file = tmp_path / "labels.tsv"
    label_file.write_text("index\tname\n1\tleft\n2\tright\n", encoding="utf-8")
    atlas = dict(name="test", kind="labels", path=str(atlas_path), labels=str(label_file),
                 template="synthetic-test", sha256=None, license="synthetic", source="generated",
                 include_labels=[], exclude_labels=[0])
    config = dict(project=dict(mode="synthetic"), paths=dict(atlases=[atlas]),
                  analysis=dict(acute=dict(template="synthetic-test", probability_scale=100, expected_spatial_unit="mm",
                                            allow_grid_resampling=True, binary_tolerance=1e-6, grid_affine_tolerance=1e-5)))
    subject = dict(subject_key="C3:0001", participant_id="0001", cohort="C3", lesion_mask=str(mask_path), lesion_template="synthetic-test")
    return subject, config


def read_output(path: Path):
    with path.open() as stream:
        return next(csv.DictReader(stream, delimiter="\t"))


def test_exact_fraction_and_volume(tmp_path):
    subject, config = inputs(tmp_path)
    result = extract_acute(subject, config, tmp_path / "out")
    row = read_output(tmp_path / "out" / "acute.tsv")
    assert result["status"] == "completed"
    assert row["participant_id"] == "0001"
    assert float(row["acute_volume_ml"]) == pytest.approx(.016)
    assert float(row["acute_test_1"]) == pytest.approx(1 / 32)
    assert float(row["acute_uncovered_fraction"]) == 0


def test_probability_weighted_not_regression(tmp_path):
    subject, config = inputs(tmp_path)
    atlas = config["paths"]["atlases"][0]
    prob = np.zeros((4, 4, 4, 2), dtype=np.float32)
    prob[..., 0] = 25
    prob[0, 0, 0, 0] = 100
    prob[..., 1] = 50
    save_image(prob, np.diag([2, 2, 2, 1]), atlas["path"])
    Path(atlas["labels"]).write_text("index\tname\n0\ttract0\n1\ttract1\n")
    atlas.update(kind="probability", exclude_labels=[])
    extract_acute(subject, config, tmp_path / "out")
    row = read_output(tmp_path / "out" / "acute.tsv")
    assert float(row["acute_test_0"]) == pytest.approx(125 / (63 * 25 + 100))
    assert float(row["acute_test_1"]) == pytest.approx(2 / 64)


def test_template_mismatch_is_not_fixed_by_resampling(tmp_path):
    subject, config = inputs(tmp_path)
    subject["lesion_template"] = "another-template"
    with pytest.raises(ValueError, match="template"):
        extract_acute(subject, config, tmp_path / "out")


def test_mask_must_be_binary(tmp_path):
    subject, config = inputs(tmp_path)
    bad = np.full((4, 4, 4), .5)
    save_image(bad, np.eye(4), subject["lesion_mask"])
    with pytest.raises(ValueError, match="binary"):
        extract_acute(subject, config, tmp_path / "out")


def test_missing_mask_waits(tmp_path):
    subject, config = inputs(tmp_path)
    subject["lesion_mask"] = None
    result = extract_acute(subject, config, tmp_path / "out")
    assert result["status"] == "waiting_input"


def test_atlas_hash_mismatch(tmp_path):
    subject, config = inputs(tmp_path)
    config["paths"]["atlases"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA256"):
        extract_acute(subject, config, tmp_path / "out")


def test_uncovered_is_union_not_sum_of_atlas_overlaps(tmp_path):
    subject, config = inputs(tmp_path)
    config["paths"]["atlases"][0]["include_labels"] = [1]
    extract_acute(subject, config, tmp_path / "out")
    row = read_output(tmp_path / "out" / "acute.tsv")
    assert float(row["acute_uncovered_fraction"]) == .5
    assert "acute_test_2" not in row


def test_grid_resampling_preserves_original_lesion_volume(tmp_path):
    subject, config = inputs(tmp_path)
    atlas = config["paths"]["atlases"][0]
    highres = np.ones((8, 8, 8), dtype=np.uint8)
    highres[4:] = 2
    save_image(highres, np.eye(4), atlas["path"])
    extract_acute(subject, config, tmp_path / "out")
    row = read_output(tmp_path / "out" / "acute.tsv")
    assert float(row["acute_volume_ml"]) == pytest.approx(.016)
    assert float(row["acute_test_1"]) == pytest.approx(1 / 32)
    config["analysis"]["acute"]["allow_grid_resampling"] = False
    with pytest.raises(ValueError, match="resampling is disabled"):
        extract_acute(subject, config, tmp_path / "out_disabled")


@pytest.mark.parametrize("unit", ["meter", "micron", "unknown"])
@pytest.mark.parametrize("target", ["lesion", "atlas"])
def test_non_mm_or_unknown_spatial_units_rejected(tmp_path, unit, target):
    subject, config = inputs(tmp_path)
    path = subject["lesion_mask"] if target == "lesion" else config["paths"]["atlases"][0]["path"]
    image = nib.load(path)
    save_image(image.get_fdata(), image.affine, path, spatial_unit=unit)
    with pytest.raises(ValueError, match="spatial units must explicitly be mm"):
        extract_acute(subject, config, tmp_path / "out")
