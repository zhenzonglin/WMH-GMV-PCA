"""图谱准备的标签编号、网络策略、半球核查及不可覆盖规则。"""
import importlib.util
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
import requests


script_path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_atlases.py"
spec = importlib.util.spec_from_file_location("prepare_atlases", script_path)
atlas = importlib.util.module_from_spec(spec)
spec.loader.exec_module(atlas)


def test_xml_offsets_for_labels_and_probability_volumes(tmp_path):
    source = tmp_path / "labels.xml"
    source.write_text('<atlas><data><label index="0">A</label><label index="1">B</label></data></atlas>')
    assert atlas.xml_labels(source, 1, True) == {0: "Background", 1: "A", 2: "B"}
    assert atlas.xml_labels(source, 0, False) == {0: "A", 1: "B"}


def test_duplicate_xml_indices_rejected(tmp_path):
    source = tmp_path / "labels.xml"
    source.write_text('<atlas><data><label index="0">A</label><label index="0">B</label></data></atlas>')
    with pytest.raises(ValueError, match="Repeated label"):
        atlas.xml_labels(source, 1, True)


def test_offline_policy_rejects_http_before_connection():
    with atlas.network_policy(False, 1):
        with pytest.raises(RuntimeError, match="downloads disabled"):
            requests.get("https://example.com/")


def test_lateralization_uses_centroid_and_reports_published_boundary_crossing():
    values = np.zeros((101, 1, 1), dtype=np.int16)
    values[:51] = 1
    affine = np.eye(4)
    affine[0, 3] = -49
    image = nib.Nifti1Image(values, affine)
    # 左侧区域允许少量图谱边界体素越过中线；不修改图谱定义。
    result = atlas.validate_left_right(image, {0: "Background", 1: "Left test"}, 1e-5, .05)
    assert result["lateralization_qc"][0]["opposite_hemisphere_voxels"] == 1
    assert result["lateralization_qc"][0]["opposite_hemisphere_fraction"] == pytest.approx(1/51)
    with pytest.raises(ValueError, match="violates world-coordinate"):
        atlas.validate_left_right(image, {0: "Background", 1: "Right test"}, 1e-5, .05)


def test_publish_is_idempotent_and_does_not_replace_different_resource(tmp_path):
    source, output = tmp_path / "source", tmp_path / "output"
    source.write_bytes(b"atlas")
    assert atlas.publish(source, output, False) == "created"
    stamp = output.stat().st_mtime_ns
    assert atlas.publish(source, output, False) == "unchanged"
    assert output.stat().st_mtime_ns == stamp
    source.write_bytes(b"different")
    with pytest.raises(FileExistsError):
        atlas.publish(source, output, False)
    assert output.read_bytes() == b"atlas"
