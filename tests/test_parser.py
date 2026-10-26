import tempfile
from pathlib import Path
from src.vcf_parser import VCFParser
from src.qc_metrics import evaluate_sample_qc, parse_info_csv


def test_vcf_parser_with_mock_file():
    """Valida el parseo completo de cabecera y variantes usando data/mock_sample.vcf."""
    mock_path = Path(__file__).parent.parent / "data" / "mock_sample.vcf"
    assert mock_path.exists(), f"No se encuentra el archivo: {mock_path}"

    parser = VCFParser(mock_path)
    header = parser.parse_header()
    variants = parser.parse_variants()

    # QC de cabecera
    assert header["mapd"] == 0.194
    assert header["fusion_qc"] == "PASS"
    assert header["total_fusion_reads"] == 254000
    assert header["workflow"] == "OPA GX5 w3.4.2 Combo DNA Fusions"

    # SNVs/indels
    assert len(variants["mutations"]) == 1
    mut = variants["mutations"][0]
    assert mut["gene"] == "KRAS"
    assert mut["protein"] == "p.Gly12Asp"
    assert mut["coding"] == "c.35G>A"
    assert mut["af"] == 0.15
    assert mut["dp"] == 650
    assert mut["filter"] == "PASS"

    # CNVs
    assert len(variants["cnvs"]) == 1
    cnv = variants["cnvs"][0]
    assert cnv["gene"] == "MET"
    assert cnv["cn"] == 5.2
    assert cnv["filter"] == "GAIN"

    # Alteraciones de ARN / Fusiones
    assert len(variants["rna"]) == 1
    rna = variants["rna"][0]
    assert rna["svtype"] == "Fusion"
    assert rna["read_count"] == 320


def test_evaluate_sample_qc_pass():
    """Valida una muestra con parámetros de calidad suficientes según el SOP[cite: 3, 6]."""
    qc_result = evaluate_sample_qc(
        mapd=0.18,
        uniformity=95.5,
        mapped_fusion_reads=220000,
        mean_read_length=75.0,
        tumor_cellularity=50.0,
        evalua_arn=True,
    )
    assert qc_result["es_suficiente"] is True
    assert qc_result["dna_ok"] is True
    assert qc_result["rna_ok"] is True
    assert qc_result["celularidad_ok"] is True
    assert len(qc_result["motivos_fallo"]) == 0


def test_evaluate_sample_qc_fails():
    """Valida que se detecten y listen los fallos de calidad (MAPD, ARN y celularidad)[cite: 3, 6]."""
    qc_result = evaluate_sample_qc(
        mapd=0.55,  # Fallo: MAPD >= 0.5[cite: 3, 6]
        uniformity=88.0,  # Fallo: Uniformidad <= 90%[cite: 3, 6]
        mapped_fusion_reads=150000,  # Fallo: < 200.000 lecturas[cite: 3, 6]
        mean_read_length=50.0,  # Fallo: < 60 bp[cite: 3, 6]
        tumor_cellularity=15.0,  # Fallo: <= 20%[cite: 3, 6]
        evalua_arn=True,
    )
    assert qc_result["es_suficiente"] is False
    assert qc_result["dna_ok"] is False
    assert qc_result["rna_ok"] is False
    assert qc_result["celularidad_ok"] is False
    assert len(qc_result["motivos_fallo"]) == 5


def test_parse_info_csv_depth_calculation():
    """Valida el cálculo de profundidad media por amplicones en OPA (263) y OCA v3 (3707)[cite: 4, 5]."""
    mock_csv_content = """DNA
Mapped Reads,526000
Uniformity Of Base Coverage,98.5
RNA
Mapped Fusion Reads,210000
Mean Read Length (bp),72.4
"""
    with tempfile.NamedTemporaryFile("w+", suffix=".csv", delete=False) as f:
        f.write(mock_csv_content)
        temp_csv = Path(f.name)

    try:
        # Panel OPA: 526000 / 263 = 2000[cite: 4]
        metrics_opa = parse_info_csv(temp_csv, panel_type="OPA")
        assert metrics_opa["mapped_reads"] == 526000
        assert metrics_opa["uniformity"] == 98.5
        assert metrics_opa["mean_depth"] == 2000
        assert metrics_opa["mapped_fusion_reads"] == 210000
        assert metrics_opa["mean_read_length"] == 72.4

        # Panel OCA v3: 526000 / 3707 ≈ 142[cite: 5]
        metrics_oca = parse_info_csv(temp_csv, panel_type="OCA_V3")
        assert metrics_oca["mean_depth"] == round(526000 / 3707)
    finally:
        temp_csv.unlink()
