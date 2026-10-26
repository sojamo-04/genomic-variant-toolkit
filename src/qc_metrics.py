"""Quality Control and Coverage Metrics Module for Targeted NGS Panels.

Based on HGUA SOP criteria and Ion Reporter Info.csv output parsing.
"""

from pathlib import Path
from typing import Any, Dict, Optional

# Recuento de amplicones por panel segun la configuracion del secuenciador
PANEL_AMPLICONS: Dict[str, int] = {
    "OPA": 263,  # OPA GX5 Combo DNA Fusions
    "OCA_V3": 3707,  # OCA v3 DNA Fusions
}


def parse_info_csv(
    csv_path: Path, panel_type: str = "OPA"
) -> Dict[str, Optional[float]]:
  """Parses DNA and RNA QC blocks from Genexus Info.csv file."""
  metrics: Dict[str, Optional[float]] = {
      "mapped_reads": None,
      "uniformity": None,
      "mean_depth": None,
      "mapped_fusion_reads": None,
      "mean_read_length": None,
  }

  path = Path(csv_path)
  if not path.exists():
    return metrics

  try:
    with open(path, "r", encoding="utf-8-sig") as f:
      lines = f.readlines()

    current_section: Optional[str] = None
    for line in lines:
      line_clean = line.strip()
      if line_clean == "DNA":
        current_section = "DNA"
        continue
      elif line_clean == "RNA":
        current_section = "RNA"
        continue

      if current_section == "DNA":
        if "Uniformity Of" in line_clean and "," in line_clean:
          parts = line_clean.split(",")
          if len(parts) >= 2:
            try:
              metrics["uniformity"] = float(parts[1].strip())
            except ValueError:
              pass
        elif line_clean.startswith("Mapped Reads,"):
          parts = line_clean.split(",")
          if len(parts) >= 2:
            try:
              mapped = int(parts[1].strip())
              metrics["mapped_reads"] = mapped
              amplicons = PANEL_AMPLICONS.get(panel_type.upper(), 263)
              metrics["mean_depth"] = round(mapped / amplicons)
            except ValueError:
              pass

      elif current_section == "RNA":
        if line_clean.startswith("Mapped Fusion Reads,"):
          parts = line_clean.split(",")
          if len(parts) >= 2:
            try:
              metrics["mapped_fusion_reads"] = float(parts[1].strip())
            except ValueError:
              pass
        elif "Mean Read Length" in line_clean and "," in line_clean:
          parts = line_clean.split(",")
          if len(parts) >= 2:
            try:
              metrics["mean_read_length"] = float(parts[1].strip())
            except ValueError:
              pass
  except Exception:
    pass

  return metrics


def evaluate_sample_qc(
    mapd: Optional[float],
    uniformity: Optional[float],
    mapped_fusion_reads: Optional[int] = None,
    mean_read_length: Optional[float] = None,
    tumor_cellularity: Optional[float] = None,
    evalua_arn: bool = True,
) -> Dict[str, Any]:
  """Evaluates sample QC according to HGUA SOP thresholds.

  - DNA sufficient: MAPD < 0.5 and Uniformity > 90%
  - RNA sufficient: Mapped Fusion Reads >= 200,000 and Mean Read Length >= 60
  - Cellularity: > 20%
  """
  motivos_fallo = []

  # Evaluacion ADN
  dna_ok = True
  if mapd is None or mapd >= 0.5:
    dna_ok = False
    motivos_fallo.append(f"MAPD >= 0.5 o ausente ({mapd})")
  if uniformity is None or uniformity <= 90.0:
    dna_ok = False
    motivos_fallo.append(f"Uniformidad <= 90% o ausente ({uniformity})")

  # Evaluacion Celularidad
  celularidad_ok = True
  if tumor_cellularity is not None and tumor_cellularity <= 20.0:
    celularidad_ok = False
    motivos_fallo.append(f"Celularidad <= 20% ({tumor_cellularity}%)")

  # Evaluacion ARN (si aplica al panel)
  rna_ok = True
  if evalua_arn:
    if mapped_fusion_reads is not None and mapped_fusion_reads < 200000:
      rna_ok = False
      motivos_fallo.append(
          f"Mapped Fusion Reads < 200.000 ({mapped_fusion_reads})"
      )
    if mean_read_length is not None and mean_read_length < 60.0:
      rna_ok = False
      motivos_fallo.append(f"Mean Read Length < 60 bp ({mean_read_length})")

  es_suficiente = dna_ok and celularidad_ok and (rna_ok if evalua_arn else True)

  return {
      "es_suficiente": es_suficiente,
      "dna_ok": dna_ok,
      "rna_ok": rna_ok,
      "celularidad_ok": celularidad_ok,
      "motivos_fallo": motivos_fallo,
  }
