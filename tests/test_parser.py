from pathlib import Path
import tempfile
from src.vcf_parser import VCFParser


def test_vcf_parser_extraction():
  sample_vcf = """##fileformat=VCFv4.1
##mapd=0.194
##TotalMappedFusionPanelReads=250000
##FusionSampleQC=PASS,DetailedQC
##IonReporterWorkflowName=OPA GX5 w3.4.2 Combo DNA Fusions
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr12\t25398284\t.\tC\tT\t100\tPASS\tAF=0.15;FDP=650;FUNC=[{'gene':'KRAS','protein':'p.Gly12Asp','coding':'c.35G>A'}]\tGT:GQ\t0/1:100
chr7\t116371891\t.\tN\t<CNV>\t.\tGAIN\tSVTYPE=CNV;END=116422166;CI=0.05:4.2,0.95:6.1;FUNC=[{'gene':'MET'}]\tCN\t5.2
"""
  with tempfile.NamedTemporaryFile("w+", suffix=".vcf", delete=False) as f:
    f.write(sample_vcf)
    temp_path = Path(f.name)

  try:
    parser = VCFParser(temp_path)
    header = parser.parse_header()
    variants = parser.parse_variants()

    assert header["mapd"] == 0.194
    assert header["fusion_qc"] == "PASS"

    assert len(variants["mutations"]) == 1
    mut = variants["mutations"][0]
    assert mut["gene"] == "KRAS"
    assert mut["protein"] == "p.Gly12Asp"
    assert mut["af"] == 0.15
    assert mut["dp"] == 650

    assert len(variants["cnvs"]) == 1
    cnv = variants["cnvs"][0]
    assert cnv["gene"] == "MET"
    assert cnv["cn"] == 5.2
  finally:
    temp_path.unlink()
