"""Audit fixed-analysis-catalogue density-recovery artifact."""
from pathlib import Path
import argparse,json
from dark_siren_robustness.fixed_analysis_catalogue_audit import audit_fixed_catalogue_artifact
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--artifact-dir",type=Path,required=True);p.add_argument("--output",type=Path)
 a=p.parse_args();r=audit_fixed_catalogue_artifact(a.artifact_dir);s=json.dumps(r,indent=2,sort_keys=True)+"\n"
 if a.output:a.output.write_text(s,encoding="utf-8")
 print(s)
if __name__=="__main__":main()
