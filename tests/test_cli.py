import os
import subprocess

def test_cli():
    mimseq_out = subprocess.run(
        "mimseq --species Hsap --cluster-id 0.97 --threads 2 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir hg38_HEK239vsK562 --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt",
        shell=True,
        text=True,
        capture_output=True,
    )
    is_success = mimseq_out.returncode == 0
    # TODO: check md5sums of output files
    assert is_success
