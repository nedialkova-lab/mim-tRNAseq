import os
import subprocess

def run_cli(command):
    out = subprocess.run(
        command,
        shell=True,
        text=True,
        capture_output=True,
    )
    if out.returncode:
        print(out)
    return out

def test_cli():
    mimseq_out = run_cli(
        "mimseq --species Hsap --cluster-id 0.97 --threads 2 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir test_cli --max-multi 4 --remap --remap-mismatches 0.05 tests/data/sampleData_subset_HEKvsK562.txt"
    )
    is_success = mimseq_out.returncode == 0
    # TODO: check md5sums of output files
    assert is_success

def test_cli_local_modomics():
    mimseq_out = run_cli(
        "mimseq --species Hsap --cluster-id 0.97 --threads 2 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir test_cli_local_modomics --max-multi 4 --remap --remap-mismatches 0.05 --local-modomics tests/data/sampleData_subset_HEKvsK562.txt"
    )
    is_success = mimseq_out.returncode == 0
    # TODO: check md5sums of output files
    assert is_success
