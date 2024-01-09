import hashlib
import os
import subprocess
import yaml

def read_yaml(filepath = 'tests/md5sums.yml'):
    with open(filepath, 'r') as infile:
        yaml_dict = yaml.load(infile.read(), Loader=yaml.Loader)
    return yaml_dict

def get_md5(filepath, encoding = 'utf-8'):
    with open(filepath, 'r') as file:
        md5sum = hashlib.md5(file.read().encode(encoding)).hexdigest()
    return md5sum

def check_snapshot_md5sums(files):
    all_equal = True
    for snapshot in files:
        observed_md5 = get_md5(snapshot['path'])
        expected_md5 = snapshot['md5sum']
        if observed_md5 != expected_md5:
            all_equal = False
            print("md5sum changed for", snapshot['path'], '\n\texpected:', expected_md5, "observed:", observed_md5)
    return all_equal

def run_cli(command):
    out = subprocess.run(
        command,
        shell=True,
        text=True,
        capture_output=True,
    )
    if out.returncode:
        print(out.stdout)
        print(out.stderr)
    return out

def test_cli():
    mimseq_out = run_cli(
        "mimseq --species Hsap --cluster-id 0.97 --threads 2 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir test_cli --max-multi 4 --remap --remap-mismatches 0.05 tests/data/sampleData_subset_HEKvsK562.txt"
    )
    is_success = mimseq_out.returncode == 0
    snapshot_files = read_yaml('tests/md5sums.yml')['test_cli']
    assert is_success and check_snapshot_md5sums(snapshot_files)

# def test_cli_local_modomics():
#     mimseq_out = run_cli(
#         "mimseq --species Hsap --cluster-id 0.97 --threads 2 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir test_cli_local_modomics --max-multi 4 --remap --remap-mismatches 0.05 --local-modomics tests/data/sampleData_subset_HEKvsK562.txt"
#     )
#     is_success = mimseq_out.returncode == 0
#     snapshot_files = read_yaml('tests/md5sums.yml')['test_cli_local_modomics']
#     assert is_success and check_snapshot_md5sums(snapshot_files)
