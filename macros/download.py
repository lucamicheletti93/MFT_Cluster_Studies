import os

# Define arguments directly
target_dir = '/data/shared/mft_pid/data/test/'  # Target directory for merged files
input_dirs = ['///alice/data/2022/LHC22o/526641/apass7/0820/o2_ctf_run00526641_orbit0282004224_tf0000587293_epn226/001']  # List of input directories
suffix = ['lhc22o_pass7_test']  # List of suffixes for each input directory
file_to_merge = ['AO2D']  # List of files to merge

# Loop over the input directories and suffixes
for i, (input_dir, suf) in enumerate(zip(input_dirs, suffix)):
    train_number = input_dir.split('/')[-2]
    for file in file_to_merge:
        os.system(f'alien.py cp -T 64 alien://{input_dir}/{file}.root file:{target_dir}/{file}_{suf}.root')
