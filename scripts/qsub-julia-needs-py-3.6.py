# NB, needs Python version >= 3.6

import os
from os.path import exists, join, basename
import argparse
from glob import glob
import subprocess


parser = argparse.ArgumentParser()

parser.add_argument('-p', '--program',
                    default='makenema.jl',
                    help='The Julia program to be executed')

parser.add_argument('-d', '--dir-in',
                    help='Directory containing input files')

parser.add_argument('-o', '--dir-out',
                    help='Output directory')

parser.add_argument('-n', '--number-of-jobs', type=int,
                    help='Split work among N jobs', required=True)

args, other_args = parser.parse_known_args()
other_args = ' '.join(other_args)

n_files = len(glob(f'{args.dir_in}/*')) # Need python >= 3.6

launch_dir = os.getcwd()

number_of_jobs = args.number_of_jobs
small_job_size, number_of_big_jobs = divmod(n_files, number_of_jobs)
number_of_small_jobs = number_of_jobs - number_of_big_jobs
first_file_in_job = 1
job_file_ranges = []
for _ in range(number_of_big_jobs):
    job_file_ranges.append((first_file_in_job, first_file_in_job + small_job_size))
    first_file_in_job += small_job_size + 1
for _ in range(number_of_small_jobs):
    job_file_ranges.append((first_file_in_job, first_file_in_job + small_job_size - 1))
    first_file_in_job += small_job_size

width = len(str(job_file_ranges[-1][-1]))

template = """#!/bin/bash
#PBS -N {basename_dir_out}-{first}-{last}
#PBS -l nodes=1:ppn=1
#cd $PBS_O_WORKDIR
cd {launch_dir}

/software/julia-1.6.1/bin/julia {program} -d {dir_in} -o {dir_out} -x evtdf-{first:0{width}}-{last:0{width}}.csv -i {first} -l {last} {other_args}
"""

output_directory = args.dir_out
version_counter = 0
while exists(output_directory):
    output_directory = f'{args.dir_out}-{version_counter}'
    version_counter += 1
os.mkdir(     output_directory)
os.mkdir(join(output_directory, 'qsub'))
os.chdir(join(output_directory, 'qsub'))

for i, l in job_file_ranges:
    qsub_script_name = f'qsub-{i:0{width}}-{l:0{width}}.sh'
    qsub_script = template.format(first=i, last=l,
                                  dir_in           = args.dir_in,
                                  dir_out          = output_directory,
                                  launch_dir       = launch_dir,
                                  basename_dir_out = basename(output_directory),
                                  width            = width,
                                  program          = args.program,
                                  other_args       = other_args,
                                  )
    with open(qsub_script_name, 'w') as file:
        file.write(qsub_script)
    print(join(os.getcwd(), qsub_script_name))
    subprocess.call(f'qsub -l mem=6gb -q long {qsub_script_name}', shell=True)
