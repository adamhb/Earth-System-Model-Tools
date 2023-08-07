#!/bin/bash


start_year=1900
end_year=2049

case_1=CZ2-bench-1951-1979_CLM-17e2acb6a_FATES-3f6effb3_frequent_fire
case_2=CZ1-bench-1951-1979_CLM-17e2acb6a_FATES-3f6effb3_frequent_fire
case_3=stan-bench-1951-1979_CLM-17e2acb6a_FATES-3f6effb3_frequent_fire

python get_output.py -c $case_1 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
tar -czf "${case_1}.tar.gz" -C dir_path .

python get_output.py -c $case_2 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
tar -czf "${case_2}.tar.gz" -C dir_path .

python get_output.py -c $case_3 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
tar -czf "${case_3}.tar.gz" -C dir_path .








