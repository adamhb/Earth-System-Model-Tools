#!/bin/bash


start_year=1901
end_year=1970

case_1=CZ2-spinup-1951-1979-HF-051823_CLM-17e2acb6a_FATES-d738ae07
case_2=stan-spinup-1951-1979-fire-suppression-051823_CLM-17e2acb6a_FATES-d738ae07
case_3=CZ2-spinup-1951-1979-fire-suppression-051823_CLM-17e2acb6a_FATES-d738ae07
case_4=stan-spinup-1951-1979-HF-051823_CLM-17e2acb6a_FATES-d738ae07


python get_output.py -c $case_1 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
#dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
#tar -czf "${case_1}.tar.gz" -C dir_path .
echo "done with case 1"

python get_output.py -c $case_2 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
#dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
#tar -czf "${case_2}.tar.gz" -C dir_path .
echo "done with case 2"

python get_output.py -c $case_3 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
#dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
#tar -czf "${case_3}.tar.gz" -C dir_path .

echo "done with case 3"

python get_output.py -c $case_4 -s ${start_year} -e ${end_year} --dbh-min 0 --pft-names True
#dir_path=$(ls -t /glade/scratch/adamhb/archive | head -n 1)
#tar -czf "${case_3}.tar.gz" -C dir_path .

echo "done with case 4"




