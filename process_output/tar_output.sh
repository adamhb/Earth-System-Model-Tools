#!/bin/bash

folder_name=$1
folder_path=/glade/scratch/adamhb/archive/$1
tar_ext='.tar.gz'
file_name="${folder_name}${tar_ext}"
tar_name="${folder_path}/${file_name}"
echo $tar_name
cd $folder_path
tar -czvf $tar_name ./

