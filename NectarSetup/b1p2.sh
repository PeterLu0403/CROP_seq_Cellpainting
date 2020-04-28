#!/bin/sh
for j in $(seq 1 64 1281)
do
  for i in $(seq $j 2 $(($j+62))); do

    # First image of the batch is our $i
    BATCH_START=$i
    # Last image is $i + 1 (2 images inclusive)
    # e.g. if first in batch is 21, last is 22 (2 images)
    # Change the 1 to match your real batch size
    BATCH_END=$((i+1))

    # Run CellProfiler for tha batch
 
    # -p specifies our batch file
    # -c disables the GUI
    # -r runs the batch immediately
    # -b do not build extensions
    # -f first image set for this batch
    # -l last image set for this batch
    # -o output folder for this batch
     cellprofiler -p /home/b1p2/Dataset/Pipelines/Analysis_hTMC_B1_20200106.cpproj -c --file-list=/home/b1p2/Dataset/hTMC_B1/filelist_2.txt -r -f ${BATCH_START} -l ${BATCH_END} -o /home/b1p2/Dataset/hTMC_B1/Output_Plate2/batch_${i}_out &
  done
  wait
done
wait

for i in $(seq 1345 2 1349); do

   # First image of the batch is our $i
   BATCH_START=$i
   # Last image is $i + 1 (2 images inclusive)
   # e.g. if first in batch is 21, last is 22 (2 images)
   # Change the 1 to match your real batch size
   BATCH_END=$((i+1))

   # Run CellProfiler for tha batch

   # -p specifies our batch file
   # -c disables the GUI
   # -r runs the batch immediately
   # -b do not build extensions
   # -f first image set for this batch
   # -l last image set for this batch
   # -o output folder for this batch
    cellprofiler -p /home/b1p2/Dataset/Pipelines/Analysis_hTMC_B1_20200106.cpproj -c --file-list=/home/b1p2/Dataset/hTMC_B1/filelist_2.txt -r -f ${BATCH_START} -l ${BATCH_END} -o /home/b1p2/Dataset/hTMC_B1/Output_Plate2/batch_${i}_out &
 done
 wait
