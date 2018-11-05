#!/bin/sh
#
#
# sample usage:  sh do_plot_forall.sh D21 1000

stnm=$1
mcrun=$2

# fetch all individual solution IDs and use as an input for do_anisocorrect.sh
for fileID in `ls 19*.log 20*.log  | awk -F ".log" '{print $1}'`; do

path1=`echo $fileID"*.SAC"`
/bin/rm -f pstr*SAC $path1


sh do_anisocorrect.sh $stnm $fileID $mcrun

done

