#!/bin/sh

stnm=$1

ofile=baz_anisopar.dat

/bin/rm -f $ofile

datadir=`echo "stat_"$stnm`


# fetch all individual solutions into a single parametric file
for elogfi in `ls 19*.log 20*.log | awk '{print $1}'`; do

sacid=`echo $elogfi | awk -F ".log" '{print $1}'`
path1=`echo $datadir"/"$sacid"*.BHE.*SAC"`
sacfi=`ls $path1 | grep "BHE" | awk '{print $1}'`

/bin/rm -f tmp.out

sac << eof
read $sacfi
setbb bazi &1,baz
getbb bazi

sc echo %bazi% > tmp.out
eof
bazv=`awk '{print $1}' tmp.out`

efpdup=`grep "Top Fast dir.:" $elogfi | awk '{print $4}'`
etdup=`grep "Top Fast dir.:" $elogfi | awk '{print $7}'`
efpdlow=`grep "Bot Fast dir.:" $elogfi | awk '{print $4}'`
etdlow=`grep "Bot Fast dir.:" $elogfi | awk '{print $7}'`


echo $efpdup $etdup $efpdlow $etdlow $bazv  >> $ofile
filename2=`echo "'"$ofile"'"`
done




/bin/rm -f tmp.out
