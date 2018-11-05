#!/bin/sh 


stnm=$1
mrun=$2

datadir=`echo "stat_"$stnm`

path1=`echo $datadir"/*SAC"`

/bin/rm -f $outfile *.hdr *.bin *.log *_mc.x *_bootstrap.x 


for sacID in `ls $path1 | awk -F "00.BH" '{print $1}' | sort -u`; do


sace=`echo $sacID"00.BHE.M.SAC"`
sacn=`echo $sacID"00.BHN.M.SAC"`

outID=`echo $sacID | awk -F "/" '{print $2}' | awk -F "." '{print $1"."$2"."$3"."$4"."$5""}'`
outfile=`echo $outID".log"`


echo "Double layer analysis for "$sacn $sace 

#multisplit -mt -double 1 0.05 5.0 -name $outID -data $sacn $sace -taper 1.0 -wint t4 -10 25 -dof 0.01 
##multisplit -mintransverse -data $sacn $sace -dof 0.01 -taper 5 -gmt5 -name $outID -wint t4 -10 25 -double 5 0.05 3.0 > $outfile
multisplit -mintransverse -data $sacn $sace -dof 0.01 -taper 5 -gmt5 -name $outID -wint t4 -15 25 -double 5 0.05 3.0 > $outfile


done

# error_stack single splitting
#
#error_stack -name single *_err.bin
#error_stack -name double *2l_err2.bin


# Generates distribution based on Monte Carlo from stacked error
#error_stack -mc 1000 -name double *_err2.bin

# Generates distribution based on Bootstrap from stacked error
/bin/rm -f screen.out final_swspar.out 
error_stack -bootstrap $mrun -name double *_err2.bin > screen.out

FPDupfinal=`grep "TopFast:" screen.out | tail -1 | awk '{print $2}'`

TDupfinal=`grep "TopSplitDly:" screen.out | tail -1 | awk '{print $2}'`

FPDbotfinal=`grep "BotFast:" screen.out | tail -1 | awk '{print $2}'`

TDbotfinal=`grep "BotSplitDly:" screen.out | tail -1 | awk '{print $2}'`

echo $FPDupfinal $TDupfinal $FPDbotfinal $TDbotfinal > final_swspar.out

