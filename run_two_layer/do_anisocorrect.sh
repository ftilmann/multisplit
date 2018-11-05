#!/bin/sh
# 
# sample usage: sh do_anisocorrect.sh D21 2007.035.20.56.59 55 0.7 120 1.2 
# sample usage: sh do_anisocorrect.sh D21 2007.035.20.56.59 1000  
#

stnm=$1
sacid=$2
mcrunno=$3


ofile1=single_indvsol.dat
ofile2=all_indvsol.dat
finalfi=final_swspar.out

finalparfi=`echo "'"$finalfi"'"`


outps1=`echo $sacid"_1.ps"`
outps2=`echo $sacid"_2.ps"`
outps3=`echo $sacid"_3.ps"`


/bin/rm -f $ofile1 $ofile2 $outps1 $outps2 $outps3

logfi=`echo $sacid".log"`		

fpdup=`grep "Top Fast dir.:" $logfi | awk '{print $4}'`
tdup=`grep "Top Fast dir.:" $logfi | awk '{print $7}'`
fpdlow=`grep "Bot Fast dir.:" $logfi | awk '{print $4}'`
tdlow=`grep "Bot Fast dir.:" $logfi | awk '{print $7}'`


echo $fpdup $tdup $fpdlow $tdlow > $ofile1

filename1=`echo "'"$ofile1"'"`


path1=`echo "stat_"$stnm"/"$sacid"*.SAC"`
sacbhe=`ls $path1 | grep ".*HE." | awk '{print $1}'`
sacbhn=`ls $path1 | grep ".*HN." | awk '{print $1}'`



bheproc=`echo $sacbhe | awk -F "/" '{print $2}'`
bhnproc=`echo $sacbhn | awk -F "/" '{print $2}'`

/bin/rm -f $bheproc $bhnproc tmprad.sac tmptan.sac

/bin/cp $sacbhe $bheproc 
/bin/cp $sacbhn $bhnproc 

# Rotates horizontal components into radial & tangential components
sac<< !
read $sacbhn $sacbhe
rotate to gcarc
w tmprad.sac tmptan.sac
!


# corrects for average upper layer
split_cor -v -data $bhnproc $bheproc -split $fpdup $tdup -pre pstr1

bheprocnew=`echo "pstr1"$bheproc`
bhnprocnew=`echo "pstr1"$bhnproc`


# corrects for average lower layer
split_cor -v -data $bhnprocnew $bheprocnew -split $fpdlow $tdlow -pre pstr2



corfile11=`ls $bheprocnew | grep "..*HE." | awk '{print $1}'`
corfile12=`ls $bhnprocnew | grep "..*HN." | awk '{print $1}'`

crfi11=`echo "'"$corfile11"'"`
crfi12=`echo "'"$corfile12"'"`

bheprocnewest=`echo "pstr2pstr1"$bheproc`
bhnprocnewest=`echo "pstr2pstr1"$bhnproc`


corfile21=`ls $bheprocnewest | grep "..*HE." | awk '{print $1}'`
corfile22=`ls $bhnprocnewest | grep "..*HN." | awk '{print $1}'`

crfi21=`echo "'"$corfile21"'"`
crfi22=`echo "'"$corfile22"'"`


sacidnm=`echo "'"$sacid"'"`


# fetch all individual solutions into a single parametric file
for elogfi in `ls 19*.log 20*.log | awk '{print $1}'`; do

efpdup=`grep "Top Fast dir.:" $elogfi | awk '{print $4}'`
etdup=`grep "Top Fast dir.:" $elogfi | awk '{print $7}'`
efpdlow=`grep "Bot Fast dir.:" $elogfi | awk '{print $4}'`
etdlow=`grep "Bot Fast dir.:" $elogfi | awk '{print $7}'`


echo $efpdup $etdup $efpdlow $etdlow >> $ofile2
filename2=`echo "'"$ofile2"'"`
done

#matlab -nodesktop -nosplash -r "[data]=do_plot_vars( $mcrunno ) ; quit; "

/bin/rm -f tmp1.ps tmp2.ps tmp3.ps


matlab -nodesktop -nosplash -r "[data]=do_plot_vars( $mcrunno , $finalparfi , $sacidnm , $filename1 , $filename2 , $crfi11 , $crfi12 , $crfi21 , $crfi22, 'tmprad.sac' , 'tmptan.sac') ; quit; " 

/bin/mv tmp1.ps $outps1
/bin/mv tmp2.ps $outps2
/bin/mv tmp3.ps $outps3

/bin/rm -f $bheproc $bhnproc $bheprocnew $bhnprocnew $bheprocnewest $bhnprocnewest
