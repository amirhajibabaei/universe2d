line=`awk '{print $4}' $1`
#
alphas=`echo $line | awk '{print $1}'`
line=`echo $line | cut -d " " -f2-`
#
tems=`echo $line | awk '{print $1}'`
line=`echo $line | cut -d " " -f2-`
#
rhos=`echo $line | awk '{print $1}'`
line=`echo $line | cut -d " " -f2-`
#
nxs=`echo $line | awk '{print $1}'`
line=`echo $line | cut -d " " -f2-`
#
cuts=`echo $line | awk '{print $1}'`
line=`echo $line | cut -d " " -f2-`
#
#echo cut=$cuts nx=$nxs rho=$rhos tem=$tems alpha=$alphas
