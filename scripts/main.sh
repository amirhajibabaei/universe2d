

# initiate variables
home=`pwd`
alias python=python2
source scripts/get_com_args.sh

# path
path=$root/c_$cut/a_$alpha/t_$tem/r_$rho/n_$nx

# functions
function lsdirs {
	for p in `ls -d $path`
	do
		echo $p
	done
}


for D in `lsdirs`
do
	source ./scripts/get_seed.sh $D/seed.txt  # defines cuts alphas rhos tems nxs
	if [ $act == "avg" ]
	then
		echo $rhos `python scripts/analyse_scalars.py $D/mc_scalars.txt`
	elif [ $act == "list" ]
	then
		echo $D
	elif [ $act == "progress" ]
	then
		if [ -f $D/mc_restart.txt ]
		then
			echo $D `sed -n "2p" $D/mc_restart.txt`  `ls $D/*.pe* | xargs -n1 basename | sed 's/.*\.pe//' `
		else
			echo $D 0
		fi
	elif [ $act == "post" ]
	then
		echo $rhos `python scripts/auto_cf.py $D/mc_scalars.txt`

	else
		echo "unknown act: $act"
	fi
done
