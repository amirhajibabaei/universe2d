

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
	else
		echo "unknown act: $act"
	fi
done
