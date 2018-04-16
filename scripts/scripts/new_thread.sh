
# parameters
cutoff=2.5
alpha=$1
T=$2
rho=$3
nx=$4

# paths
path=c_$cutoff/a_$alpha/t_$T/r_$rho/n_$nx
file=$path/seed.txt

# main
if [ -f $file ]; then
	echo "$file already exists"
else

	mkdir -p $path

cat << EOF > $file
variable        alpha   equal    $alpha
variable        T       equal    $T
variable        rho     equal    $rho
variable        nx      equal    $nx
variable        cutoff  equal    $cutoff
variable        dmax    equal    0.15
variable        tdamp   equal    1000
variable        dtime   equal    0.005
EOF

	echo "$file created"
	here=`pwd`
	cp mc_submit.sh $path
	cd $path
	qsub mc_submit.sh
	cd $here
	echo "returned to $pwd"


fi
