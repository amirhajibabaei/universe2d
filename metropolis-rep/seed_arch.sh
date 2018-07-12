
# parameters
pow=12
cutoff=1.8
rho=$1
nx=256

# paths
path=p_$pow/c_$cutoff/r_$rho/n_$nx
file=$path/seed.txt

# main
if [ -f $file ]; then
	echo "$file already exists"
else

	mkdir -p $path

cat << EOF > $file
variable        pow     equal    $pow
variable        cutoff  equal    $cutoff
variable        rho     equal    $rho
variable        nx      equal    $nx
variable        dmax    equal    0.3
variable        tdamp   equal    1000
EOF

	echo "$file created"
	here=`pwd`
	cp seed_submit.sh $path
	cd $path
	qsub seed_submit.sh
	cd $here
	echo "returned to $pwd"

fi
