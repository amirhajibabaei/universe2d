#!/bin/bash
#
#    -R=root  -c=cutoff  -a=alpha  -t=tem  -n=nx  -r=rho -act=act
#
for i in "$@"
do
	case $i in 
		-R=*|--root=*)
		root="${i#*=}"
		shift
		;;
	 	-c=*|--cutoff=*)
		cut="${i#*=}"
		shift
		;;
		-a=*|--alpha=*)
		alpha="${i#*=}"
		shift
		;;
		-t=*|--tem=*)
		tem="${i#*=}"
		shift	# past argument = value
		;;
		-n=*|--nx=*)
		nx="${i#*=}"
		shift	# past argument = value
		;;
		-r=*|--rho=*)
		rho="${i#*=}"
		shift	# past argument = value
		;;
		-act=*|--action=*)
		act="${i#*=}"
		shift	# past argument = value
		;;
		--default)
		DEFAULT=YES 
		shift 	# past argument no value
		;;
		*)
		echo "ERROR: unknown option"; exit 1 
		;;
	esac
done

#echo root=$root cut=$cut alpha=$alpha tem=$tem nx=$nx rho=$rho 
