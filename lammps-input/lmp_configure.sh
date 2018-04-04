
#
if [ -f md_restart ]; then
	echo "read_restart    md_restart" > lmp_temp
elif [ -f md_restart_alt ]; then
	echo "read_restart    md_restart_alt" > lmp_temp
else
	cat lmp_box > lmp_temp
fi
#
if [ -f md_vectors.txt ]; then
	cat md_vectors.txt >> md_vectors_hld.txt
fi
#
if [ -f md_snaps.txt ]; then
	cat md_snaps.txt >> md_snaps_hld.txt
fi
#
exit 0
