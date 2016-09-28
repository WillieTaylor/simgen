foreach n ( `awk '{print $1}' $argv[1]` )
	@ len = `~/util/topol models/model$n.cas | grep len | tail -1 | awk '{print $3}'`
	echo model$n.cas $len
	echo $len >> knot.len
end
paste -d ' ' rms.dat knot.len | grep -v ' 2$' > knot.dat
