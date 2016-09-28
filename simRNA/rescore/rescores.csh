# 1 = models dir, 2 = links.dat

rm rescores.$argv[2]
foreach d ( 7 8 9 10 11 12 )
	foreach s ( 2 4 6 8 10 12 15 20 )
		foreach n ( 10 20 30 40 50 60 70 80 90 100 110 120 )
			set gap = `tcsh rescore.csh $argv[1] $d $s $n $argv[2]`
			echo $d $s $n $gap >> rescores.$argv[2]
		end
		echo >> rescores.$argv[2]
	end
end
