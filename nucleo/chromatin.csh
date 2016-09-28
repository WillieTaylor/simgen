# 1 = number of nucleosomes

#REBOND atomid  1141  1959	1959-1141 =  818
#REBOND atomid  1988  1142      1988-1142 =  846
#				2294-1141 = 1153
#REBOND atomid  2294  3112	3112-2294 =  818
#REBOND atomid  3141  2295	3142-2295 =  846

cat head.dat

@ n = $argv[1]

echo GROUP 0 $n
@ m = 1
while ( $m < $n )
	echo "	INPUT nucseg.dat"
	@ t = -30 * $m
	echo "	TRANS  0  0  $t"
	@ t = 110 * $m
	echo "	TWIST    0  90  -10     0  90  10    $t"
	@ m++
end
echo "	INPUT nucseg.dat"
@ m = 1
@ p = 818
@ q = 846
@ r = 1153
@ a = 1141
while ( $m < $n )
	@ b = $a + $p
	echo "REBOND atomid   $a  $b"
	@ d = $a + 1
	@ c = $d + $q
	echo "REBOND atomid   $c  $d"
	@ a += $r
	@ m++
end
echo RETERM atomid   806   835
echo END
