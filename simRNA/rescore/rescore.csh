# 1 = xxx(.model) dir, 2 = ideal dist, 3 = give, 4 = N links, 5 = links.dat
#
# tcsh rescore.csh sams 8 4 50 gremlin.dat

#gnuplot> plot 'rescore.true' u 2:3 not, 'rescore.fake' u 2:3 not
 
if ( -e rescore.dat ) rm rescore.dat
head -$argv[4] $argv[5] > links.tmp

# score fake models
@ n = 1
foreach model (`ls $argv[1].models/model*.cas`)
	home/scoreWrms links.tmp $model $argv[2] $argv[3] > score.log
	set bum = `head -1 score.log | awk '{print $3}'`
	set rms = `tail -1 score.log | awk '{print $6}'`
	set sco = `tail -2 score.log | head -1 | awk '{print $12}'`
	echo $model $rms $sco $bum >> rescore.dat
	@ n++
end
set fake = `awk '{n++; s+=$3; print s/n}' rescore.dat | tail -1`
set vari = `awk -v m=$fake '{n++; s=$3-m; ss+=s*s; print ss/n}' rescore.dat | tail -1`
mv rescore.dat rescore.fake
# score true models
@ n = 1
foreach model (`ls true.models/model*.cas`)
	home/scoreWrms links.tmp $model $argv[2] $argv[3] > score.log
	set bum = `head -1 score.log | awk '{print $3}'`
	set rms = `tail -1 score.log | awk '{print $6}'`
	set sco = `tail -2 score.log | head -1 | awk '{print $12}'`
	echo $model $rms $sco $bum >> rescore.dat
	@ n++
end
set true = `awk '{n++; s+=$3; print s/n}' rescore.dat | tail -1`
mv rescore.dat rescore.true
# calc Z score
echo $true $fake $vari | awk '{print ($1-$2)/sqrt($3)}'
echo ranked by
echo raw score
awk '{print $0, $3}' rescore.fake | sort -nr -k5 |  head
echo "$3*$3/($4+10)"
awk '{print $0, $3*$3/($4+10)}' rescore.fake | sort -nr -k5 |  head
echo "$3*$3-($4*10)"
awk '{print $0, $3*$3-($4*10)}' rescore.fake | sort -nr -k5 |  head
