# 1 = xxx.model dir, 2 = ideal, 3 = give, 4 = N links, 5 = links.dat
 
head -$argv[4] $argv[5] > links.tmp
@ n = 0
rm rescore.dat
cp true.pdb $argv[1].models/model0.cas
foreach model (`ls $argv[1].models/model*.cas`)
	./scoreWrms links.tmp $model $argv[2] $argv[3] > score.log
	set bum = `head -1 score.log | awk '{print $3}'`
	set rms = `tail -1 score.log | awk '{print $6}'`
	set sco = `tail -2 score.log | head -1 | awk '{print $12}'`
	echo $n $rms $sco $bum >> rescore.dat
	@ n++
end
grep '^0 ' rescore.dat > true.score.dat
set true = `awk '{print $3}' true.score.dat`
set mean = `awk '{n++; s+=$3; print s/n}' rescore.dat | tail -1`
set vari = `awk -v m=$mean '{n++; s=$3-m; ss+=s*s; print ss/n}' rescore.dat | tail -1`
echo $true $mean $vari | awk '{print ($1-$2)/sqrt($3)}'
