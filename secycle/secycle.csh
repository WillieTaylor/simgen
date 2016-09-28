# 1 = pdb, 2 = beta bond error (0.3)

cp $argv[1] model.pdb
@ n = 0
while (1)
	echo Cycle $n
	sims/setsec/setsec model.pdb 1.3 1.0 $argv[2] | grep 'Aa-'
	./simgen > simgen.log
	@ n++
end
