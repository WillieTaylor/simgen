rm mull.models/*
foreach pdb ( `ls -1 mullout/batch*/models/model*.cas ` )
	@ n++
	cp $pdb mull.models/model$n.cas
end
