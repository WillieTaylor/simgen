# 1 = speed

make
rm run.log

# 50:20 + subunit
echo >> run.log
cp bump1.model temp.model
sed -i 's/HH/50/g' temp.model
sed -i 's/SS/20/g' temp.model
echo "hard:soft = 50:20" >> run.log
sed -i 's/CC/ 0/' temp.model
echo "subunit" >> run.log
# no NC link
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/50/g' temp.run
	echo "miss    " >> tmp.log
	./simprot temp.run $argv[1] > sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/25/g' temp.run
	echo "clip    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/15/g' temp.run
	echo "half    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/ 5/g' temp.run
	echo "full    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
# with NC link
echo >> run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/50/g' temp.run
	echo "miss    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/25/g' temp.run
	echo "clip    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/15/g' temp.run
	echo "half    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/ 5/g' temp.run
	echo "full    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log

# 50:20 + domains
echo >> run.log
cp bump1.model temp.model
sed -i 's/HH/50/g' temp.model
sed -i 's/SS/20/g' temp.model
echo "hard:soft = 50:20" >> run.log
sed -i 's/CC/ 1/' temp.model
echo "domains" >> run.log
# no NC link
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/50/g' temp.run
	echo "miss    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/25/g' temp.run
	echo "clip    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/15/g' temp.run
	echo "half    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	sed -i 's/LINK/#LINK/g' temp.run
	echo "	no link  " >> tmp.log 
	sed -i 's/YY/ 5/g' temp.run
	echo "full    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
# with NC link
echo >> run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/50/g' temp.run
	echo "miss    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/25/g' temp.run
	echo "clip    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/15/g' temp.run
	echo "half    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
	rm tmp.log
	cp bump1.run temp.run
	echo "	link on  " >> tmp.log 
	sed -i 's/YY/ 5/g' temp.run
	echo "full    " >> tmp.log
	./simprot temp.run $argv[1] >sim.log
	grep aa sim.log | sort -nr -k12 | head -1 >> tmp.log
	cat tmp.log | tr -d "\n" >> run.log
	tcsh domrms.csh >> run.log
tail -1 run.log
