echo speed 1

echo 50:20
tcsh bump1.50:20.csh 1
cp run.log run1.50:20.log

echo 100:20
tcsh bump1.100:20.csh 1
cp run.log run1.100:20.log


echo speed 5

echo 50:20
tcsh bump1.50:20.csh 5
cp run.log run5.50:20.log

echo 100:20
tcsh bump1.100:20.csh 5
cp run.log run5.100:20.log
echo speed 5


echo speed 10

echo 50:20
tcsh bump1.50:20.csh 10
cp run.log run10.50:20.log

echo 100:20
tcsh bump1.100:20.csh 10
cp run.log run10.100:20.log
