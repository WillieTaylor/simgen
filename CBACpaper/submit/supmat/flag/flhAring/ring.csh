cp head.run ring.run
awk '{r=40; a=$1*3.141/180; print "SPINX 10\n" "SPINY 25\n" "SPINZ",$1,"\n" "TRANS", r*sin(a),r*cos(a),"\n" "INPUT model15.cas"}' ang.dat >> ring.run
~/newsims/main ring.run
