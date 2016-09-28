# for simple shift to orered full atom file
# awk '{printf("%s%4d     %7.3f %7.3f %7.3f  1.00  0.00\n",substr($0,1,22),$6-237,$7,$8,$9)}' $argv[1] > $argv[2]
# for complete renumber of a shuffled full atom file
# awk '{if($6!=last){n++}; printf("%s%4d     %7.3f %7.3f %7.3f  1.00  0.00\n",substr($0,1,22),n,$7,$8,$9); last=$6}' $argv[1] > $argv[2]
# for simple renumbering of a CA file
# awk '/ATOM/ {n++; printf("ATOM %6d  CA  %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", n,$4,n,$7,$8,$9,$10,$11)}' $argv[1] > $argv[2]

awk '/ATOM/ {n++; s=0.5; printf("ATOM %6d  CA  %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", n,$4,n/10,s*$7,s*$8,s*$9,$10,$11)}' final.out > shrink.out
