
/pggb-latest.simg /pggb \
    -i /chr${i}.fa \
    -s 5000 -l 25000 -p 90 -c 1 -K 19 -F 0.001 -g 30 \
    -k 23 -f 0 -B 10M \
    -n 111 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
    -o /output/chr${i} \
    -Y "#" --threads 20 --poa-threads 20

