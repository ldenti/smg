# SMG

### Compilation
```
git clone https://github.com/ldenti/smg.git
cd smg
mkdir build ; cd build
cmake ..
make
cd ..
```

### Prerequisites
* [vg](https://github.com/vgteam/vg)
* python3 (+ intervaltree seaborn pandas matplotlib)


### Example
```
cd example
tar xvfz 4.tar.gz
vg construct --progress --threads 4 --reference 4.fa --vcf 4.vcf.gz > 4.pg
vg rna --progress --threads 4 --add-ref-paths --transcripts 4.gtf 4.pg | vg mod --unchop - | vg ids -s - > 4.spliced.pg
vg view 4.spliced.pg > 4.spliced.gfa
python3 ../scripts/reduce_graph.py 4.spliced.gfa > 4.spliced.reduced.gfa
../build/smg 4.spliced.reduced.gfa > minimizers.txt
python3 ../cluster.py minimizers.txt 4.spliced.pruned.gfa > tree.txt
# python3 ../scripts/analyze_mmers.py minimizers.txt
```

### TODO
