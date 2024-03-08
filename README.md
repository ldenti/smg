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
* python3
* python3-intervaltree

### Example
```
cd example
tar xvfz 4.tar.gz
vg construct --progress --threads 4 --reference 4.fa --vcf 4.vcf.gz > 4.pg
vg rna --progress --threads 4 --add-ref-paths --transcripts 4.gtf 4.pg | vg mod --unchop - | vg ids -s - > 4.spliced.pg
vg view 4.spliced.pg > 4.spliced.gfa
python3 ../prune_graph.py 4.spliced.gfa > 4.spliced.pruned.gfa
../build/smg 4.spliced.pruned.gfa
```

### TODO
