## Experimental notes on developing code to predict off target primer binding

# intial steps

Initially i explored methods for finding seeds tho start the search, and just ises hammin/levenschtein distance adter finding the targets

A  few different approaches wer taken `main.rs.7mer, main.rs.pidgenhole`. I settled on using the aho_corasick multi string  tatching methods with a kmer size of 7 to find seed regions.  the 7mer script could run the human T2T genome in about 30 seconds with computationally expensive leven distance.

# first thermodynamic scripts

That inital method was used to create a simple script to  find regions over a deltaG threshold.  Initlaly i modified the main.rs.7mer  method to  then calulate  the deltad and Tm using a simplified santalicua method that neglected salt ajustments and end correction. 

model | primers |  runtime | cores | hits | command
------|---------|----------|-------|------|--------
`main_simpletest` | 1 | 4.594 s | 12 | 8493 | time time cargo run --release --bin main_simplesalt -- --file ~/Downloads/ncbi_dataset\ \(5\)/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --threshold=-20.0 --patterns patterns.txt > ss1.txt
`main_simpletest` | 10 | 11.176 s | 12 | 856620 |  time time cargo run --release --bin main_simplesalt -- --file ~/Downloads/ncbi_dataset\ \(5\)/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --threshold=-20.0 --patterns patterns10.txt > ss2.txt
`main_simpletest` | 100 | 125.12 s | 12 |  12899133 | time time cargo run --release --bin main_simplesalt -- --file ~/Downloads/ncbi_dataset\ \(5\)/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --threshold=-20.0 --patterns patterns100.txt > ss3.txt

I reran writing to dev/null t osee if writing was a big part of the speed. 

`time time cargo run --release --bin main_simplesalt -- --file ~/Downloads/ncbi_dataset\ \(5\)/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --threshold=-20.0 --patterns patterns100.txt > /dev/null`

skipping file writing reduced runtime from 125 s to 101 s.


I then created a baseline brute force script to benchmark against `main_bruteforce_simplesalt`.  After 10 minutes of running wit ha single probe is stopped the incomplete run/

the simplesalt method is a  sinplification from the full santalucia model which accounts for salt concentrations using either the Santalucie 2004 or the Oransku 2008 model.  it also skipped a  modification for gase pairs near the end of the sequence.

in one sample sequence and its reverse complement the main_simplesalt.rs score was:

`pattern_1rc     0       -21.13  58.99   CGATCGATCGATCGATCGAT`


while primer3_py predicted:

```
In [1]: import primer3
In [2]: primer3.calc_heterodimer("ATCGATCGATCGATCGATCG", "CGATCGATCGATCGATCGAT")
Out[2]: ThermoResult(structure_found=True, tm=57.29, dg=-21065.89, dh=-160400.00, ds=-449.25)
```

Close but not the same (units are cal/mole for primer3-py and kcal/mol for my script).


     Running `target/release/main_fullsalt --file pattern_rc.fna --threshold=-1.0 --patterns patterns.txt`
pattern_1rc     0       -21.13  58.99   CGATCGATCGATCGATCGAT
        0.46 real         0.04 user         0.04 sys
time cargo run --release --bin main_fullsalt -- --file pattern_rc.fna     0.05s user 0.05s system 20% cpu 0.471 total