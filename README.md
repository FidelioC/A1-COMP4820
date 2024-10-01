# A1-COMP4820

## Experiment 1

### Result

![Experiment 1](/Experiment1.png)

### Commands for The Experiment

#### Initialize Files

```
head -47000 Sorangium_cellulosum.fasta > Sorangium_cellulosum_25percent.fasta

head -93000 Sorangium_cellulosum.fasta > Sorangium_cellulosum_50percent.fasta
```

#### Brute Force

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum_25percent.fasta --algorithm brutefoce" "python3 main.py --pattern ATGA --file Sorangium_cellulosum_50percent.fasta --algorithm brutefoce" "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm brutefoce"
```

#### BMH

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum_25percent.fasta --algorithm bmh" "python3 main.py --pattern ATGA --file Sorangium_cellulosum_50percent.fasta --algorithm bmh" "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm bmh"
```

#### BMH-Pattern-Alphabet

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum_25percent.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern ATGA --file Sorangium_cellulosum_50percent.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet"
```

#### Aho-Corasick

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum_25percent.fasta --algorithm aho-corasick" "python3 main.py --pattern ATGA --file Sorangium_cellulosum_50percent.fasta --algorithm aho-corasick" "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm aho-corasick"
```

## Experiment 2
