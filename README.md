# A1-COMP4820

## Experiment 1

### Result

![Experiment 1](/Experiment1.png)

### Experiment 1 Q&A

`How does the performance of your algorithm implementations change as you increase the size of the text
T? You should evaluate for at least the following fractions of lines of the original genome file: 25%, 50%, and 100%.
Does the performance of your implementations match the expected worst-case performance for the algorithm as you increase the size of the text
T? Why or why not?`

- According from the result above, the performance of each algorithm will be slower as we increase the size of the text. This happened because each algorithm will have to iterate through a bigger text, thus consuming more time.

- Yes, the implementations of each algorithm matches the expected worst-case performance.
  - We know that the performance of Brute-force is O(n<sup>2</sup>),from the graph, we see that the time increases exponentially as we increase the text, therefore, it matches the expected worst-case performance.
  - For this experiment, both BMH and BMH-Pattern-Alphabet, has a similar runtime benchmark. And as we know, the worst case is for BMH is O(m x n). Seeing from the graph, this is true, since the time increases as we increase the length of the text.
  - For Aho-corasick, the worst case for this algorithm is O(m + n). Here, compared to the other algorithms, since the worst is O(m + n), we see that Aho-corasick has the fastest run time. And yes, this performance does match the expected worst case performance because as we increase the length of T, the time doesn't increase that much based on the graph above.

### Commands Used for Experiment 1

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

### Result

![Experiment 1](/Experiment1.png)

### Experiment 2 Q&A

`How does the performance of your algorithm implementations change as you increase the number of patterns in the set? You should evaluate for at least the following number of patterns: 1, 5, 10.
Does the performance of your implementations match the expected worst-case performance for the algorithm as you increase the number of patterns? Why or why not?`

#### Brute-force

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm brutefoce" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --file Sorangium_cellulosum.fasta --algorithm brutefoce" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --pattern ATAT --pattern AAAA --pattern GGGG --pattern CCCC --pattern TTTT --file Sorangium_cellulosum.fasta --algorithm brutefoce"
```

#### BMH

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm bmh" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --file Sorangium_cellulosum.fasta --algorithm bmh" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --pattern ATAT --pattern AAAA --pattern GGGG --pattern CCCC --pattern TTTT --file Sorangium_cellulosum.fasta --algorithm bmh"
```

#### BMH-Pattern-Alphabet

```
hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --pattern ATAT --pattern AAAA --pattern GGGG --pattern CCCC --pattern TTTT --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet"
```
