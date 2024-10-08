# A1-COMP4820

## How to Run

`python3 main.py --file Sorangium_cellulosum.fasta --pattern AAA --algorithm brute-force --output count`

- A filename (the .fasta file) (--file; required, no default value).
- A list of patterns to search for (--pattern; required at least once, no default value).
- The algorithm to use (--algorithm; optional, default value is brute-force, other values are bmh, bmh-pattern-alphabet, and aho-corasick).
- What to print: count of matches, locations, or both (--output; optional, default value is count, other values are locations or both).

_Expected output_:

```
BRUTE FORCE:
Pattern: AAA
Count: 24811
```

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
  - We know that the worst case of Brute-force is O(n x m),from the graph, we see that the time increases expectedly as we increase the text, therefore, it matches the expected worst-case performance.
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

![Experiment 1](/Experiment2.png)

### Experiment 2 Q&A

`How does the performance of your algorithm implementations change as you increase the number of patterns in the set? You should evaluate for at least the following number of patterns: 1, 5, 10.
Does the performance of your implementations match the expected worst-case performance for the algorithm as you increase the number of patterns? Why or why not?`

- As we can see from the result above, as we increase the number of patterns, the time taken to finish the program also increases. Though, we didn't see a lot of time increase in the Aho-corasick algorithm, since it's an algorithm that that has the main purpose of detecting multiple different patterns.
- Yes, the performance of the implementations match the expected worst case performance. We can see that Brute-force has an exponential growth. Both BMH and BMH-pattern-alphabet also has a decent growth, though it is not as extreme as brute-force since the worst case run time is O(m x n). Lastly, we can see how fast the Aho-corasick algorithm is compared to the others, since aho-corasick is using a trie that can handle multiple different patterns easily, where it has the worst case of O(m + n).

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

#### Aho-Corasick

```

hyperfine "python3 main.py --pattern ATGA --file Sorangium_cellulosum.fasta --algorithm aho-corasick" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --file Sorangium_cellulosum.fasta --algorithm aho-corasick" "python3 main.py --pattern AGCT --pattern AGCA --pattern GACG --pattern GAAC --pattern ATGC --pattern ATAT --pattern AAAA --pattern GGGG --pattern CCCC --pattern TTTT --file Sorangium_cellulosum.fasta --algorithm aho-corasick"

```

## Experiment 3

### Result

![Experiment 1](/Experiment3.png)

### Experiment 3 Q&A

`How does the performance of your algorithm implementations change as you increase the length of the pattern(s) that you are using? You should evaluate for at least the following lengths of patterns: 1, 3, 7.
Does the performance of your implementations match the expected worst-case performance for the algorithm as you increase the length of patterns? Why or why not?`

- For brute-force and Aho-corasick, we can see that there's not much of a difference when we increase the length of the patterns since the increase in the length is not drastic. On the other hand, with BMH and BMH pattern alphabet, we can see the drastic decrease in time as we increase the pattern.
- For brute-force and Aho-corasick, yes, the performance does match the expected worse case performance of the algorithm. On the other hand, for both BMHs, the performance also match the expected worst case but with length = 1 something needs to be analyzed. Why? So, for BMH, as we only use pattern with length = 1, this means that we could only skip or shift the pattern by 1 (this is similar to Brute-force). But why does brute force has a much better performance with length 1? This is the part where BMH has an issue, since in BMH, we would need to keep track of the index and also retrieve the shift value from the table for each iteration, this is the reason why BMH is so much worse than Brute-force when length = 1. Furthermore, as we increase the length of the pattern, now we can see the real potential of BMH, where it would shift and skip the pattern further than what bruteforce can do. Thus, makes the performance better.

#### Brute-force

```

hyperfine "python3 main.py --pattern A --file Sorangium_cellulosum.fasta --algorithm brutefoce" "python3 main.py --pattern ATG --file Sorangium_cellulosum.fasta --algorithm brutefoce" "python3 main.py --pattern ATGATAA --file Sorangium_cellulosum.fasta --algorithm brutefoce"

```

#### BMH

```
hyperfine "python3 main.py --pattern A --file Sorangium_cellulosum.fasta --algorithm bmh" "python3 main.py --pattern ATG --file Sorangium_cellulosum.fasta --algorithm bmh" "python3 main.py --pattern ATGATAA --file Sorangium_cellulosum.fasta --algorithm bmh"

```

#### BMH-Pattern-Alphabet

```

hyperfine "python3 main.py --pattern A --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern ATG --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet" "python3 main.py --pattern ATGATAA --file Sorangium_cellulosum.fasta --algorithm bmh-pattern-alphabet"

```

#### Aho-Corasick

```

hyperfine "python3 main.py --pattern A --file Sorangium_cellulosum.fasta --algorithm aho-corasick" "python3 main.py --pattern ATG --file Sorangium_cellulosum.fasta --algorithm aho-corasick" "python3 main.py --pattern ATGATAA --file Sorangium_cellulosum.fasta --algorithm aho-corasick"

```

## Experiment 4

### Result

![Experiment 1](/Experiment2.png)

### Experiment 4 Q&A

`How does the performance of your implementations of Boyer-Moore-Horspool compare?
Is either implementation observably faster or slower? If one is faster or slower, why is this happening? If they are both the same overall speed, why do you think this is happening?`

- The performance of BMH-pattern-alphabet is faster in all cases. For instance, the BMH-pattern-alphabet appears to be faster when we did an experiment where we increase the number of patterns. The reason why pattern-alphabet is faster is because when we did the pattern processing, this means that we eliminate any potential characters from the text that are not going to appear in the pattern. This results in a reduced time for searching that specific character in the table. Thus, BMH-pattern-alphabet is faster than the bmh with a text processing.
