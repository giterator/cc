# Concentric Circle (CC) Compression Algorithm

A new run length encoding algorithm for lossless data compression that exploits positional redundancy by representing data in a two-dimensional model of concentric circles is presented. This visual transform enables detection of runs (each of a different character) in which runs need not be contiguous and hence, is a generalization of run length encoding. Its compression factor is benchmarked against TurboRLE's by running it on the Silesia Compression Corpus. The results are shown in the bar chart below.

![Alt text](performance.png?raw=true "CC vs TurboRLE")



CC's overall compression factor: 1.10

TurboRLE's overall compression factor: 1.07



CC's compression factor is greater than TurboRLE's by 2.8%

