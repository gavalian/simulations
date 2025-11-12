# simulations
Code for generating final state for given decays. 

# running the code 

Use ROOT to run the example codes:

```
root -l generate.C
```

# code output

The code outputs the simlated particles and their decays.
The header prints out the generated angles for production 
process and decay process.

first two numbers are the theta and phi angles in the rest 
frame of the reaction, for proton. where the ep decays to a 
proton and a vector meson.

the second two numbers are the theta and phi angles in of the 
decay products in the rest frame of the vector meson.

The header is followed by the printoutof all the vectors in the 
reaction:

```
header: production [  2.96388,   2.01860], decay [  1.23560   1.61794]
e  in:   0.00000   0.00000  10.60000 [  0.00050]
e out:   1.18395   0.00000   6.40598 [  0.00050]
p out:  -0.10248  -0.12368   0.58169 [  0.93827]
decay:  -1.08148   0.12368   3.61233 [  1.02000]
dgt 1:  -0.58369   0.18683   1.97271 [  0.49368]
dgt 2:  -0.49778  -0.06315   1.63962 [  0.49368]
```

# running for LUND output

run command:

```
root -l -q -b rhoproduce.C
```

or 

```
root -l -q -b phiproduce.C
```

the script prints out events in LUND format on the screen. the header contains
q2,xb, production theta and phi, and decay theta and phi, starting from column 11;
