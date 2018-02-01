# Dracula
DRACULA: Decay, Recombination And Collisions compUtation for Low-temperature Antihydrogen

Results published in Phys. Rev. A 90, 032704 (2014), B. Radics, D.J. Murtagh, Y. Yamazaki and F. Robicheaux
DOI:https://doi.org/10.1103/PhysRevA.90.032704

Requirements: sqlite3 (tested 3.20.1) and Gnu Scientific Library (tested 1.16)

The code reads in rate coefficients from an external sqlite3 data base file (~60 MB),
which is available for download at: https://drive.google.com/open?id=17kysiRw3CseA4OU4iRft3aLnoH4rMVae

1. To build the code basically use "compile.sh"

2. Running the code is possible in 2 ways: command line options or steering from the code.
- For the command line option the "main" executable can be used e.g. as follows:

```
$ ./main -f data/hydrogenDB.db -t 100 -d 1e+13 -an 1e+06 -b 1.0 -ot 5 -r 300 -otp ascii -of blah.txt

SQLite filename : data/hydrogenDB.db
positron temp : 100 K
positron dens : 1e+13 /m3
number of antiprotons : 1e+06
magnetic field : 1 T
overlap time : 5 us
black-body radiation temp : 300 K
Reading in scat.rates from sqlite3 file data/hydrogenDB.db...
Total tbr rate: 2.09077 /s
Reading in rad. recombination rates  from sqlite3 file data/hydrogenDB.db...
Adding stimulated absorption/emission due to Black-Body radiation...
Reading in decay rates  from sqlite3 file data/hydrogenDB.db...
Starting integration...
        1.000000E-09	         1.843444E-08
        6.000000E-09             1.419992E-07
        2.964199E-08             8.732732E-07
        1.071688E-07             3.685034E-06
        3.281334E-07             1.274194E-05
        9.227646E-07             3.983970E-05
        2.404090E-06             1.148360E-04
        5.000000E-06             2.627611E-04
Integration time: 1.0 seconds
Total number of bound states:         4.280536E+00, burned antiprotons:         4.280536E+00
Saving result to ascii file: blah.txt
```

The output file "blah.txt" contains two columns: n principal q.n. and level population.

- An example for steering the computation from a code is added "mainFlight.cxx"
```
./mainFlight data/hydrogenDB.db blah.txt
flight time: 1000 us
Reading in decay rates  from sqlite3 file data/hydrogenDB.db...
Starting integration...
        1.000000E-09	         4.225419E-08
        6.000000E-09             2.221629E-07
        3.100000E-08             9.637613E-07
        1.441408E-07             3.644229E-06
        6.358317E-07             1.252652E-05
        2.583181E-06             3.855787E-05
        9.662858E-06             1.116572E-04
        3.440850E-05             3.325520E-04
        1.167149E-04             1.016909E-03
        3.892226E-04             3.085021E-03
        1.000000E-03             6.996756E-03
Integration time: 0.0 seconds
Total number of bound states:         4.107787E+00, burned antiprotons:         0.000000E+00
Saving result to ascii file: blah.txt_flight.dat
```
