# Example of p-n junction with circular boundary #
15 # Boundary point section #
 0: 0    0    .05  1
 1: .5   0    .01  1
 2: .75  0    .02  2
 3: 1.0  0    .02  2
 4: 1.0  .15  .02  2
 5: 1.0  .25  .01  1
 6: 1.0  1.0  .05  3
 7: 0    1.0  .05  3
 8: .5086 .0773 .01 4 # Circular arc #
 9: .5334 .1469 .01 4
10: .5721 .2023 .01 4
11: .6209 .2378 .01 4
12: .6750 .2500 .01 4
 # define interior points for two materials #
13: .5   0.5   0   1  # Material 1 marker #
14: .75  0.1   0   2  # Material 2 marker #
14 # Boundary sides section #
 0: 0  1  1
 1: 1  2  1
 2: 2  3  2
 3: 3  4  2
 4: 4  5  1
 5: 5  6  1
 6: 6  7  3
 7: 7  0  1
 8: 1  8  400 # p-n junction boundary #
 9: 8  9  400
10: 9 10  400
11: 10 11 400
12: 11 12 400
13: 12 5  400
