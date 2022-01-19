# Example of p-n junction #
11 # Boundary point section #

 0: 0    0    .05  1
 1: .5   0    .01  1
 2: .75  0    .02  2
 3: 1.0  0    .02  2
 4: 1.0  .15  .02  2
 5: 1.0  .25  .01  1
 6: 1.0  1.0  .05  3
 7: 0    1.0  .05  3
 8: 0.5  0.25 .01  4

 # define interior points for two materials #
 9: .5   0.5   0   1  # Material 1 marker #
 10: .75  0.1  0   2  # Material 2 marker #

10 # Boundary sides section #

 0: 0  1  1
 1: 1  2  1
 2: 2  3  2
 3: 3  4  2
 4: 4  5  1
 5: 5  6  1
 6: 6  7  3
 7: 7  0  1

 8: 5  8  400 # p-n junction boundary #
 9: 8  1  400