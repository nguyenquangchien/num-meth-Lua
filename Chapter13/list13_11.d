# Example of capacitor with two dielectric materials #
8  # Boundary points section #

 0:   0    0    .04   1 # Boundary markers 1 to 4 #
 1:   2.0  0    .04   1
 2:   2.0  0.5  .04   2
 3:   2.0  1.0  .04   3
 4:   0.0  1.0  .04   3
 5:   0.0  0.5  .04   4 

# Define internal points for two materials #
 6:   1.0  0.25 0     1 # Material 1 marker point #
 7:   1.0  0.75 0     2 # Material 2 marker point #
 
7 

 0:   0   1   1 # Boundary markers 1 to 4, consistent with points #
 1:   1   2   2
 2:   2   3   2
 3:   3   4   3
 4:   4   5   2
 5:   5   0   2

# Define internal boundary line between two materials #
 6:   2   5   400 # Note, use marker number > 100 for internal boundaries!!! #   
