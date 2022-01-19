--/* File list4_2.lua */

require "Matrix"

A = Matrix.new{ -- Coefficient matrix
	{10, -7, 0},
	{-3, 2, 6},
	{5, -1, 5}
}
B = Matrix.new{{7}, {4}, {6}} -- Column matrix
X = A^-1*B -- Solve -- or X = Matrix.inv(A)*B or A:inv()*B
print("Solution is : ")
print(X)
