--/* File list4_22.lua */

require "Matrix" -- Load matrix manipulation code

--Tests of matrix operations	
ax = Matrix.new({{1,1,1},{2,2,1},{3,3,1}})
ay = Matrix.new({{3,1,0},{1,-2,1},{-1,2,-1}})
az = ax+ay -- Matrix addition
print("Type of az is: ",type(az))

az = ax*ay -- Matrix multiplication
print('Matrix az is:'); print(az)
a = Matrix.new{ -- New matrix definition
	{1, 0.1, 0, 0, 0},
	{0.1, 1, 0.1, 0, 0},
	{0, 0.1, 1, 0.1, 0},
	{0 ,0, 0.1, 1, 0.1},
	{0, 0, 0, 0.1, 1}}
m,n = Matrix.size(a) -- Obtain size of matrix a
print("Size of matrix a is: ",m,n)
-- Take inverse of matrix
b = a^-1 -- Or b = Matrix.inv(a) or b = a:inv()
print(b,'\n')
print(b*a) -- Should be unity diagonal matrix