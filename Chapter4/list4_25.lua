 --/* File list4_25.lua */ -- Several LU methods for solving linear equations
 
require"Matrix"; require"gauss"
-- define some test matrices
aa = Matrix.new{{2,1,1,3,2},{1,2,2,1,1},{1,2,9,1,5},{3,1,1,7,1},{2,1,5,1,8}}
bb = Matrix.new{{-2},{4},{3},{-5},{1}} -- Single column matrix

a,b = aa:new(), bb:new() -- make copy of then
print('Determinant of a = ',a:det(),'\n') -- Modifies a

a = aa:new() -- Need new copy 
ainv = a^-1 -- One form of obtaining inverse
print(ainv*b,'\n') -- #1 solve equation set

-- a and b are unchanged by above operations
ainv = a:inv() -- Second form for inverse -- changes a matrix
print(ainv*b,'\n') -- #2 solve equation set

a = aa:new() -- New copy for a, b is still unchanged
fi = a:LUfunction() -- Function for solving with various b arrays
sol = fi(b) -- #3 Call function for solution
table.foreach(sol,print); print('\n') -- Solve and print results

a = aa:new() -- New copy for a, b is still unchanged
sol = Matrix.LUsolve(a,b) -- #4 Another method of solving equation set
table.foreach(sol,print); print('\n')

a = aa:new() -- New copy, b still unchanged
b = {-2, 4, 3, -5, 1} -- Single table needed for gauss
sol = gauss(a,b) -- #5 Gauss elimination, does not use LU 
table.foreach(sol,print)