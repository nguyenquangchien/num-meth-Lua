--/* File list4_3.lua */

require "gauss"; require"spgauss"

-- Example of use
A = {
	{10, -7, 0},
	{-3, 2, 6},
	{5, -1, 5}
}
B = {7, 4, 6}
gauss(A,B) -- or gauss(A,B,3)
table.foreach(B,print)
-- Example of use with spgauss()
A = {{10, -7, 0}, {-3, 2, 6}, {5, -1, 5}} -- Compact form
B = {7, 4, 6}
sol = spgauss(A,B) -- Optional form of call
table.foreach(sol,print)
