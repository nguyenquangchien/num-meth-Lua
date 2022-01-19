 --/* File list4_12.lua */ -- Equations for chemical reaction
 
require"nsolv"

P = 20 -- Single parameter
eqs = function(f,x)
	f[1] = x[1]/2 + x[2] + x[3]/2 -x[6]/x[7]
	f[2] = x[3] + x[4] + 2*x[5] - 2/x[7]
	f[3] = x[1] + x[2] + x[5] -1/x[7]
	f[4] = x[1] + x[2] + x[3] + x[4] + x[5] - 1
	f[5] = P^2*x[1]*x[4]^3 - 1.7837e5*x[3]*x[5]
	f[6] = x[1]*x[3] - 2.6058*x[2]*x[4]
	f[7] = -28837*x[1] - 139009*x[2] - 78213*x[3] + 18927*x[4] +
		8427*x[5] + 13492/x[7] - 10690*x[6]/x[7]
end

x = {.5, 0, 0, .5, 0, .5, 2} -- Initial guesses
step = {0,0,0,0,0,0,0} -- Try different values
nmx,a = nsolv(eqs,x,step) -- Solve equations
print('Number of iterations, error = ',nmx,a)
table.foreach(x,print)
