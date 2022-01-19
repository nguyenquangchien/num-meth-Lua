-- File list4_6.lua --nsolv() with linear equations
require "nsolv"

lineqs = function(f,x)
	f[1] = 10*x[1] - 7*x[2] - 7
	f[2] = -3*x[1] + 2*x[2] + 6*x[3] - 4
	f[3] = 5*x[1] - x[2] + 5*x[3] - 6
end

x = {0,0,0} -- Initial guess
getfenv(nsolv).linear=1 -- One iteration needed for linear equations
nmx,a = nsolv(lineqs,x)
print('For lineqs : ',x[1],x[2],x[3])
print('Number of iterations, error = ',nmx,a)
