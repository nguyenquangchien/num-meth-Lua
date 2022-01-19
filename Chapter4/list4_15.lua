 --/* File list4_15.lua */ -- Example of large number of equations
 
require"nsolv"

N = 1000 -- Number of equations
eqs = function(f,x)
	f[1] = 3*x[1] - 2*x[1]^2 -2*x[2] + 1
	for i=2,N-1 do f[i] = 3*x[i]-2*x[i]^2-x[i-1]-2*x[i+1]+1 end
	f[N] = 3*x[N] - 2*x[N]^2 - x[N-1] + 1
end

x = {}; for i=1,N do x[i] = -1 end -- try -0.1 to -10
nmx,a = nsolv(eqs,x)
print(nmx,a); table.foreach(x,print)
