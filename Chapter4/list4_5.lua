-- File list4._5.lua  --Examples of use of nsolv()
require "nsolv"

eqs1 = function(f,x)
	f[1] = x[1]^2 + 50*x[1] + x[2]^2 + x[3]^2 -200
	f[2] = x[1]^2 + 20*x[2] + x[3]^2 - 50
	f[3] = -x[1]^2 - x[2]^2 + 40*x[3] + 75
end

x = {0,0,0} -- Initial guess
nmx,a = nsolv(eqs1,x)
print('For eqs1 : ',x[1],x[2],x[3])
print('Number of iterations, error = ',nmx,a)

eqs2 = function(f,x)
	f[1] = x[1]*x[1] + math.exp(x[2]) + 1.0/x[3] - 41.
	f[2] = x[1] + 1./x[2] + x[3] - 6.
	f[3] = x[1]*x[1] - 24.0
end
x = {4,4,4}
nmx,a = nsolv(eqs2,x)
print('For eqs2 : ',x[1],x[2],x[3])
print('Number of iterations, error = ',nmx,a)

