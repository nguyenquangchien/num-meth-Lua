-- /* list11_24.lua */
-- Solution of nonlinear Chemical Engineering BV problem
require"odefd"
getfenv(odefd).nprint = 1

M,EdR,b,Q,ce,Te = 2, 22000, .5e8, 1000, 0.07, 1250

f = function(eqs,x,u,up,upp) -- Differntial equation
	eqs[1] = upp[1]/M - up[1] - b*u[1]^2*math.exp(-EdR/u[2])
	eqs[2] = upp[2]/M - up[2] +Q*b*u[1]^2*math.exp(-EdR/u[2])
end
fl = function(eqs,u,up) -- Left boundary conditions
	eqs[1] = u[1] - up[1]/M - ce
	eqs[2] = u[2] -up[2]/M - Te
end
fr = function(eqs,u,up)
	eqs[1] = up[1] -- Zero slopes
	eqs[2] = up[2]
end

u = {{0,0},{1250,1250}} -- Constant initial guesses, end points only
x = {0,48,1000} -- Linear grid spacing
s,se,nn = odefde({f,fl,fr},x,u)
plot(s[1],s[2]); plot(s[1],s[3])
write_data('list11_24.dat',s,se)

