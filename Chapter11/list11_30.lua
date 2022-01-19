-- /* list11_30.lua */
-- Solution of nonlinear third order BV problem by two methods
require "odebv"; require"odefd"
getfenv(odefd).nprint=1
getfenv(odefd).umin = {1.e-3,1.e-3}

L = 10
f = function(eqs,x,u,up,upp) -- Equations for FD approach
	eqs[1] = upp[1] - up[2]
	eqs[2] = 2*upp[2] + u[1]*up[2]
end
fbl = function(eql,u,up) -- Left boundary conditions for FD
	eql[1] = u[1] 
	eql[2] = u[2]
end
fbr = function(eqr,u,up) -- Right boundary conditions for FD
	eqr[1] = up[1] - u[2]
	eqr[2] = u[2] - 1
end

fst = function(eqs,x,u,up) -- Equations for shooting approach
	eqs[1] = up[1] - u[2]
	eqs[2] = up[2] - u[3]
	eqs[3] = 2*up[3] + u[1]*u[3]
end
fb2 = function(bv,uL,uR) -- Boundary values for shooting method
	bv[1] = uL[1]; bv[2] = uL[2] -- Left boundary values
	bv[3] = uR[2] -1 -- Right boundary values
end	

x = {0,L,2000}; u = {{0,10},{0,1}} -- FD method
s1,n1 = odefd({f,fbl,fbr},x,u); plot(s1)
x = {0,L,2000}; u = {0,0,1} -- Shooting method
s2,n2 = odebvst({fst,fb2},x,u); plot(s2); print(n1,n2)
write_data('list11_30.dat',s1,s2)