-- /* File list11_7.lua */
-- Shooting method for boundary value problem with nonlinear DE
require"odebv" -- Boundary value by shooting method solver
odeiv = odeivs -- Use adaptive step IV solver

nx,ubL,ubR = 2000, 0, 2 -- #points, Left, Right Boundary values

f = function(eqs,x,u,up,upp) -- Differntial equation
	eqs[1] = upp[1] + 4*up[1]^2
end	
-- Define Left and Right boundary equations
fb = function(bv,uL,uR) 
	bv[1] = uL[1] - ubL -- Left boundary
	bv[2] = uR[1] - ubR -- Right boundary
end

s1,n1,n2,err = odebvst({f,fb},{0,1},{0},{0}) -- Solve BV problem
print(n1,n2,err)
plot(s1[1],s1[2]) 
write_data("list11_7.dat",s1)
