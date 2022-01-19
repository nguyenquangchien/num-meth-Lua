-- /* File list11_2.lua */
-- Shooting method for boundary value problem with mixed conditions
require"odebv" -- Boundary value by shooting method solver

nx,xmax = 2000,1 -- Independent variable parameters
a1,b1,c1 = 1,1,0; a2,b2,c2 = -1,1,-3 -- Boundary condition parameters

f = function(eqs,x,u,up) -- Differntial equation
	eqs[1] = up[1] - u[2]
	eqs[2] = up[2] - (2*x/(x^2+1))*u[2] + (2/(x^2+1))*u[1] - x^2 - 1
end

-- Define Left and Right boundary equations
fbound = function(bv,uL,uR) -- uL,uR, Left, Right values 
	bv[1] = a1*uL[1] + b1*uL[2] + c1 -- Left, x = 0 condition
	bv[2] = a2*uR[1] + b2*uR[2] + c2  -- Right, x = L condition
end

s,nm,ns,err = odebvst({f,fbound},{0,xmax,nx},{0,0})
print(nm,ns,err); plot(s[1],s[2])
write_data("list11_2.dat",s)

