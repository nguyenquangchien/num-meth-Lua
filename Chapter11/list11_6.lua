-- /* File list11_6.lua */
-- Shooting method for boundary value problem with mixed conditions
require"odebv" -- Boundary value by shooting method solver

E,I,w,nx,xmax = 1.e7,500,100,1000,100 -- Independent variable parameters
EI = E*I

f2 = function(eqs,x,u,up,upp) -- Two second-order Differntial equations
	eqs[1] = upp[1] - u[2] 
	eqs[2] = upp[2] - w/EI
end	
fb2 = function(bv,uL,uR,upL,upR) -- Left, Right values and derivatives
	bv[1] = uL[1]; bv[2] = upL[1] -- Second-order, u[1](Left) = u'[1](Left) = 0
	bv[3] = uR[1]; bv[4] = upR[1] -- u[1](Right) = u'[1](Right) = 0
end

-- Solution with estimated error
s2,err = odebvste({f2,fb2},{0,xmax,nx},{0,0},{0,0})
write_data(20,"list11_6.dat",s2,err) ;plot(s2);plot(err)


