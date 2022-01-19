-- /* File list11_4a.lua */
-- Shooting method for fourth order boundary value problem 
require"odebv" -- Boundary value by shooting method solver

E,I,w,nx,xmax = 1.e7,500,100,1000,100 -- Independent variable parameters
EI = E*I

f2 = function(eqs,x,u,up,upp) -- Two second-order Differntial equations
	eqs[1] = upp[1] - u[2] 
	if x>=40 and x<=60 then wx = w 
	else wx = 0.0 end
	eqs[2] = upp[2] - wx/EI
end	
fb2 = function(bv,uL,uR,upL,upR) -- Left, Right values and derivatives
	bv[1] = uL[1]; bv[2] = upL[1] -- Second-order, u[1](Left) = u'[1](Left) = 0
	bv[3] = uR[1]; bv[4] = upR[1] -- u[1](Right) = u'[1](Right) = 0
end

f4 = function(eqs,x,u,up) -- Four first-order Differential equations
	eqs[1] = up[1] - u[2] 
	eqs[2] = up[2] - u[3]
	eqs[3] = up[3] - u[4]
	if x>=40 and x<=60 then wx = w 
	else wx = 0.0 end
	eqs[4] = up[4] - wx/EI
end
fb4 = function(bv,uL,uR) -- Left, Right values, no derivatives
	bv[1] = uL[1];  bv[2] = uR[1]  -- First-order, u[1](Left) = u[2](Left) = 0
	bv[3] = uL[2]; bv[4] = uR[2] -- u[1](Right) = u[2](Right) = 0
end

s2 = odebvst({f2,fb2},{0,xmax,nx},{0,0},{0,0})
s4 = odebvst({f4,fb4},{0,xmax,nx},{0,0,0,0})

plot(s2[1],s2[2]); plot(s4[1],s4[2])

nx,dmax = #s2[1], 0
for i=1,nx do 
	x = math.abs(s2[2][i] - s4[2][i])
	if x>dmax then dmax = x end
end
print('maximum difference = ',dmax)
print('maximum deflection =',s2[2][math.floor(nx/2)])
write_data("list11_4.dat",s2,s4)


