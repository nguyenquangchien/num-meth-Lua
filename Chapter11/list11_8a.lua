-- /* File list11_8a.lua */
-- Shooting method for boundary value plus eigenvalue problems
require"odebv"
			
f = function(eqs,E,x,u,up,upp) -- Differntial equation
	eqs[1] = upp[1] + E[1]*u[1]
end

nx,xmin,xmax = 2000,0,1
Ei = 0; E = {Ei} -- Guess at an eigenvalue
-- Set initial value to zero and derivative to 1.0
s,ns,nm,err = odebvev({f,E},{xmin,xmax,nx},{0},{1})

plot(s); print(E[1],ns,nm,err)
print('Eigenvalue error =',E[1]-math.pi^2)
write_data('list11_8.dat',s)
print('number of zeros = ',nzeros(s))
E1 = E[1]
nx = 1000
s,ns,nm,err = odebvev({f,E},{xmin,xmax,nx},{0},{1})
Eext = (4*E1 - E[1])/3
print(E[1],Eext, Eext-math.pi^2)


