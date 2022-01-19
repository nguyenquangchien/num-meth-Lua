-- /* File list11_10.lua */
-- Shooting method for boundary value plus eigenvalue problems
require"odebv"
			
f = function(eqs,E,x,u,up,upp) -- Differntial equation
	eqs[1] = upp[1] + E[1]*u[1]
end

nx,xmin,xmax = 2000,0,1
Ei = 210; E = {Ei} -- Guess at an eigenvalue
-- Set initial value to zero and derivative to 1.0
s,err,Ex = odebveve({f,E},{xmin,xmax,nx},{0},{1})
plot(s,err) -- Solution and estimated error
s,sq = sqnorm(s,1.0); plot(sq) -- Normalize square integral to 1.0
nz = nzeros(s)+1; Et = (nz*math.pi)^2
write_data('list11_10.dat',s,err)
print('Eigenvalue Order = ',nz)
print('EV, Extrapolated EV, exact EV =',E[1],Ex[1],Et)
print('EV error, Extrapolated EV error =',E[1]-Et,Ex[1]-Et)

