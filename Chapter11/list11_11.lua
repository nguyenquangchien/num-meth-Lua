-- /* File list11_11.lua */
-- Shooting method for boundary value plus eigenvalue problems

require"odebv"

f = function(eqs,E,x,u,up,upp) -- Quantum Harmonic oscillator
	eqs[1] = upp[1] + (E[1] - x^2)*u[1]
end

nx,xmax,Ei = 2000, 5.0, 5.1 -- Try different ranges
E = {Ei}
-- Try different initial conditions for shooting method
s,err,Ex = odebveve({f,E},{-xmax,xmax,nx},{0},{1})
--s,err,Ex = odebveve({f,E},{-xmax,xmax,nx},{1},{0})

plot(s,err)
sn,sqn,fac = sqnorm(s,1.0)
if fac>1e-5 then 
	print('Spatial range is not large enough for high accuracy')
end
nz = nzeros(sn)
print(xmax,nx,E[1],Ex[1],E[1]-2*nz-1,Ex[1]-2*nz-1)

plot(sn); plot(sqn) 
print('number of zeros = ',nz,fac)
write_data('list11_11.dat',s,sn,sqn)


