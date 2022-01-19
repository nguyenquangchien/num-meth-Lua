-- /* File list11_12.lua */
-- Shooting method for fourth power potential well

require"odebv"

f = function(eqs,E,x,u,up,upp) -- Fourth order potential
	eqs[1] = upp[1] + (E[1] - x^4/10)*u[1]
end

nx,xmax,Ei = 2400, 5.5,30 -- Typical parameters
E = {Ei}
s,err,Ex = odebveve({f,E},{-xmax,xmax,nx},{0},{1})
sn,sp,fac = sqnorm(s,1.0)
if fac>1.e-5 then -- Test for spatial interval increase
	print('Need to increase spatial interval')
end

print(nx,xmax,Ex[1],E[1],nzeros(s))
plot(sn); plot(sp)
write_data('list11_12.dat',s,sn,sp,err)
