-- /* File list8_9.lua */

require"prob"
local Pigamma = elemfunc.Pigamma
x,y = {},{} -- Display Incomplete Gamma function
j,gam = 1,10
for i=0,22,.1 do x[j],y[j],j = i,Pigamma(i,gam),j+1 end
plot(x,y); write_data('list8_9.dat',x,y) -- plot &save data e

iCDF = function(CDF,pbv,xinit,...) -- General inverse CDF function
	local arg = {...}
	local iCDFx =  function(x) -- Function for Newton's method
		return CDF(x,unpack(arg)) - pbv
	end
	return newton(iCDFx,xinit) -- Calculate inverse value
end

xinit = gam-1 -- Peak slope ocurs here
print(iCDF(Pigamma,.05,xinit,gam))
print(iCDF(Pigamma,.95,xinit,gam))
		