-- /* File list8.3.lua */

require"prob"; require"deriv" -- rderiv used to handle derivative at x = 0

x,y,yd = {},{},{}
v = {1,1.5,2,3,4,5,6}; nv = #v -- Define v values for Weibull function
-- Probability density calculated from derivative of distribution function
f = function(x) return Pweibull(x,vv) end -- Function for derivative

for i=1,nv do y[i],yd[i] = {},{} end
i = 1
for xv=0,5,.01 do -- Step over x values
		x[i] = xv
		for j=1,nv do -- Step over v values
			vv = v[j]
			y[j][i] = Pweibull(xv,vv) -- Accumulate distribution values
			yd[j][i] = rderiv(f,xv) -- Accumulate derivative values
		end
	i = i+1
end
plot(x,unpack(y))
plot(x,unpack(yd))
write_data('list8_4.dat',x,yd)


