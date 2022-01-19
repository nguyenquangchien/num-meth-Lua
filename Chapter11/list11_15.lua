-- /* File list11_15.lua */
-- Use of boundary value solver for second order equation

require"odebvfd"; getfenv(odebv1fd).nprint=1

-- Parameters
E,I,w,L = 1.e7, 500, 100, 100; EI = E*I

f = function(x,y,yp,ypp) -- Differntial equation
	return  ypp - w*x*(x-L)/(2*EI)
end
fl = function(y,yp) return y end -- Left boundary, y=0 
fr = function(y,yp) return y end -- Right boundary, y=0

x,y,y1,nx = {},{},{}, 501; dx = L/(nx-1)
for i=1,nx do
	x[i], y[i] = (i-1)*dx, 0
	y1[i] = w/(24*EI)*(x[i]^4-2*L*x[i]^3+L^3*x[i])
end
print(odebv1fd({f,fl,fr},x,y))
plot(x,y,y1); write_data('list11_15.dat',x,y,y1)
