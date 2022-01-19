require "nsolv"; require"deriv"
exp = math.exp

x0,x1 = 8,20
x3 = x0/1.25;x4 = x0/10
local cx
eqs = function(f,c)
	x10,x90,xm = c[4],c[4]+x3,c[5]
	x50 = c[4] -x4 + x1
	f[1] = fx(x10,c) -.1 -- 10% point, c[4] = time
	f[2] = fx(x90,c) -.9 -- 90% point c[4]+x3 = time
	f[3] = fx(x50,c) - .5 -- 50% point x1 = time
	f[4] = deriv(fx,xm,c) -- Derivative at peak = 0
	f[5] = fx(xm,c) -1 -- Peak value = 1 at c[5] = time
end
fx = function(x,c)
	return c[3]*x*exp(-x/c[2]-(x/c[1])^2)
end

c = {x1,3*x1,1,x0/10,1.5*x0}
getfenv(nsolv).nprint=1
nmx,a = nsolv(eqs,c,step)
print('at ',x0,x1)
print('c[1],c[2],c[3] = : ',c[1],c[2],c[3],c[4],c[5])
print('Number of iterations, error = ',nmx,a)
i = 1; xv,yv = {},{}
for x = 0,40,.1 do
	xv[i],yv[i] = x, fx(x,c)
	i = i+1
end
plot(xv,yv)
--write_data('tmp.dat',xv,yv)
print('10% point at ',c[4],fx(c[4],c))
print('90% point at',x3+c[4],fx(x3+c[4],c))
print('peak point at',c[5],fx(c[5],c))
print('50% point at',c[4]+x1-x4,fx(c[4]+x1-x4,c))
print('10% to 90% times 1.25 =',x3*1.25)
print(fx(20,c))