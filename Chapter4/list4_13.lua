--/* File 4_13.lua */ -- Parameters for time pulse
require "nsolv"; require"deriv"
exp = math.exp

--t1,t2 = 10,1000 -- Select pulse as desired
--t1,t2 = 10,350
t1,t2 = 5,50
--t1,t2 = 1.2,50
t3,t4 = t1/1.25,t1/10

eqs = function(f,c)
	t10,t90 = c[4],c[4]+t3,c[5]
	tm,t50 = c[5],c[4] -t4 + t2
	f[1] = fx(t10,c) -.1 -- 10% point, c[4] = time
	f[2] = fx(t90,c) -.9 -- 90% point x90 = time
	f[3] = fx(t50,c) - .5 -- 50% point x50 = time
	f[4] = deriv(fx,tm,c) -- Derivative at peak = 0
	f[5] = fx(tm,c) -1 -- Peak value = 1 at c[5] = time
end
fx = function(t,c)
	return c[3]*(1-exp(-t/c[1]))*exp(-t/c[2])
end

c = {t1/2,1.5*t2,1,t1/10,1.5*t1}; nc = #c
nn,err = nsolv(eqs,c,step)
print('at ',t1,t2)
print('Number of iterations, error = ',nn,err)
for i=1,nc do printf('c[%d] = %12.4e\n',i,c[i]) end
i = 1; xv,yv = {},{}
for t = 0,2*t2,t2/1000 do
	xv[i],yv[i] = t, fx(t,c); i = i+1
end
plot(xv,yv); write_data('list4_13.dat',xv,yv)
print('10% point at ',c[4],fx(c[4],c))
print('90% point at',t3+c[4],fx(t3+c[4],c))
print('peak point at',c[5],fx(c[5],c))
print('50% point at',c[4]+t2-t4,fx(c[4]+t2-t4,c))
