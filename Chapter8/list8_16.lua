-- /* File list8_16.lua */ -- Weibull distribution analysis, steel data

require"prob"; require"nlstsq"

life = {}
read_data('steel_yield.txt',life); print(stats(life))
life,F = makeCDF(life); nd = #life -- Make distribution function data

Fwloglog = function(x,c) -- log-log Weibull function
	local xx = x[2] - c[3]
	if xx<0 then return -55 -x[1] end
	return c[1]*math.log(xx) - c[2] - x[1]
end

xx,yy = {},{}
for i=1,nd do -- Convert F values to log-log Weibull data
	xx[i] = life[i]
	yy[i] = math.log(-math.log(1-F[i]))
end

c = {4,xx[nd]/2,xx[1]/1.1} -- Now non-zero c[3]
del,err,nn = nlstsq({yy,life}, fw, Fwloglog, c); print('Number of Newton steps =',nn)
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
print('m,s,mu = ',c[1],math.exp(c[2]/c[1]),c[3])

xp,yp = {},{}
xs = {0,0}
for i=1,nd do -- Theoretical Weibull function 
	xp[i] = math.log(life[i])
	yp[i] = Fwloglog({0,life[i]},c)
end
scatterplot(xp,yp,yy)
write_data('list.8_16.dat',life,yy,yp)