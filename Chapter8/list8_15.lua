-- /* File list8_15.lua */ -- Weibull distribution analysis

require"prob"; require"nlstsq"

qbd,F = {},{}
read_data('qbd.dat.txt',qbd,F) -- Data in form of F already
nd = #qbd

Fwloglog = function(x,c) -- log-log Weibull function
	local xx = x[2] - c[3]
	if xx<0 then return -55 -x[1] end
	return c[1]*math.log(xx) - c[2] - x[1]
end

xx,yy = {},{}
for i=1,nd do -- Convert F values to log-log Weibull data
	xx[i] = qbd[i]
	yy[i] = math.log(-math.log(1-F[i]))
end

c = {4,xx[nd]/2,xx[1]/2} -- Initial guesses
del,err,nn = nlstsq({yy,qbd}, fw, Fwloglog, c); print('Number of Newton steps =',nn)
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
print('m,s,mu = ',c[1],math.exp(c[2]/c[1]),c[3])

xp,yp = {},{}
xs = {0,0}
for i=1,nd do -- Theoretical Weibull function 
	xp[i] = math.log(qbd[i])
	yp[i] = Fwloglog({0,qbd[i]},c)
end
scatterplot(xp,yp,yy)
write_data('list.8_15.dat',qbd,yy,yp)
