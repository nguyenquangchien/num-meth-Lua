-- /* File list7_12.lua */ -- Reversed independent and dependent variables

require"nlstsq"
xd,yd,ycalc = {},{},{}
read_data('list7_1.dat',xd,yd)

ft = function(x,c) -- Define function to fit data
	return c[1]*(1 - c[2]*math.exp(c[3]*x[2])) - x[1]
end
frev = function(x,c) -- Reverse variables to fitted function
	return ft({x[2],x[1]},c)
end

c = {1,1,-.2} -- Initial guess at coefficients
del,err,nmax = nlstsq({xd,yd},fw,frev,c) -- Call fiting 
print(err,nmax) -- print results
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
for i=1,#xd do ycalc[i] = ft({0,xd[i]},c) end
write_data('list7_12.dat',xd,yd,ycalc)
plot(xd,yd,ycalc)