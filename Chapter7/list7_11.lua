-- /* File list7_11.lua */ -- Data fitting for Figure 7.1

require"nlstsq"

xd,yd = {},{}
read_data('list7_1.dat',xd,yd)

ft = function(x,c) -- Define function to fit data
	return c[1]*(1 - c[2]*math.exp(c[3]*x[2])) - x[1]
end

c = {1,1,-.2} -- Initial guess at coefficients
actv = {1,0,1}
del,err,nmax = nlstsq({yd,xd},fw,ft,c,actv) -- Call fiting 
print(err,nmax) -- print results
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
ycalc = {}
for i=1,#xd do 
	x = {0,xd[i]}
	ycalc[i] = ft(x,c)
end
write_data('list7_11.dat',xd,yd,ycalc)
plot(xd,yd,ycalc)