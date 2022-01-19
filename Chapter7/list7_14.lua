-- /* File list7_14.lua */ -- data fitting for Figure 7.1 with nlstsq()

require"nlstsq"

xd,yd = {},{}
read_data('list7_1.dat',xd,yd)

ft = function(x,c) -- Define function to fit data
	local xx = x[2]
	return c[1]*math.sin(c[2]*xx) +c[3]*math.sin(c[4]*xx)- x[1]
end

c = {1,2*math.pi/80,.1,6*math.pi/80} -- Initial guess at coefficients
del,err,nmax = nlstsq({yd,xd},fw,ft,c) -- Call fiting 
print('RMS error =',err,'#Iter =',nmax) -- print results
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
print('periods = ',2*math.pi/c[2],2*math.pi/c[4])
ycalc,y1,y2 = {},{},{}
for i=1,#xd do 
	y1[i] = c[1]*math.sin(c[2]*xd[i])
	y2[i] = c[3]*math.sin(c[4]*xd[i])
	ycalc[i] = y1[i] + y2[i]
end
write_data('list7_14.dat',xd,yd,ycalc,y1,y2)
plot(xd,yd,ycalc)