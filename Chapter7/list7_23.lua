-- /* File list7_23.lua */

require"nlstsq"

fft = function(yx,c)
	x = yx[2]
	return (c[1]+ c[2]*x+c[3]*x^2+c[4]*x^3)/(1+c[5]*x+c[6]*x^2+c[7]*x^3) - yx[1]
end

xd,yd={},{}
nd = read_data('thurber.dat',yd,xd) -- Read data
c = {0,0,0,0,0,0,0}
nc = #c
actv,step = {},{}
for i=1,nc do actv[i],step[i] = 1,0 end
actv[5],actv[6],actv[7] = 0,0,0
yx = {yd,xd} -- Data values

del,err,nmax = nlstsq(yx,fw,fft,c,actv,step) -- Call fiting - step 1
print(err,nmax,'-- Numerator only') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
c[5],actv[5],step[5] = .3,1,1.2 -- Include c[5]
del,err,nmax = nlstsq(yx,fw,fft,c,actv,step) -- Call fiting - step 2
print(err,nmax,'-- Add c[5]') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
c[6],actv[6],step[6] = .1,1,1.2 --Include c[6]
del,err,nmax = nlstsq(yx,fw,fft,c,actv,step) -- Call fiting - step 3
print(err,nmax,'-- Add c[6]') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
c[7],actv[7],step[7] = .03,1,1.2 -- Include c[7]
del,err,nmax = nlstsq(yx,fw,fft,c,actv,step) -- Call fiting - step 4
print(err,nmax,'-- Add c[7]') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
xmin,xmax,nx = xd[1],xd[nd],1000
dx = (xmax-xmin)/nx
xx,yy={},{}
for i=1,nx do 
	xx[i] = xmin+dx*(i-1)
	yy[i] = fft({0,xx[i]},c)
end
plot(xx,yy)
write_data("list7_23.dat",xx,yy) -- Save data 


