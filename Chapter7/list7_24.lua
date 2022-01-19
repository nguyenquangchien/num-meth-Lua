-- /* File list7_24.lua */ -- Two step fitting to rational polynomial

require"nlstsq"

fft = function(yx,c) -- Final fitting function, Nonlinear in C's
	x = yx[2]
	return (c[1]+ c[2]*x+c[3]*x^2+c[4]*x^3)/(1+c[5]*x+c[6]*x^2+c[7]*x^3) - yx[1]
end
fft1 = function(yx,c) -- First approximation, Linear in C's
	x,xy = yx[2],yx[2]*yx[3]
	return (c[1]+ c[2]*x+c[3]*x^2+c[4]*x^3) - xy*(c[5]+c[6]*x+c[7]*x^2) -yx[1]
end

xd,yd={},{}
nd = read_data('thurber.dat',yd,xd) -- Read data
c = {0,0,0,0,0,0,0}; nc = #c
actv,step = {},{}
for i=1,nc do actv[i],step[i] = 1,0 end -- No limits on steps
yx = {yd,xd,yd} -- Data values, pass y values as third argument

del,err,nmax = nlstsq(yx,fw,fft1,c,actv,step) -- Step 1 -- Linear in C's
print(err,nmax,'-- First Approximation') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

for i=5,nc do step[i] = 1.2 end -- Set limits on steps for denominator coefficients
del,err,nmax = nlstsq(yx,fw,fft,c,actv,step) -- Step 2 -- Nonlinear in C's
print(err,nmax,'-- Final Calculation') -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

xmin,xmax,nx = xd[1],xd[nd],1000
dx = (xmax-xmin)/nx; xx,yy={},{}
for i=1,nx do 
	xx[i] = xmin+dx*(i-1)
	yy[i] = fft({0,xx[i]},c)
end
plot(xx,yy)
write_data("list7_24.dat",xx,yy) -- Save data 


