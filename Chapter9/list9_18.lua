-- /* File list9.18.lua */ -- Piecewise diode model 

requires("nlstsq","mcpar","DataFit","cltwodim")

infile = 'ivdata'
xd,yd = {},{}; read_data(infile..'.txt',xd,yd) -- x stored first
nd = #xd
for i=1,nd do yd[i] = math.log(yd[i]) end; plot(xd,yd)

n,m,Vt = 4,4,.0259 -- Experiment with values for best fit
ft = function(x,c) -- Define function to fit data
	v = (x[2] - c[4]*math.exp(x[1]))/Vt
	i12 = ((c[1]*(math.exp(v/2)-1))^n + (c[2]*math.exp(v))^n)^(1/n)
	return math.log(((i12)^-m + (c[3]*math.exp(v/2))^-m)^-(1/m)) - x[1]
end

c = {5e-9,5.6e-13,8e-9,4}; nc = #c -- Initial approximations

del,err,nn = nlstsq({yd,xd},fw,ft,c) -- Call fiting, print max iterations 
print('Newton iterations, err =',nn,err)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

ycalc,res = {},{} -- Calculate fitted values and residuals
for i=1,nd do ycalc[i] = newtonfc(ft,{yd[i],xd[i]},c); res[i] = yd[i] - ycalc[i] end
fres1,fres2 = datafit({res,xd}), datafit({res,yd}) -- Fit residuals
yres1,yres2,xyp,xypn = {},{},{},{}
for i=1,nd do yres1[i],yres2[i] = fres1(xd[i]),fres2((ycalc[i])) end
xyp[1],xyp[2] = makeCDF(res); xypn[1],xypn[2] = normalCDF(res) -- CDF data for residuals
plot(xd,yd,ycalc); plot(xd,yres2,res); plot(ycalc,yres2,res); plot(xyp,xypn) 
write_data(infile..'e.dat',xd,yd,ycalc,res,yres1,yres2,xyp,xypn)

