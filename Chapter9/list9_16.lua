-- /* File list9_16.lua */ -- Parameter estimation for diode equation

requires("nlstsq","mcpar","DataFit")

-- Section 1. Input data file and define fitting function
-- Change lines 6 through 12 for different data sets
infile = 'ivdata'
xd,yd = {},{}; read_data(infile..'.txt',xd,yd) -- x stored first
nd = #xd; for i=1,nd do yd[i] = math.log(yd[i]) end

ft = function(x,c) -- Define function to fit data
	return  math.log(c[1]*(math.exp(x[2]/(0.0259*c[2])) - 1)) - x[1]
end

c = {1.e-9,2}; nc = #c; step ={0,1.2} -- Initial approximations. End of changes for most data sets
fw = {}; for i=1,nd do fw[i] = 1 end

-- Section 2. Perform data fit and print fitted parameters
del,err,nn = nlstsq({yd,xd},fw,ft,c,actv,step) -- Call fiting, print max iterations 
print("Number of iterations, err =",nn,err)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

-- Section 3. Generate 6-Plot Data, Only Plots 1,2,3 and 6 are used here
ycalc,res = {},{} -- Calculate fitted values and residuals
for i=1,nd do ycalc[i] = newtonfc(ft,{yd[i],xd[i]},c); res[i] = yd[i] - ycalc[i] end
fres1,fres2 = datafit({res,xd}), datafit({res,yd}) -- Fit residuals
yres1,yres2,xyp,xypn = {},{},{},{}
for i=1,nd do yres1[i],yres2[i] = fres1(xd[i]),fres2((ycalc[i])) end
xyp[1],xyp[2] = makeCDF(res); xypn[1],xypn[2] = normalCDF(res) -- CDF data for residuals
plot(xd,yd,ycalc); plot(xd,yres2,res); plot(ycalc,yres2,res); plot(xyp,xypn) 
write_data(infile..'c.dat',xd,yd,ycalc,res,yres1,yres2,xyp,xypn)

