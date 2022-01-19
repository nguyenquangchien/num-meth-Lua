-- /* File list9_1.lua */ -- Data fitting and 6 plots for Figure 9.1

requires("nlstsq", "prob", "DataFit")

xd,yd = {},{} -- Data arrays
read_data('list7_1.dat',xd,yd); nd = #xd -- Read in data set to be fitted

ft = function(x,c) -- Define function to fit data
	return c[1]*(1 - math.exp(c[2]*x[2])) - x[1]
end

c = {1,-.2}; nc = #c -- Initial guess at coefficients
del,err,nmax = nlstsq({yd,xd},fw,ft,c) -- Call fiting 
print('err, nmax =',err,nmax) -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
-- Generate 6-Plot data
ycalc,res,yres1,yres2 = {},{},{},{}
for i=1,nd do -- Calculate fitted values and residuals
	ycalc[i] = ft({0,xd[i]},c)
	res[i] = yd[i]-ycalc[i]
end
fres1,fres2 = datafit({res,xd}), datafit({res,ycalc})  -- Fit to residual data
for i=1,nd do yres1[i],yres2[i] = fres1(xd[i]), fres2(ycalc[i]) end
xp,yp = makeCDF(res); xpn,ypn = normalCDF(res) -- CDF data for residuals
xh,yh = hist(res,15); xhn,yhn = histnorm(res,15) -- Histogram data
xl,yl = lag(res) -- lag data for residuals
write_data('list9_1.dat',xd,yd,ycalc,res,yres1,yres2,xp,yp,ypn,xl,yl,xh,yh,yhn)
-- Now show 6 plots
plot(xd,yd,ycalc) -- Plot 1 -- Predicted values and data points
plot(xd,res,yres1) -- Plot 2 -- Residuals vs x
plot(ycalc,res,yres2) -- Plot 3 -- Residuals vs y
scatterplot(xl,yl) -- Plot 4 -- Lag plot of residuals
plot(xh,yh,yhn) -- Plot 5 -- Histogram plot of residuals 
plot(xp,yp,ypn) -- Plot 6 -- Distribution plot of residuals

