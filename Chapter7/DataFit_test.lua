-- /* File Data_fit.lua */
-- Program to fit data with noise to a smooth curve for plotting

require "DataFit"; require"deriv"

f = function(x) -- Function to generate some data
	return x/(x^3 + (10^3))^(1/3)
end

xd,yd={},{}
j = 1
for i=0,20,.4 do -- Generate data with noise
	xd[j] = i
	yd[j] = f(i) + (-.5+math.random())/20
	j = j+1
end
yd[1] = 0
nd = #xd -- Number of data points
fxc,xcc,err = datafit({yd,xd}) -- Now do fit with floating knots, 6 knots
-- fxc(x) is a returned function which is the fitting function with fixed knots set by nlstsq.
-- xcc = {xc,c,del} are the floating knot locations and values from least squares fitting
-- err is the mean square error per data point.
-- Can call with just fxc = datafit{yd,xd} or datafit(yx) if other parameters are not important

for i=1,#xcc[2] do printf('c[%d] = %12.4e +/- %12.4e\n',i,xcc[2][i],xcc[3][i]) end
print("err = ",err)

ycalc,yexact,ycalc1 = {},{},{}
for j=1,nd do
	ycalc[j] = fxc(xd[j])
	yexact[j] = f(xd[j])
end
plot(xd,yd,ycalc)
fxc,xcc,err = datafitn({yd,xd},2)
xc,c,del = unpack(xcc)
for i=1,#c do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

--for i=1,#xcc[2] do printf('c[%d] = %12.4e +/- %12.4e\n',i,xcc[2][i],xcc[3][i]) end
print("err = ",err)

for j=1,nd do
	ycalc1[j] = fxc(xd[j])
end
plot(xd,yd,ycalc)
write_data("test.dat",xd,yd,yexact,ycalc,ycalc1)


