-- /* File list8_10.lua */

require"prob"

x,y = {},{}
read_data('steel_yield.txt',y)	-- Change to read in any desired file
nd = #y; ds = 'kpsi'
-- Do 4-Plots
for i=1,nd do x[i] = i end -- Generate run number
plot(x,y) -- Run plot
xx,yy = makeCDF(y) -- Generate distribution data
xn,yn = normalCDF(y)
plot(xx,yy,yn); write_data('list8_10.dat',xx,yy,yn) -- Distribution plot
scatterplot(lag(y)) -- Lag plot
xh,yh = hist(y,15); xph,yph = histnorm(y,15) 
plot(xh,yh,yph) -- Histogram plot
-- Do mean analysis with bounds
m,std = stats(y)
print('number of data points = ',nd); print('mean, std = ',engr(m,ds),engr(std,ds))
stderr = std/math.sqrt(nd-1)
probv = {.6827,.9545,.9,.95,.99,.999}
for i=1,#probv do
	ts = ttable(nd-1,probv[i]); dm = ts*stderr
	print(engr(m-dm,ds)..' < mean <'..engr(m+dm,ds)..' at'..engr(100*probv[i])..'percent confidence level')
end
print(engr(m-stderr,ds)..' < mean <'..engr(m+stderr,ds)..' for standard error bounds ')
print(engr(m-2*stderr,ds)..' < mean <'..engr(m+2*stderr,ds)..' for 2 standard error bounds ')
-- Do std analysis with bounds
stdf = std*math.sqrt(nd-1)
for i=1,#probv do
	tchL,tchH = tchisq(nd-1,probv[i])
	varL,varH = stdf/math.sqrt(tchL),stdf/math.sqrt(tchH)
	print(engr(varH,ds)..' <   std  <'..engr(varL,ds)..' at'..engr(100*probv[i])..'percent confidence level')
end
fdf = 1/math.sqrt(2*nd)
print(engr(std*(1-fdf),ds)..' <   std  <'..engr(std*(1+fdf),ds)..' for standard error bounds ')
print(engr(std*(1-2*fdf),ds)..' <   std  <'..engr(std*(1+2*fdf),ds)..' for 2 standard error bounds ')
