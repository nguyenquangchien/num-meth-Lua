-- /* File list7_18.lua */  Simple example with specified knots

require "DataFit"

xd,yd,ycalc={},{},{}
read_data('list7_1.dat',xd,yd) -- Read data
xc,c = {0,6,9,11,16,20},{} -- Set up location of knots
c[1],c[6] = 0,0.97 -- Specify two end points
fxc = datafit({yd,xd},{xc,c}) -- Now do fit with floating knots
for j=1,#xd do ycalc[j] = fxc(xd[j]) end -- Calculate fit values
plot(xd,yd,ycalc) -- Plot data and fit
write_data("list7_18.dat",xd,yd,ycalc) -- Save data 

