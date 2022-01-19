-- /* File list7_17.lua */ -- Simple fit to data using datafit()

require "DataFit"

xd,yd,ycalc={},{},{}
read_data('list7_1.dat',xd,yd) -- Read data
fxc = datafit{yd,xd} -- Now do fit with floating knots
for j=1,#xd do ycalc[j] = fxc(xd[j]) end -- Calculate fit values
plot(xd,yd,ycalc) -- Plot data and fit
write_data("list7_17.dat",xd,yd,ycalc) -- Save data 

