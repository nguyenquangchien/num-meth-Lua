-- /* File list7_19.lua */  Fit to gauss3.dat with selected knots

require "DataFit"

xd,yd,xcalc,ycalc={},{},{},{}
read_data('Gauss3.dat',yd,xd) -- Read data
xc = {0,35,70,75,85,95,103,112,117,125,135,140,148,
162,180,190,210,250}
c = {}
fxc = datafit({yd,xd},{xc,c})
for i=1,#xd do ycalc[i] = fxc(xd[i]) end
plot(xd,yd,ycalc)
write_data("list7_19.dat",xd,yd,ycalc) -- Save data 

