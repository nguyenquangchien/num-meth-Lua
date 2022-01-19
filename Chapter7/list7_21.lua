-- /* File list7_21.lua */ -- Example of datafitn() function

require "DataFit"

xd,yd,xcalc,ycalc={},{},{},{}
read_data('thurber.dat',yd,xd) -- Read data
fxc = datafitn({yd,xd},1,5)
nd = #xd
for i=1,nd do ycalc[i] = fxc(xd[i]) end
xx,yy = {},{}
xmin,xmax,nx = xd[1],xd[nd],1000
dx = (xmax-xmin)/nx
for i=1,nx do 
	xx[i] = xmin+dx*(i-1)
	yy[i] = fxc(xx[i])
end
plot(xx,yy)
write_data("list7_21.dat",xx,yy,xd,yd) -- Save data 
