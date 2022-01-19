-- /* File list7_22.lua */ -- datafit() to transistor data

require "DataFit"
xd,yd={},{{},{},{},{},{},{},{},{},{}}
nd = read_data('mosidd1.5.0.6.dat',xd,yd)
xx,yy = {},{{},{},{},{},{},{},{},{},{}}
xmin,xmax,nx = xd[1],xd[nd],1000
dx = (xmax-xmin)/nx
for i=1,nx do xx[i] = xmin+dx*(i-1) end
ncrv,fxc = #yd, {}
for j=1,ncrv do
	fxc[j] = datafit{yd[j],xd}
end
for j=1,ncrv do	
	for i=1,nx do yy[j][i] = fxc[j](xx[i]) end
end
plot(xx,unpack(yy))
write_data("list7_22.dat",xx,yy) -- Save data 

