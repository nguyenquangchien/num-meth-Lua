-- /* File list6_2.lua */ -- Linear interpolation using lintpf()
	
require'intp'
xd = {0,1,2,3,4,5}; yd = {3,4,4.5,4,3,2.5}
ff = lintpf(xd,yd) --- Convert to function

x,y = {},{}
i,xv =1, -1
while xv<6.001 do
	x[i],y[i] = xv, ff(xv) -- Call function
	i,xv  = i+1, xv+.01
end 
write_data('list6_2.dat',x,y); plot(x,y)