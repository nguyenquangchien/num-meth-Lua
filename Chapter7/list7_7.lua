-- /* File list7_7.lua */ -- data fitting for Figure 7.1

require"Fourier"

xd,yd = {},{}
read_data('list7_1.dat',xd,yd) -- Data for Figure 7.1

n1 = #xd
i = n1+1
for j=n1-1,1,-1 do -- Expand to half cycle
	xd[i],yd[i] = xd[j+1]-xd[j]+xd[i-1], yd[j]
	i = i+1
end
n1 = #xd
for j=2,n1 do -- Expand to full cycle
	xd[i] = xd[i-1]+xd[j+1]-xd[j]
	yd[i] = -yd[j]
	i = i+1
end
fx,fy,fz = fourier(yd,40) -- Fourier coefficients -- 40 components
plotFourier(fx) -- plot Fourier coefficients

c1,c3 = fx[2][1],fx[4][1]
print('c1,c3 = ',c1,c3)
ycalc = {}
for i=1,n1 do
	ycalc[i] = c1*math.sin(2*math.pi*xd[i]/80) + c3*math.sin(6*math.pi*xd[i]/80)
end
write_data('list7_7.dat',xd,yd,ycalc)
