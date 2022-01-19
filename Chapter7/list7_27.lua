-- /* File list7_27.lua */ -- Some tests of elementary functions

require"elemfunc"
gamma = elemfunc.gamma -- Gamma function
prob = elemfunc.Pn -- Probability function
J0 = elemfunc.J0 -- Bessel function
i,xt,y1t,y2t,y3t = 1,{},{},{},{} -- Now test function
for x=-5,5.0001,.01 do 
	xt[i] = x
	y1t[i] = gamma(x)
	y2t[i] = prob(x)
	y3t[i] = J0(x)
	i = i+1
end	
plot(xt,y1t,{'gamma plot','x','y','set yrange [-5:5]'})
plot(xt,y2t); plot(xt,y3t)
write_data("list7_27.dat",xt,y1t,y2t,y3t)
