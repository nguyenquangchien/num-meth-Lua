-- /* File list6_7.lua */ -- Function inversion by Interpolation 
require"intp"
  
x,i = 0,1
xd,yd = {},{}
while x<18. do -- Generate table of values
	xd[i] = x
 	yd[i] = math.log(1+x^2)/(1+.2*math.atan(x))
 	i,x = i+1, x+.5
end
fx = intpf(yd,xd) -- Define interpolation function
xf,yv = {},{}
i=1
for y=0,4,.01 do
	yv[i] = y
	xf[i] = fx(y) -- Use the new inverse function
	i = i+1
end
print(fx(1),fx(2),fx(3)) -- Print selected values 	
write_data('test1.dat',xd,yd)
write_data('test2.dat',yv,xf)
plot(yv,xf)
