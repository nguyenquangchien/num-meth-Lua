 --/* File list3_11.lua */ -- Example of periodic function
 
require "newton"

function f(x) return .7 - math.sin(x*math.pi/4)^4 end

i = 1;x,y = {},{}
for xx=0,10,.002 do 
	x[i],y[i] = xx,f(xx)
	i = i+1
end
plot(x,y); write_data('list3.11.dat',x,y)	
print(newton(f,0.1))
print(newton(f,1))
print(newton(f,3))