-- /* File list5_8.lua */ -- Simple integration 

function sintg(xmin,xmax,f)
	--local n = 4000
	local dx = (xmax-xmin)/n
	local sum = 0.5*(f(xmin)+f(xmax))
	for i=1,n-1 do sum = sum + f(xmin+i*dx) end
	return sum*dx
end
iv1 = 25*4^1.5/1.5
iv2 = 50*(1-math.cos(4))
iv3 = 10*(math.atan(37)-math.atan(-3))+5*(math.atan(15.5)-
	math.atan(-4.5))-24
	
fx = function(x) return 25*x^.5 end
fy = function(x) return 50*math.sin(x) end
humps = function(x) -- Has 2 maxima and 1 minimum between 0 and 1
	return 1/((x-.3)^2+.01) + 1/((x-.9)^2+.04) -6
end	
	

n = 40; dx = 4/n; i = 1
y = {{},{},{}}; dy = {{},{},{}}
FCT = 2	
while dx>3.e-8 do
	y[1][i] = sintg(0,4,fx)
	dy[1][i] = (y[1][i]-iv1)/iv1
	y[2][i] = sintg(0,4,fy)
	dy[2][i] = (y[2][i]-iv2)/iv2
	y[3][i] = sintg(0,4,humps)
	dy[3][i] = (y[3][i]-iv3)/iv3
	n = FCT*n; dx = 4/n; i = i+1
	print('n =',n,dx);io.flush()
end
	
write_data('list5_8a.dat',y,dy)




