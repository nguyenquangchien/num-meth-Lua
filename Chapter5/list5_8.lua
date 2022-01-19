-- /* File list5_8.lua */ -- Simple integration 

function sintg(xmin,xmax,f)
	local n = 4000
	local dx = (xmax-xmin)/n
	local sum = 0.5*(f(xmin)+f(xmax))
	for i=1,n-1 do sum = sum + f(xmin+i*dx) end
	return sum*dx
end

fx = function(x) return 25*x^.5 end
iv = sintg(0,4,fx)
ivv = 25*4^1.5/1.5
print(iv,ivv,(iv-ivv)/ivv)

fy = function(x) return 50*math.sin(x) end
iv = sintg(0,4,fy)
ivv = 50*(1-math.cos(4))
print(iv,ivv,(iv-ivv)/ivv)

humps = function(x) -- Has 2 maxima and 1 minimum between 0 and 1
	return 1/((x-.3)^2+.01) + 1/((x-.9)^2+.04) -6
end	
iv = sintg(0,4,humps)
ivv = 10*(math.atan(37)-math.atan(-3))+5*(math.atan(15.5)-
	math.atan(-4.5))-24
print(iv,ivv,(iv-ivv)/ivv)



