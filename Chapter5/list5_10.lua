-- /* File list5_10.lua */ -- Simple integration 

require"intg"

fx = function(x) return 25*x^.5 end
iv,n,err = intg(0,4,fx)
ivv = 25*4^1.5/1.5
print(n,iv,ivv,err,(iv-ivv)/ivv)

fy = function(x) return 50*math.sin(x) end
iv,n,err = intg(0,4,fy)
ivv = 50*(1-math.cos(4))
print(n,iv,ivv,err,(iv-ivv)/ivv)

humps = function(x) -- Has 2 maxima and 1 minimum between 0 and 1
	return 1/((x-.3)^2+.01) + 1/((x-.9)^2+.04) 
end	
iv,n,err = intg(0,4,humps)
ivv = 10*(math.atan(37)-math.atan(-3))+5*(math.atan(15.5)-
	math.atan(-4.5))
print(n,iv,ivv,err,(iv-ivv)/ivv)

fz = function(x) return x*math.exp(-x) end
iv,n,err = intg(0,4,fz)
ivv = 1 - 5*math.exp(-4)
print(n,iv,ivv,err,(iv-ivv)/ivv)

