-- /* File list5_8.lua */ -- Simple integration 

require'intg'

fx = function(x) return x^.5 end
iv,n = intg(0,4,fx)
ivv = 2*4^1.5/3
print(iv,n,ivv,(iv-ivv)/ivv)

fy = function(x) return math.sin(x) end
iv,n = intg(0,4,fy)
ivv = 1-math.cos(4)
print(iv,n,ivv,(iv-ivv)/ivv)

humps = function(x) -- Has 2 maxima and 1 minimum between 0 and 1
	return 1/((x-.3)^2+.01) + 1/((x-.9)^2+.04) -6
end	
t1 = os.clock()
for k=1,100 do
iv,n = intg(0,4,humps)
end
print('time =',os.clock()-t1)
ivv = 10*(math.atan(37)-math.atan(-3))+5*(math.atan(15.5)-
	math.atan(-4.5))-24
print(iv,n,ivv,(iv-ivv)/ivv)



