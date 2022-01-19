-- /* File list5_4.lua */ -- Finding maximum and minimum values

require"deriv"; require"newton"

pos = 6.3324567
f1 = function(x) return math.exp((x-pos)^2) end
f2 = function(x) return -(x+pos)^4 end
f3 = function(x) return 5*x^4-4*x^2+2*x-3 end

fderiv = function(f)
	return function(x) return deriv(f,x) end
end
f1deriv = fderiv(f1)

sol,n,err = newton(f1deriv,0)
print(sol,n,err,f1(sol),(sol-pos)/pos)
if f1(sol)<f1(1.01*sol) then print("solution is a minimum at ",sol)
else print("solution is a maximum at ",sol) end

sol,n,err = newton(fderiv(f2),0)
print(sol,n,err,f2(sol),(sol+pos)/pos)
if f2(sol)<f2(1.01*sol) then print("solution is a minimum at ",sol)
else print("solution is a maximum at ",sol) end

sol,n,err = newton(fderiv(f3),0)
print(sol,n,err,f3(sol))
if f2(sol)<f2(1.01*sol) then print("solution is a minimum at ",sol)
else print("solution is a maximum at ",sol) end
