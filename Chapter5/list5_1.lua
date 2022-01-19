function deriv(f,x)
	local dx = x*FACT
	if dx==0 then dx = FACT end
	return (f(x+dx) - f(x-dx))/(2*dx)
end
setfenv(deriv,{FACT=1.e-6})

f1 = function(x) return x^5 end
f2 = math.exp
print(deriv(f1,1.e8),5*(1.e8)^4)
print(deriv(f2,10),math.exp(10))

