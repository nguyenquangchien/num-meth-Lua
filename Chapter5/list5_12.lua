-- /* File list5_12.lua */ Test of integrals with infinite limit

require"intg"; makeglobal(math)
getfenv(intg).ERROR=1.e-8 -- Try commenting out

f1 = function(x) -- First test integral --> sqrt(pi)/2
	return exp(-x^2)
end
f2 = function(x) -- Second test integral --> pi
	return 2/(1+x^2)
end
f3 = function(x) -- Third test integral --> 2.0
		return x^2*exp(-x)
end
f4 = function(x) -- Fourth test integral --> -pi*ln(10)/20
	return log(x)/(1 + 100*x^2)
end	
f5 = function(x) -- Fifth test integral --> .0674467742
	return x^.2/(1 + 10*x)^2
end	

a,b,c,d,e = sqrt(pi)/2,pi,2.0, -pi*log(10)/20, -10^(-6/5)*pi/(5*sin(6*pi/5))
aa,ab,ac = intg_inf(f1), intg_inf(f2), intg_inf(f3)
ad,ae = intg_inf(f4), intg_inf(f5)
print('value of f1 = ',aa,(a-aa)/a)
print('value of f2 = ',ab,(b-ab)/b)
print('value of f3 = ',ac,(c-ac)/c)
print('value of f4 = ',ad,(d-ad)/d)
print('value of f5 = ',ae,(e-ae)/e)
