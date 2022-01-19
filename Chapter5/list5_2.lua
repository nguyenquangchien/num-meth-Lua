-- /* File list5_2a.lua */ -- Function x^10 and x = 1e10 and 1e-10
f1= function(x)
	return math.exp(x)
end

n2 = 5
n2 = 10
f2 = function(x)
	return x^n2
end


a = 1

sa = {}
sf1 = {}
ssf1 = {}
sf2 = {}
ssf2 = {}
f = 10^0.1
k = 1
while a>1.e-14 do
		a = a/f
		x1 = 1e-10
		x2 = 1e10
		dx1 = x1*a
		dx2 = x2*a
		ssf1[k] = n2*x1^(n2-1)
		ssf2[k] = n2*x2^(n2-1)
		--df1 =((f2(x1+dx1) - f2(x1))/(dx1))/ssf1[k] -1
		--df2 = (f2(x2+dx2)-f2(x2))/(dx2)/ssf2[k] -1
		df1 =((f2(x1+dx1) - f2(x1-dx1))/(2*dx1))/ssf1[k] -1
		df2 = (f2(x2+dx2)-f2(x2-dx2))/(2*dx2)/ssf2[k] -1
		sa[k] = a
		sf1[k] = math.abs(df1)
		sf2[k] = math.abs(df2)
		k = k+1
end

write_data('list5_2d.dat',sa,sf1,ssf1,sf2,ssf2)		
		