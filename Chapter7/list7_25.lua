-- /* File list7_25.lua */ -- Rational approximation to prob fumction
require"intg"; require"nlstsq"

x1,y1,x2,y2 = {},{},{},{}
f1 = function(x) return math.exp(-x^2/2) end -- Normal curve
f2 = function(x) return f1(1/x)/x^2 end -- Folded normal curve
x1[1],y1[1],i = 0,0,2
for xi=.05,1.01,.05 do -- Integral from 0 to 1
	x1[i],y1[i] = xi, intg(xi-.05,xi,f1)/math.sqrt(2*math.pi)+y1[i-1]
	i = i+1
end
x2[1],y2[1],i = 1,y1[#x1],2
for xi=.95,-.01,-.05 do -- Integral from 1 to infinity
	x2[i],y2[i] = xi, intg(xi,xi+.05,f2)/math.sqrt(2*math.pi)+y2[i-1]
	i = i+1
end
n = #x2; j = n+1
for i=1,n do y1[i],y2[i] = .5*y1[i]/y2[n]+.5,.5*y2[i]/y2[n]+.5 end -- Exactly .5 max
for i=2,n-1 do -- Combine results into table from x=0 to x=20
	x1[j],y1[j] = 1/x2[i],y2[i] 
	j = j+1
end -- x1,y1 now contains P(x) values

fft = function(yx,c) -- Function to fit data
	x = yx[2]
	return 1. - 0.5/(1+c[1]*x+c[2]*x^2+c[3]*x^3+c[4]*x^4+
		c[5]*x^5+c[6]*x^6)^16 - yx[1]
end

c = {0,0,0,0,0,0} -- Six possible coefficients
yx = {y1,x1} -- Data values

del,err,nmax = nlstsq(yx,fw,fft,c) -- Call fiting
print(err,nmax) -- print results
for i=1,#c do printf('c[%d] = %14.6e +/- %12.4e\n',i,c[i],del[i]) end
