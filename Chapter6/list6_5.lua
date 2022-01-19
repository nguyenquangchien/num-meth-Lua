-- /* File list6_5.lua */ -- Tests of interpolation functions

require"intp"; require"spline"; require"deriv"
makeglobal(math) -- So we can use exp() in place of math.exp()

xd,yd = {},{} -- Arrays for tabular data
f1 = function(x) return  exp(-x^2) end -- Five test functions
f2 = function(x) return x*exp(-x) end
f3 = function(x) return 1/(1+25*x^2) end
f4 = function(x) return log(x) end
f5 = function(x) return sin(x) end
f1d = function(x) return -2*x*exp(-x^2) end
f2d = function(x) return (1-x)*exp(-x) end
f3d = function(x) return -50*x/(1+25*x^2)^2 end
f4d = function(x) return 1/x end
f5d = function(x) return cos(x) end
-- Select function to test -- Change for different functions
--funct,functd,xmin,xmax,dx1,dx2 = f1,f1d,-2,2,0.5,0.01 -- f1 tests
--funct,functd,xmin,xmax,dx1,dx2 = f2,f2d,0,5,.5,.01 -- f2 tests
funct,functd,xmin,xmax,dx1,dx2 = f3,f3d,-1,1,.25,.01 -- f3 tests
--funct,functd,xmin,xmax,dx1,dx2 = f4,f4d,0.5,5,.5,.01 -- f4 tests
--funct,functd,xmin,xmax,dx1,dx2 = f5,f5d,0,2*pi,.5,.01 -- f5 tests
i,x = 1,xmin
while x<=xmax+.001 do -- Define tabular data
	xd[i],yd[i] = x, funct(x)
	x,i = x + dx1, i+1
end -- Select one of following two lines
--fintp,fspline = intpf(xd,yd), splinef(xd,yd)
fintp,fspline = intp4f(xd,yd), splinef(xd,yd); getfenv(spline).stype = 2

xi,y_lcb,y_csp,y_fun={},{},{},{} -- Tables for function and interpolated data
y1_lcb,y1_csp,y1_fun,y2_lcb,y2_csp,y2_fun = {},{},{},{},{},{} -- Derivative values
x = xmin; i = 1
while x<=xmax+.001 do -- Generate interpolated data values
	xi[i] = x
	y_fun[i],y1_fun[i],y2_fun[i] = funct(x), functd(x), deriv(functd,x)
	y_lcb[i],y1_lcb[i],y2_lcb[i] = fintp(x), fintp(x,1), fintp(x,2)
	y_csp[i],y1_csp[i],y2_csp[i] = fspline(x), fspline(x,1),fspline(x,2)
	x = x+dx2; i = i+1
end
ndat = #xi
er0lcb,er0csp,er1lcb,er1csp,er2lcb,er2csp = 0,0,0,0,0,0
for i=1,ndat do -- Calculate average errors
	er0lcb = er0lcb + (y_fun[i] - y_lcb[i])^2
	er0csp = er0csp + (y_fun[i]-y_csp[i])^2
	er1lcb = er1lcb + (y1_fun[i]-y1_lcb[i])^2
	er1csp = er1csp + (y1_fun[i]-y1_csp[i])^2
	er2lcb = er2lcb + (y2_fun[i]-y2_lcb[i])^2
	er2csp = er2csp + (y2_fun[i]-y2_csp[i])^2
end
print("Average mean square errors are:")
print("Functions, LCB, CSP = ",sqrt(er0lcb/ndat), sqrt(er0csp/ndat))
print("First derivatives, LCB, CSP = ",sqrt(er1lcb/ndat),sqrt(er1csp/ndat))
print("Second derivatives, LCB, CSP = ",sqrt(er2lcb/ndat),sqrt(er2csp/ndat))
--write_data("intp_org.dat",xd,yd)
write_data("list6_5.dat",xi,y_fun,y_lcb,y_csp,y1_fun,y1_lcb,y1_csp,
	y2_fun,y2_lcb,y2_csp,xd,yd)
plot(xi,y_fun,y_lcb,y_csp)


