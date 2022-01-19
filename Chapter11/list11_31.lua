-- /* list11_31.lua */
-- Estimation of parameters for first order differential equation
require"odeiv"; require"nlstsq"; require"intp"
getfenv(nlstsq).ylinear = 1 -- Not required 

t,ct = {},{}; read_data('conc.txt',t,ct) -- Data points

c = {.3, 1.0, 24} -- Guess at Differential eqn.  coefficients
nt = 200 -- Number of time intervals for integration
xmin,xmax,nd = 0.0, 20, #t -- number of data points
xx = {xmin, xmax, nt} -- nt time intervals
step = {2,1.1,0} -- nlstsqi parameters
--actv = {1, 0, 1} -- Fix exponent at 1.0 -- Try it !!!

ft = function(eqt,x,u,up)  -- Differential equation
	eqt[1] = up[1] + c[1]*u[1]^c[2]
end

feq = function(yx,c,id) -- Equation for nlstsq()
	xv = yx[2]
	if xv==xmin then -- Some coefficient has changed
		s = odeiv(ft,xx,{c[3]}) -- Solve differential equation
		xt,con = s[1],s[2]
	end -- New solution now available
	return intp(xt,con,xv) - yx[1] -- Return solved DE value
end
del,err,nmax = nlstsq({ct,t},fw,feq,c,actv,step) 
for i=1,#c do printf('c[%d] = %12.4e +/- %12.3e\n',i,c[i],del[i]) end
s = odeiv(ft,xx,{c[3]}) -- Final solution of DE
write_data('list11_31.dat',s,t,ct)
plot(s,{t,ct})
