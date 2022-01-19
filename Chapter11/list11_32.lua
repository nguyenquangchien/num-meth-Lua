-- /* list11_32.lua */
-- Fitting of data to BV problem for parameter estimation 
require"odefd"; require"nlstsq"; require"intp"; require"DataFit"
exp = math.exp
x,u1,u2 = {},{},{}
read_data('list11_32_in.dat',x,u1,u2)
_,_,std1 = datafit{u1,x}; _,_,std2 = datafit{u2,x} -- Est. standard dev.
print('std1, std2 = ',std1, std2)
nd = #x -- number of data points
utp = {1}; for i=2,nd do utp[i] = 2 end -- flag for u data 
u = {}; for i=1,nd do u[i] = u1[i] end
for i=1,nd do
	x[i+nd] = x[i] -- repeat x values
	u[i+nd] = u2[i] -- make single array
	utp[i+nd] = 3 -- flag for u2 data
end -- 3 data values, one of which is data flag
yx = {u,x,utp} -- data arrays, dependent values first
c = {0.5, 0.5e8, 22000, .5e11} -- Guess at Differential eqn.  coefficients
nt = 200 -- Number of points for integration
nx,nc,xx = #yx, #c, {0,50,nt} 
nd = #x -- 2*number of data points
actv,step,del,fw = {},{},{},{} -- nlstsq parameters
for i=1,nd do 
	if i<= nd/2 then fw[i] = 1/std1^2 
	else fw[i] = 1/std2^2 end -- 1/std^2 weighting 
	--fw[i] = 1 -- Unity weighting -- Try these
	--fw[i] = 1/u[i] -- 1/Y weighting
	--fw[i] = 1/u[i]^2 -- 1/y^2 weighting
end 
ce, te = 0.07, 1250
ui = {{0,0},{1250,1250}} -- Initial approximations to solutions
actv = {0,1,1,0} -- Tru any two combinations
ft = function(eqt,x,u,up,upp)  -- Differential equation
	eqt[1] = c[1]*upp[1] - up[1] - c[2]*u[1]^2*exp(-c[3]/u[2])
	eqt[2] = c[1]*upp[2] - up[2] + c[4]*u[1]^2*exp(-c[3]/u[2])
end
fbl = function(eqb,u,up) -- Left boundary conditions
	eqb[1] = u[1] - c[1]*up[1] - ce
	eqb[2] = u[2] - c[1]*up[2] - te
end
fbr = function(eqb,u,up) -- Right boundary conditions
	eqb[1] = up[1]; eqb[2] = up[2] -- Zero derivatives
end
feq = function(yx,c,new) -- Equation for nlstsq()
	flg = yx[3]
	if flg==1 then
		s = odefd({ft,fbl,fbr},xx,ui) -- Solve differential equation
		xt, flg = s[1], 2 
		ff = {intpf(s[1],s[2]), intpf(s[1],s[3])}
	end -- New solution now available
	return ff[flg-1](yx[2]) - yx[1]
end
del, err, nmax = nlstsq(yx,fw,feq,c,actv,step)
print('RMS error =',err)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.3e\n',i,c[i],del[i]) end
s = odefd({ft,fbl,fbr},xx,ui) -- Solve differential equation
write_data('list11_32.dat',s,x,u); plot({s[1],s[2]},{x,u1}); plot({s[1],s[3]},{x,u2})
