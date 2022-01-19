-- /* File list9_13.lua */ -- Data fitting, Parameter estimation for transfer function

requires("nlstsq","Polynomial")

-- Section 1. Input data file and define fitting function
infile = 'amplifier2' -- or 'amplifier3'
xd,yd,_ = {},{},{}; read_data(infile..'.dat',xd,_,_,yd) -- columns 1 and 4 of fule

c1,f1,f2,f3,f4,f5,f6,f7 = 22,70,7,10,300,3e7,1e8,1e9 -- Initial guesses
c1 = c1*f1/(f2*f3*f4)
p1 = Polynomial.new{1,1/f2}*{1,1/f3}*{1,1/f4}
p2 = Polynomial.new{1,1/f5}*{1,1/f6}
c = {c1,1/f1,p1[2],p1[3],p1[4],p2[2],p2[3],1/f7}; nc = #c

ft = function(x,c) -- Define function to fit data
	w = j*x[2]
	h = c[1]*w^2*(1+c[2]*w)/(1+c[3]*w+c[4]*w^2+c[5]*w^3)*
		(1+c[8]*w)/(1+c[6]*w+c[7]*w^2)
	return  20*math.log10(Complex.abs(h)) - x[1]
end

--- Section 2. Perform data fit and print model parameters
del,err,nx =nlstsq({yd,xd},fw,ft,c) -- Call fiting, print max iterations 
print('Number of iterations,err =',nx,err)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
print('Zeros of Transfer Function at (in Hz):')
print(-1/c[2]); print(-1/c[6])
roots = Polynomial.roots{1,c[3],c[4],c[5]}
print('Poles of Transfer Function at (in Hz):')
print(roots[1]); print(roots[2]); print(roots[3])
roots = Polynomial.roots{1,c[6],c[7]}
print(roots[1]); print(roots[2])
