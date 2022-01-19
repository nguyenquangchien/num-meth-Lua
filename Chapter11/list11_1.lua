-- /* File list11_1.lua */
-- Shooting method for boundary value problem

require"odeiv"; require"newton"

-- Parameters
T,E,I,w,L = 500, 1.e7, 500, 100, 100
EI = E*I; print('w,EI =',w,EI)
y0, yL,nx = 0, 1.e-15, 2000 -- Boundary values,  #x values

f = function(eq,x,u,up) -- Differntial equations
	eq[1] = up[1] - u[2]
	eq[2] = up[2] - T*u[1]/EI - w*x*(x-L)/(2*EI)
end

bvalue = function(up) -- Boundary function, up is derivative at x = 0
	s1 = odeiv(f,{0,L,nx},{0,up}) -- Solve initial value problem
	return s1[2][nx+1] - yL -- Return difference in boundary value
end

yp = 0 -- Initial guess at derivative
yp,nm,err = newton(bvalue,yp) -- Use Newton's method to satisfy boundary value
print('Initial derivative, #iterations, errors =',yp,nm,err,bvalue(yp))
s2 = odeiv(f,{0,L,nx},{0,yp*1.1}) -- Larger derivative
s3 = odeiv(f,{0,L,nx},{0,yp/1.1}) -- Smaller derivative
plot(s1[1],s1[2],s2[2],s3[2])
write_data("list11_1.dat",s1,s2,s3)

