-- /* File list11_16.lua */
-- Shooting method for boundary value problem with nonlinear DE
require"odebvfd" -- Use ode1fd() function

nx,L,ubL,ubR = 1000, 1,0, 2 -- #points, Left, Right Boundary values

f = function(x,u,up,upp) -- Differntial equation
	return upp + 4*up^2
end	
-- x range values, boundary values specified
x,u = {0,L,nx},{ubL,ubR}
u,nn,err1,err2 = ode1fd(f,x,u) -- Solve BV problem

print(nn,err1,err2); plot(u) 
C2,uex = math.exp(8),{} -- Exact solution parameters
for i=1,nx do uex[i] = math.abs(0.25*math.log(u[1][i]*(C2-1)+1) - u[2][i]) end
write_data("list11_16.dat",u,uex)

