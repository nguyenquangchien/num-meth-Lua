-- /* File list10_2.lua */
-- Programs to explore stability and accuracy for first order differential equations

require"nsolv"; exp = math.exp
typ = 'TP' -- Trapeziodal rule -- Sekect as desured
--typ = 'BD' -- Backwards difference
--typ = 'FD' -- Forward difference

feqs = function(eq, t, y, yp) -- Define differential equations
	eq[1] = yp[1] + 1001*y[1] - 999*y[2]
	eq[2] = yp[2] - 999*y[1] + 1001*y[2]
end

h,kmax = 1.e-3, 1000 -- Select desired time scale
--h,kmax = 1.e-5, 1000
neq,t = 2,0

y = {}; sol = {{t},{1},{1},{},{}} -- Solution array with Initial values
for i=1,neq do y[i] = sol[i+1][1] end -- initial y value array

fderiv = function(eqs,yp) -- Function for calculating derivatives - yp
	feqs(eqs,t,y,yp) -- Add t and y to arguments
end

fnext = function(eqs,y) -- Function for next y values
	local yp,h2 = {}
	if typ=='TP' then h2=h/2 else h2=h end
	for i=1,neq do -- TP or BD algorithms
		yp[i] = (y[i] - yn[i])/h2 -- trapezoidal rule
	end
	feqs(eqs,t,y,yp) -- Add t and yp to arguments
end

yp,yn = {},{}; for i=1,neq do yp[i] = 0 end -- Derivative and step array 

for k=1,kmax do -- Main time loop
	nsolv(fderiv,yp) -- Update y' values
	if typ=='TP' then  for i=1,neq do yn[i] = y[i] + 0.5*h*yp[i] end  -- Trapezoidal rule
	else for i=1,neq do yn[i] = y[i] end end -- Backwards differencing
	for i=1,neq do y[i] = y[i] + h*yp[i] end --Predict new value, final for Forward diff
	t = t+h -- Update time 
	if typ~='FD' then nsolv(fnext,y) end -- Calculate new values
	sol[1][k+1] = t -- Save calculated time values
	for i=1,neq do sol[i+1][k+1] = y[i] end
end
for i=1,#sol[1] do 
	f1,f2 = exp(-2*(i-1)*h), exp(-2*(i-1)*h*1000)
	sol[neq+2][i],sol[neq+3][i] = f1 + f2, f1 - f2
end
write_data('list10.2'..typ..'.dat',sol); plot(sol)