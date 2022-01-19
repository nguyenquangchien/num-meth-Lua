-- /* list11_29.lua */
-- Solution of first order diff. eqn. with FD method
require"odefd"
getfenv(odefd).nprint=1
getfenv(odefd).umin = {1.e-5,1.e-5}

a = 1 -- or 0.5 or 1 or 4 or any number
f1 = function(eqs,x,u,up,upp)
	eqs[1] = u[1]*up[1] + (1-a)*u[1] - a
end
f2 = function(eqs,x,u,up,upp)
	eqs[1] = u[1]*upp[1] +up[1]^2 + (1-a)*up[1]
end
fbl = function(eql,u,up)
	eql[1] = u[1] 
end
fbr = function(eqr,u,up)
	eqr[1] = u[1]*up[1] +(1-a)*u[1] - a
end

x = {0,100,2000}; u = {{1,500}}
s1,er1,n1 = odefde({f1,fbl,fbr},x,u)
x = {0,100,2000}; u = {{1,500}}
s2,er2,n2 = odefde({f2,fbl,fbr},x,u)
plot(s1,er1); plot(s2,er2)
print('Number of iterations =',n1,n2)
write_data('list11_29.dat',s1,s2,er1,er2)