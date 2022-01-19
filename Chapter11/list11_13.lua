-- /* File list11_13.lua */
-- Eigenvalue problem for screened Coulomb potential
require"odebv";require"intp"

f1 = function(eqs,E,x,u,up,upp) -- Hydrogen atom
	eqs[1] = upp[1] + (2*math.exp(-x/Ld)/x + E[1] - L*(L+1)/x^2)*u[1]
end

L = 0; Ld = 2
nx,xmax,Ei = 2000, 20, -1
E = {Ei}

s = odebvev({f1,E},{0,xmax,nx},{0},{1})
print('Energy eigenvalue = ', E[1])
sn,sqn = sqnorm(s,1.0)
plot(sn); plot(sqn) 
write_data('list11_13.dat',s,sn,sqn)



