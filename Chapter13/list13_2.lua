-- File list13_2.lua --

require"spgauss"; require'pde2fe' 

nodes,triangles,sides = readnts('list13_2')

nnds,nel = #nodes, #triangles
print('number nodes =',nnds)
print('number elements =',nel)

pi = math.pi; pi2 = 2*pi^2; sin = math.sin
function feq(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
	return uxx + uyy  + pi2*(sin(pi*x)*sin(pi*y))
end
function fb(nd,u,un,nbs)
	return u
end

u = {}; for i=1,nnds do u[i] = 0.0 end
a,b = setup2feeqs({feq,fb},{nodes,triangles,sides},u)

getfenv(spgauss).nprint = 1
spgauss(a,b)

x,y = {},{}; NT = 21
for j=1,NT do x[j] = (j-1)/(NT-1); y[j] = x[j] end
sol,solxy = toxysol(b,x,y,nodes,triangles)

write_data('list13_2a.dat',sol);write_data('list13_2b.dat',solxy)
splot(solxy); cplot(solxy)