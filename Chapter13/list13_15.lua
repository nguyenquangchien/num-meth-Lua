-- File list13_15.lua -- Velocity in air duct

require'pde2fe' 
nts ={readnts('list13_15')} -- square area, 2996 nodes 
Lx,Ly = 2,1 -- Ratio of lengths squared

feqs = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		return uxx/Lx^2 + uyy/Ly^2  + 1.0 
	end,
	function(nd,u,un,nbs)
		return u
	end
}

getfenv(pde2bvcoe).nprint=1 -- Observe progress
getfenv(pde2fe).nprint = 1
SPM,COE,SOR = 1,2,3 -- Try different methods

u = pde2fe(feqs,nts,u,COE) -- Solve equations

x,y = {},{}; NT = 41-- x,y arrays for plotting
for j=1,NT do x[j] = (j-1)/(NT-1); y[j] = x[j] end
sol,soxy = toxysol(u,x,y,nts[1],nts[2])
write_data('list13_15.dat',sol)
splot(soxy);cplot(soxy) -- popup plots
