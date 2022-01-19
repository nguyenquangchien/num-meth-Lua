--/* File nsolv.lua */
require"gauss"; require"spgauss" -- Make sure gauss and spgauss are loaded

nsolv = function(f,x,step,nx) 
	local nend,wk1,rowptr,b = NMAX,{},{},{} -- work, A, B
	local c1,c2,dx,cmax,ngauss,stp
	if full then ngauss=gauss else ngauss=spgauss end
	if nx==nil then f(b,x); nx = #b end -- Find number of equations
	for ilp=1,NMAX do -- iterate for at most NMAX steps
		for i=1,nx do rowptr[i] = {} end
		f(b,x) -- Values of equations returned in b[]
		for i=1,nx do -- Set up matrix equations
			c1 = x[i]; dx = FACT*c1
			if dx==0.0 then dx = FACT end
			x[i] = x[i]+dx -- Probe x[i] factor
			f(wk1,x); x[i] = c1 -- Results returned in wk1[]
			for j=1,nx do 
				c1 = b[j]-wk1[j]
				if full then rowptr[j][i] = c1/dx
				elseif c1~=0.0 then rowptr[j][i] = c1/dx end
			end
		end
		ngauss(rowptr,b,nx) -- Solve for correction terms -- Solution returned in b[]
		cmax = 0.0
		for i=1,nx do
			if step~=nil and step[i]~=0.0 then -- Skip if no step sizes specified
				stp = step[i]
				if stp<0.0 then -- Limits on absolute change
					if abs(b[i])>-stp then 
					b[i] = -stp*abs(b[i])/b[i] end
				else -- Limits on percentage change
					if stp<1 then stp=1/stp end
					c1 = 1. + b[i]/x[i]
					if c1>stp then 
						b[i] = (stp-1)*x[i] 
					elseif c1<1/stp then 
						b[i] = (1/stp-1)*x[i] 
					end
				end
			end -- End of step{} adjustments, find largest relative error
			c1 = abs(b[i]); c2 = abs(x[i]) + c1
			if c2~=0.0 then c1=c1/c2 end
			if c1>cmax then cmax=c1 end
			x[i] = x[i] + b[i] -- Update solution set
		end
		if linear then nend=ilp; break end -- No iterations needed if linear equations
		if nprint then -- Print iterative solutions
			fprintf(stderr,'Solutions at iteration %d are:\n',ilp)
			for i=1,nx do fprintf(stderr,'%s  ',tostring(x[i])) end
			fprintf(stderr,'\n'); stderr:flush()
		end
		if cmax < ERROR then nend=ilp; break end
	end
	return nend,cmax -- Return number of iterations and error, solution returned in x array
end
setfenv(nsolv,{stderr=io.stderr,gauss=gauss,spgauss=spgauss,tostring=tostring,fprintf=fprintf,
	abs=(Complex or math).abs,ERROR=2.e-6,FACT=1.e-6,NMAX=100,nprint=nil,linear=nil,full=nil})
