--/* File nlstsq.lua */

require "gauss"; require"Minv" -- Matrix inversion needed
nlstsq = function(x, fw, f, c, actv, step)
	actv,step,fw = actv or {},step or {},fw or {}
	local nx,nd,nc = #x,#x[1],#c
	local iend,nend,cv,cmax,err = 0, NMAX -- Maximum Number of iterations
	local del,p,wk1,wk2,rowptr,b,fx = {},{},{},{},{},{},{} -- Work arrays
	for i=1,nc do 
		rowptr[i],step[i],actv[i] = {},step[i] or 0,actv[i] or 1 -- Row pointers
		if step[i]>1.0 then step[i] = 1./step[i] end
	end
	for ilp=1,NMAX do -- Main loop for Newton iterations
		ix = 1
		for i=1,nc do -- Loop over C coefficients
			if actv[i]~=0 then -- An active coefficient
				b[ix],row,wk2[ix] = 0.0, rowptr[ix], {}
				for j=1,nc do row[j] = 0.0 end -- Initialize row of matrix to zero
				dx = FACT*c[i] -- Set up dC factors for partial derivatives
				if dx==0.0 then dx = FACT end
				wk1[ix], ix = dx, ix+1 -- Store dC factors for later use
			end
		end
		for k=1,nd do -- Evaluate function at data points
			for i=1,nx do p[i] = x[i][k] end
			fx[k] = newtonfc(f, p, c)
		end -- fx[k] has function values at data points
		ix, jx = 1, 1
		for i=1,nc do -- Loop over active coefficients
			if actv[i]~=0 then
				c1, dx = c[i], wk1[ix]
				c[i] = c1+dx
				pwk2 = wk2[ix] -- Pointer to work array
				for k=1,nd do -- Loop over data points
					for i=1,nx do p[i] = x[i][k] end
					pwk2[k] = (newtonfc(f, p, c) -fx[k])/dx -- Partial derivative
				end
				c[i], ix = c1, ix+1
			end
		end -- Partial derivatives stored in wk2[i][k] arrays
		for j=1,nc do -- Loop over active rows in matrix
			if actv[j]~=0 then
				row, pwk2 = rowptr[jx], wk2[jx] -- pointer to matrix & wk2 row
				for k=1,nd do -- Loop over data points
					b[jx] = b[jx] + (x[1][k] - fx[k])*pwk2[k]*(fw[k] or 1)
				end
				ix = 1
				for i=1,nc do -- Accumulate row of a[][] matrix
					if actv[i]~=0 then 
						pwk3 = wk2[ix] -- Another pointer to partial derivatives
						for k=1,nd do -- Loop over data points
							row[ix] = row[ix] + pwk3[k]*pwk2[k]*(fw[k] or 1)
						end
						ix = ix+1
					end
				end
				jx = jx+1
			end
		end -- Back for another coefficient
		if iend==0 then -- Update values except for last time through
			gauss(rowptr, b) -- Now solve for Coefficient corrections
			ix, cmax = 1, 0.0
			for i=1,nc do -- Loop over active coefficients, implement step change limitations
				if actv[i]~=0 then -- and update active coefficients
					if step[i]~=0.0 then 
						if step[i]<0.0 then -- +/- step change limits
							if abs(b[ix])>-step[i] then
								if b[ix]>0.0 then b[ix] = -step[i] 
								else b[ix] = step[i] end
							end
						else -- percentage change limits
							c1 = 1. + b[ix]/c[i]
							if c1<step[i] then b[ix] = (step[i] - 1.)*c[i] end
							if c1>1./step[i] then b[ix] = (1./step[i] - 1.)*c[i] end
						end
					end
					c1 = abs(b[ix]) -- Check on convergence
					c2 = abs(c[i]) + c1
					if c2~=0.0 then c1 = c1/c2 end
					if c1>cmax then cmax = c1 end
					c[i], ix = c[i] + b[ix], ix+1 -- Updated C's here
				end
			end
			if nprint~=0 then -- Print iteration values if desired
				printf('Coefficients at iteration %d are:\n',ilp)
				for i=1,nc do printf(' %e ',c[i]) end; printf('\n'); flush()
			end
		end
		if iend==1 then nend=ilp; break end -- End now
		if cmax <ERROR then iend = 1 end -- Converged but once more to calculate matrix
	end
	local cv = minv(rowptr,jx-1) --Covariance matrix -- needed to bypass gauss for this 
	ix, err = 1, 0.0 -- Accumulate error
	for k=1,nd do -- Sum over data points
		for i=1,nx do p[i] = x[i][k] end
		c1 = x[1][k] - newtonfc(f, p, c) -- Error evaluation 
		err = err + c1*c1 -- Square and accumulate
	end
	for j=1,nc do -- Calculate standard deviations of coefficients
		if actv[j]==0 then del[j] = 0
		else del[j],ix = sqrt(err*((type(cv)=='table' and cv[ix][ix]) or 0)/(nd-nc)),ix+1 end
	end
	return del, sqrt(err/nd), nend-1 -- return RMS error and # of iterations
end
newtonfc = function(f, x, ...) -- Newton's method for function with C coefficient array
	-- Form of f must return f(x,c) = 0 when solution is found, x is an array of values,x[1] = solved value
	local y,dy,fy,cy,ysave
	if ylinear~=0 then -- Bypass Newton iterations if linear function in variable
		y,x[1] = x[1],0.0 -- First variable is considered variable being solved for
		fy = f(x, ...) -- f must be of the form --- return cal_val - x[1]
		x[1] = y
		return fy -- Return value of function
	else
		y = x[1] -- First variable is considered the variable we are solving for, x[1] has initial guess
		ysave = y -- Save initial guess
		for i=1,NT do -- Main Newton loop
			dy = y*FACT
			if dy==0.0 then dy = FACT end -- Increment variable
			fy = f(x, ...)
			x[1] = x[1] + dy
			cy = fy - f(x, ...) -- Difference for derivative
			if cy~=0.0 then cy = fy*dy/cy end -- Correction term using numerical derivative
			y = y + cy -- New value of variable being solved for
			x[1] = y; dy = abs(cy) + abs(y)
			if dy==0.0 or abs(cy/dy)<ERR then break end -- Return if accuracy achieved
		end
		x[1] = ysave -- Restore initial guess, it might be a necessary data point
		return y -- Rerurn solved value
	end
end 
setfenv(nlstsq,{ERROR=5.e-4,FACT=1.e-8,NMAX=50,ERR=1.e-3,NT=50,
	nprint=1,ylinear=0,abs=math.abs,sqrt=math.sqrt,flush=io.flush,
	printf=printf,gauss=gauss,newtonfc=newtonfc,type=type,print=print,
	minv=minv}) -- Set parameters for nlstsq()
setfenv(newtonfc,getfenv(nlstsq)) -- Same parameters for newtonfc()



	