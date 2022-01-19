--/* File mcpar.lua */

require"prob"; require"nlstsq"

MCpar = function(yx, fw, ft, c, actv, step, nmx) 
	local cvar,cc = {},{}
	local xd = yx[2]
	local nc,nd,y,yv = #c,#xd,{},{}
	local k,kp,nnmx,nnm = 1,1,0
	local NTM = getfenv(nlstsq).NMAX-1 -- Max Newton iterations
	for i=1,nc do cvar[i] = {} end
	for i=1,nd do 
		y[i] = newtonfc(ft,{yx[1][i],xd[i]},c) 
		yv[i] = y[i] - yx[1][i]
	end
	local _,sdev = stats(yv)
	getfenv(nlstsq).nprint=0 -- Don't print iteration values
	nmx = nmx or NMAX -- Default of 1000 iterations
	while k<=nmx do -- Now MC loop
		for i=1,nc do cc[i] = c[i] end -- Initial guess at C's
		for i=1,nd do yv[i] = y[i] + sdev*gnormal() end -- New generated data set
		_,_,nnm = nlstsq({yv,xd},fw,ft,cc,actv,step) -- Analysis for C's
		if nnm>=NTM then 
			print('Convergence not achieved for MC #',k); io.flush() 
		else 
			for i=1,nc do cvar[i][k] = cc[i] end -- Save C values
			k = k+1
		end
		if nnm>nnmx then nnmx = nnm end
		if nprint~=nil then
			if kp==NP then print('Completed MC simulation #',k-1)
			io.flush(); kp = 0 end
		end
		kp = kp+1
	end
	return cvar, nnmx -- Return array of C values as {{c[1]'s},{c[2]'s},...{c[n]'s}}
end
setfenv(MCpar,{NMAX=1000,NP=100,fstp=.2,stats=stats,print=print,io=io,newtonfc=newtonfc,
	getfenv=getfenv,nlstsq=nlstsq,gnormal=gnormal,nprint=nil})

