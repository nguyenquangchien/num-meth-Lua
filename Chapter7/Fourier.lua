-- /* Fourier Series Analysis */

require"Complex"

fourier = function(fk,nc,thmin) -- Evaluate Fourier Cos-Sin  and Cos-Ang Components for time points fk
	nc = nc or 20
	thmin = thmin or 0
	local nd,fcs,fca,fexp = #fk,{},{},{} 
	local pi2,sum,th = 2*math.pi/nd
	for i=0,nc do
		fac,th = pi2*i,i*thmin
		sum1,sum2 = 0,0
		for k=1,nd do
			sum1 = sum1 + fk[k]*math.cos(fac*(k-1)+th)
			sum2 = sum2 + fk[k]*math.sin(fac*(k-1)+th)
		end 
		fcs[i+1] = {2*sum1/nd, 2*sum2/nd} -- Cos-Sin Components -- i=0 has DC component
	end -- Next line does Cos-Ang components
	fcs = fourier_trunc(fcs) 
	for i=1,nc+1 do fca[i] = {math.sqrt(fcs[i][1]^2+fcs[i][2]^2),math.atan2(-fcs[i][2],fcs[i][1])} end
	for i=1,nc+1 do fexp[i] = {fca[i][1]/2,fca[i][2]} end -- Exponential form in magnitude-angle form
	return setmetatable(fca,Fourierca_mt),setmetatable(fexp,Fourierexp_mt),setmetatable(fcs,Fouriercs_mt),nc
end

-- Assumed forms are f(t) = a[1] + sum(over n>1){a[n]*Cos(nwt+theta[n])} = a[1] + sum(over n>1){a[n]*Cos(nwt) + b[n]*Sin(nwt)}
-- Or for exponential form f(t) = sum(over all -/+ n){alfa[n]*Exp(jnwt)}, Only positive n values calculated; alfa(-n) = conjugate(alfa(n))
Fourier = function(f,tmin,tmax,nc) -- Evaluat nc terms of Fourier Series
	local nd,fd,t = 1024,{}
	nc,tmax,tmin = nc or 20,tmax or 1,tmin or 0
	dt = (tmax-tmin)/nd
	for i=1,nd do fd[i] = f(tmin + (i-1)*dt) end -- 1024 samples of time function
	return fourier(fd,nc,tmin/(tmax-tmin)*2*math.pi) -- Returns Cos-Ang and Cos-Sin components
	--return fca,fexp,fcs
end -- Only positive harmonics returned for Exp form 

iFourier = function(fcs,nt,nc) -- Inverse Fourier Series for time function
	nt,nc = nt or 512, nc or 1.e20 -- Number of time points, number of Fourier components
	local nm,typ,fl = #fcs, type(fcs), {}
	if nc>nm then nc = nm end -- Can't use more components than available 
	if typ=='CosSin' then for i=1,nc do fl[i] = fcs[i] end 
	elseif typ=='CosAng' then for i=1,nc do fl[i] = {fcs[i][1]*math.cos(fcs[i][2]),
			-fcs[i][1]*math.sin(fcs[i][2])} end
	elseif typ=='ExpAng' then for i=1,nc do fl[i] = {2*fcs[i][1]*math.cos(fcs[i][2]),
			-2*fcs[i][1]*math.sin(fcs[i][2])} end
	end
	local ft,tv = {},{}
	local pi2,sum,erun,efac = 2*math.pi/nt
	local fdc = fl[1][1]/2 -- a[0]/2 is DC term
	for i=1,nt do
		fac,sum = (i-1)*pi2,0
		for k=2,nc do -- Cos-Sin form is faster than Exp form with Complex numbers
			sum = sum + fl[k][1]*math.cos((k-1)*fac) + fl[k][2]*math.sin((k-1)*fac)
		end
		ft[i], tv[i] = sum + fdc, (i-1)/nt
	end
	return tv,ft -- Time, between 0 and 1 and Function values returned
end

Fouriercs_mt = {__type = 'CosSin', -- Type definition for Cos-Sin form
__tostring = function(x)
	local st = ''
	for i=1,#x do st = st..'(('..x[i][1]..')Cos('..(i-1)..') + ('..x[i][2]..')Sin('..(i-1)..'))\n' end
	return st
end }
Fourierca_mt = {__type = 'CosAng', --Type defintion for Cos-Ang form
__tostring = function(x)
	local st = ''
	for i=1,#x do st = st..'(('..x[i][1]..')Cos('..(i-1)..' + ('..x[i][2]..' Rad))\n' end
	return st
end}
Fourierexp_mt = {__type='ExpAng', -- Type defition for Exp-Ang form
__tostring = function(x)
	local st = ''
	for i=1,#x do st = st..'(('..x[i][1]..')Exp('..(i-1)..' + ('..x[i][2]..' Ang))\n' end
	return st
end}

function fourier_trunc(x) -- Eliminate small roundoff errors from calculations
	local n = #x -- x is assumed to be in Cos-Sin form
	local fact,fm = 0,0
	for i=1,n do 
		fact = math.sqrt(x[i][1]^2 + x[i][2]^2)
		if fact>fm then fm = fact end
	end
	fact = fm*1.e-12
	for i=1,n do
		if math.abs(x[i][1])<fact then x[i][1] = 0 end
		if math.abs(x[i][2])<fact then x[i][2] = 0 end
	end
	return x
end

plotcossin = function(fx) -- Plot Cos-Sin representation
	if type(fx)~='CosSin' then print('Wrong format in plotcossin.'); return end
	local na,ca,sa = {},{},{}
	local n = #fx
	for i=1,n do na[i],ca[i],sa[i] = i-1,fx[i][1],fx[i][2] end
	stem(na,ca,{'Cos Coefficients','n\'th Harmonic','Cos(n) value'})
	stem(na,sa,{'Sin Coefficients','n\'th Harmonic','Sin(n) value'})
end

plotcosang = function(fx) -- Plot Cos-Ang representation
	if type(fx)~='CosAng' then print('Wrong format in plotcosang.'); return end
	local na,ca,aa = {},{},{}
	local n = #fx
	for i=1,n do na[i],ca[i],aa[i] = i-1,fx[i][1],fx[i][2] end
	stem(na,ca,{'Cos Coefficients','n\'th Harmonic','Cos(n+theta)'})
	stem(na,aa,{'Angle theta','n\'th Harmonic','Theta in Radians'})
end

plotexpang = function(fx) -- Plot Complex Exp representation -- Mag and angle
	if type(fx)~='ExpAng' then print('Wrong format in plotexpang.'); return end
	local na,em,ea = {},{},{}
	local n = #fx
	for i=1,n do na[i],em[i],ea[i] = i-1,fx[i][1],fx[i][2] end
	stem(na,em,{'Exp Coefficient','n\'th Harmonic','Exp(j*n) value'})
	stem(na,ea,{'Angle','n\'th Harmonic','Radians'})
end

plotFourier = function(fx)
	local tp = type(fx)
	if tp=='CosSin' then return plotcossin(fx) 
	elseif tp=='CosAng' then return plotcosang(fx)
	else return plotexpang(fx) end
end
		
expandFt = function(f,tmin,tmax,np) -- Expand a periodic function and 
	np,nd = np or 3,#f -- set time limits for a period
	local ft,tv = {},{}
	local t,dt = tmin,(tmax-tmin)/(np*nd) -- Time limits on expansion
	k = 1
	for i=1,np do -- Number of periods to expand
		for j=1,nd do
			tv[k], ft[k] = t, f[j]
			t,k = t+dt,k+1
		end
	end
	return tv,ft -- Time and Function tables returned
end
	
