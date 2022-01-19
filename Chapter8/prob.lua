-- /* File prob.lua */ -- Some basic probability and statistical functions
require"elemfunc" -- needed for P(x) in histnorm
require"newton" -- needed in iAstud()

local sort_mt = {__type = "sorted"} -- Type for sorted table
local Pn,Qn = elemfunc.Pn,elemfunc.Qn

hist = function(y,n,min,max) -- Converts table into histogram of values
	local nyv = #y
	n = n or 10 -- 10 default intervals
	ymin,ymax = minmax(y,nyv)
	max,min = max or ymax, min or ymin
	if max==min then max,min=ymax,ymin end
	local dy,yh,ny = (max-min)/n, {}, {}
	if (ymax>max) or (ymin<min) then -- Warn about missing data in histogram
		print("Not all data points are included in histogram table.") end
	for i=1,n do
		yh[i],ny[i] = min + (i-.5)*dy, 0
	end
	for i=1,nyv do -- Ignores points outside specified ranges
		local ix = math.floor((y[i] - min +dy)/dy)
		if ix>=1 and ix<=n then ny[ix] = ny[ix] + 1 end
	end
	setmetatable(yh,sort_mt)
	return yh,ny
end

histnorm = function(y,n,min,max) -- Return normal distribution histogram corresponding
	local nyv = #y -- to the data set y[]
	n = n or 10 -- 10 default intervals
	ymin,ymax = minmax(y,nyv)
	max,min = max or ymax, min or ymin
	if max==min then max,min=ymax,ymin end
	local dy,yh,ny = (max-min)/n, {}, {}
	local mean,std = stats(y) -- mean, std for data set
	for i=1,n do -- Evaluate histogram for normal distribution
		yh[i] = min + (i-.5)*dy
		ny[i] = nyv*(Pn((min+i*dy-mean)/std) - Pn((min+(i-1)*dy-mean)/std))
	end
	setmetatable(yh,sort_mt)
	return yh,ny
end	

stats = function(yh,nh) -- Mean, variance, std deviation, skew
	local ny = #yh
	local mean,var,stdev,skew,nt=0,0,0,0,0
	if nh==nil then 
		nh = {}; for i=1,ny do nh[i] = 1 end
	end
	for i=1,ny do nt = nt + nh[i] end
	for i=1,ny do
		mean,var = mean + nh[i]*yh[i], var + nh[i]*yh[i]^2
	end
	mean = mean/nt; var = (var - nt*mean^2)/(nt-1)
	stdev = math.sqrt(var)
	for i=1,ny do skew = skew + nh[i]*(yh[i] - mean)^3 end
	skew = skew/(nt-1)
	return mean,stdev,var,skew
end

makeODF = function(x,nx) -- Make order distribution function 
	local xx,y = {}, {[0]=0} -- Set y[0]=0
	local nd = #x
	nx = nx or {}
	for i=1,nd do xx[i] = x[i] end
	if type(x)~="sorted"  then table.sort(xx) end
	for i=1,nd do y[i] = y[i-1] + (nx[i] or 1) end
	return setmetatable(xx,sort_mt),y -- Return sorted array and cumulative distribution table
end

correlation = function(x,y) -- Linear correlation function
	local xm,ym,xv,yv,cor = 0,0,0,0,0
	local n = #x
	for i=1,n do
		xm,ym = xm + x[i], ym + y[i]
	end
	xm,ym = xm/n,ym/n
	for i=1,n do
		xv,yv = xv + (x[i]-xm)^2, yv + (y[i]-ym)^2
		cor = cor + (x[i]-xm)*(y[i]-ym)
	end
	return cor/math.sqrt(xv*yv) -- Correlation coefficient
end

minmax = function(y,ny) -- Find mnimum and maximum values in table
	local ny,ymin = ny or #y,y[1]
	local ymax,imin,imax = ymin, 1, 1
	for i=1,ny do 
		if y[i]<ymin then ymin,imin=y[i],i end
		if y[i]>ymax then ymax,imax=y[i],i end
	end
	return ymin,ymax,imin,imax
end

median = function(x) -- Median from a sorted table
	local nx,xx = #x
	if type(x)=="sorted" then xx = x
	else -- sort if not already sorted
		xx = {}; for i=1,#x do xx[i] = x[i] end
		table.sort(xx); setmetatable(xx,sort_mt)
	end
	local n2 = math.floor(nx/2)
	if n2*2==nx then return (xx[n2]+xx[n2+1])/2, xx
	else return xx[n2+1], xx end
end

clinear = function(x,y) -- Linear least squares coefficients, c[1] + c[2]*x
	local ysq,xsq,xsum,ysum,xy = 0,0,0,0,0
	local nx = #x
	for i=1,nx do
		xsq,ysq = xsq + x[i]*x[i], ysq + y[i]*y[i]
		xsum,ysum = xsum + x[i], ysum + y[i]
		xy = xy + x[i]*y[i]
	end
	local den = nx*xsq - xsum^2
	return (nx*xy-xsum*ysum)/den, (ysum*xsq-xsum*xy)/den,
		(nx*xy-xsum*ysum)^2/((nx*xsq-xsum^2)*(nx*ysq-ysum^2))
end

datasort = function(x) -- Sort a table if not already sorted
	table.sort(x)
	return setmetatable(x,{__type="sorted"})
end

climits = function(carr,val)
	val = (val or 95)/100
	local car,nv
	if type(carr[1])=='table' then 
		local nc,cret = #carr,{}
		for i=1,nc do
			car,nv = {},#carr[i]
			for j=1,nv do car[i] = carr[i][j] end
			table.sort(car[i])
			cret[i] = {car[i][math.floor(nv*(1-val)/2+1)],car[i][math.ceil(nv*(1+val)/2)]}
		end
		return cret
	else
		nv = #carr
		if type(carr)~="sorted"  then 
			car = {}; for j=1,nv do car[j] = carr[j] end
			table.sort(car)
		else car = carr end
		return car[math.floor(nv*(1-val)/2+1)],car[math.ceil(nv*(1+val)/2)]
	end
end

makeCDF = function(x,mtd) -- Make distribution function from data set
	local xx,yy = makeODF(x)
	local nx = #xx
	if mtd=='mean' then for i=1,nx do yy[i] = i/(nx+1) end -- mean rank
	elseif mtd=='sym' then for i=1,nx do yy[i] = (i-.5)/nx end -- symmetrical
	else for i=1,nx do yy[i] = (i-.3)/(nx+.4) end end -- median rank
	return xx,yy
end

normalCDF = function(x,max) -- Return normal CDF with same mean, var as input data
	local mean,std = stats(x)
	local xx,yy = makeODF(x)
	max = max or 1 -- Normally max = 1, other values for order distribution function
	for i=1,#xx do yy[i] = max*Pn((xx[i]-mean)/std) end
	return xx,yy
end

makeCDFG = function(x)
	local x1,y1,y2
	x1,y1 = makeCDF(x)
	x1,y2 = normalCDF(x)
	return x1,y1,y2
end

Punif = function(x) -- Uniform distribution
	if x<=0 then return 0
	elseif x>=1 then return 1
	else return x end
end

Pchisq = function(x,a)
	if x<=0 then return 0
	else	
		if a>300 then return 1-Qn(math.sqrt(2*x)-math.sqrt(2*a+1))
		else return elemfunc.Pigamma(x/2,a/2) end
	end
end

Qchisq = function(x,a)
	if x<=0 then return 1
	else return Qigamma(a/2,x/2) end
end

tchisq = function(n,prob)
	local al1,t1 = (1-prob)/2
	if n>2 then t1 = n-2 else t1 = 0 end
	local icf = function(x)
		return Pchisq(x,n) - al1
	end
	local v1 = newton(icf,t1)
	al1 = 1-al1
	return v1,newton(icf,t1)
end

Pcauchy = function(x)
	return 0.5 + math.atan(x)/math.pi
end

Atstud = function(t,v) -- Student's A function
	if v>300 then return 2*Pn(t*(1-0.25/v)/math.sqrt(1+0.5*t^2/v)) -1
	else	return 1 - elemfunc.ibeta(v/(v+t^2),v/2,0.5) end
end

ttable = function(n,prob)
	local Atf = function(x)
		return Atstud(x,n) - prob
	end
	return newton(Atf,0)
end

Ptdist = function(x,v) -- Students t distribution
	local val = 0.5*ibeta(v/(v+x^2),v/2,0.5)
	if x>=0 then return 1 -  val 
	else return val end
end

Pfdist = function(x,v1,v2) -- F distribution P function
	if x<=0 then return 0
	else	return 1 - ibeta(v2/(v2+v1*x),v2/2,v1/2) end
end

Qfdist = function(x,v1,v2) -- F distribution Q function
	return ibeta(v2/(v2+v1*x),v2/2,v1/2)
end

Pbeta = function(x,p,q) -- Beta distruibution function
	if x<=0 then return 0
	elseif x>=1 then return 1
	else	return ibeta(x,p,q) end
end

Pweibull = function(x,v) -- Weibull distribution function
	if x<=0 then return 0
	else	return 1 - math.exp(-x^v) end
end

weibulldata = function(dat) -- Transform rank table into Weibull plot data
	local xx,yy = dist(dat) --Generate distribution function
	local nx = #yy
	for i=1,nx do -- Uses median rank for yy[i]
		yy[i] = math.log(-math.log(1-(i-.3)/(nx+.4))) -- Y axis values
		xx[i] = math.log(xx[i]) -- X axis values
	end
	return xx,yy -- return Weibull data ready for plotting
end

Plognormal = function(x,s)
	if x<=0 then return 0
	else	return Pn(math.log(x)/s) end
end

iCDF = function(CDF,prob,xinit,a,b,c) -- General inverse CDF function
	local iCDFx =  function(x) -- Function for Newton's method
		return CDF(x,a,b,c) - prob
	end
	return newton(iCDFx,xinit) -- Calculate inverse value
end

Qks = function(x)
	local EPS1,EPS2=.001,1.e-8
	local fac,sum,termbf,a2,term=2,0,0,-2*x*x
	if x<.2 then return 1 end
	for j=1,100 do
		term = fac*math.exp(a2*j*j)
		sum = sum + term
		if math.abs(term)<=EPS1*termbf or math.abs(term)<EPS2*sum then return sum end
		fac = -fac
		termbf = math.abs(term)
	end
	print('Accuracy not achieved in Qks')
	return 0
end

KS_test = function(data1,data2,...) -- Includes possible arguments to function
	local nd1,nd2,imax
	local x1,x2,d1,d2
	local D,dt,dt1,dt2=0,0,0,0
	x1,d1 = makeCDF(data1) -- Make distribuition function
	if type(data2)=='function' then -- Comparison of data to functio()
		nd = #data1
		for i=1,nd do -- Step over data set
			dt = math.abs(d1[i] - data2(x1[i],unpack(arg)))
			if dt>D then D,imax = dt,i end
		end
	else -- Comparison of two distributions in table form
		nd1,nd2 = #data1,#data2
		nd = nd1*nd2/(nd1+nd2)
		x2,d2 = makeCDF(data2)
		local j1,j2,xv1,xv2 = 1,1
		while j1<=nd1 and j2<=nd2 do -- Step over both data sets
			xv1,xv2 = x1[j1],x2[j2]
			if xv1<=xv2 then dt1,j1 = d1[j1],j1+1 end -- next step in data1
			if xv2<=xv1 then dt2,j2 = d2[j2],j2+1 end -- next step in data2
			dt = math.abs(dt1-dt2)
			if dt>D then D,imax = dt,j1 end
		end
	end
	local kstest = D*math.sqrt(nd)
	print('\n\t\t Kolmogorov-Smirnov Goodness-of-Fit Test\n')
	if nd==nd1 then 
		print('Null Hypothesis HO:		Distribution Fits the Data')
		print('Alternate Hypothesis HA:	Distribution does not Fit the Data')
	else 
		print('Null Hypothesis HO:		Two Distributions are the Same')
		print('Alternate Hypothesis HA:	Two Distribution are not the Same')
	end
	print('Number of observations = ',engr(nd))
	print('\nK-S test statistic: D, sqrt(N)*D = '..engr(D)..engr(kstest)..'(Test value)')
	print('\nConfidence Level\tCutoff\t\tConclusion')
	prob = {.1,.05,.01,.001}
	for j=1,#prob do
		lv = iCDF(Qks,prob[j],1)
		if lv>kstest then print('',engr((1-prob[j])*100,'%'),'',engr(lv),'\t Accept HO')
		else print('',engr((1-prob[j])*100,'%'),'',engr(lv),'\t Reject HO') end
	end
	print('Accept Above Confidence Level of ',engr((1-Qks(kstest))*100,'%\n'))
end	

ran0 = function() -- Very simple random number generator
	jran = mod(jran*ia+ic,im)
	return jran/im
end
setfenv(ran0,{jran=885937,im=1771875,ia=2416,ic=374441,mod=math.fmod})

ran1 = function(mseed) -- Better random number generator
	mseed = mseed or minit
	local temp,j
	if mseed then 
		minit = nil
		while mseed>IC1 do mseed = mseed/2 end
		ix1 = mod(IC1-mseed,M1); ix1 = mod(IA1*ix1+IC1,M1)
		ix2 = mod(ix1,M2); ix1 = mod(IA1*ix1+IC1,M1)
		ix3 = mod(ix1,M3)
		for j=1,97 do 
			ix1,ix2 = mod(IA1*ix1+IC1,M1), mod(IA2*ix2+IC2,M2)
			r[j] = (ix1+ix2*RM2)*RM1
		end
	end 
	ix1,ix2,ix3 = mod(IA1*ix1+IC1,M1), mod(IA2*ix2+IC2,M2), mod(IA3*ix3+IC3,M3)
	j = floor(1 + (97*ix3/M3))
	temp,r[j] = r[j], (ix1+ix2*RM2)*RM1
	return temp
end 
setfenv(ran1,{M1=259200,IA1=7141,IC1=54773,M2=134456,IA2=8121,IC2=28411,
		M3=243000,IA3=4561,IC3=51349,RM1=1/259200,RM2=1/134456,ix1,ix2,ix3,
		r={},minit=0,mod=math.fmod,floor=math.floor})

ran2 = function(mseed) --Another random number generator
	local j,temp
	mseed = mseed or minit
	if mseed then
		minit = nil
		while mseed>IC do mseed = mseed/2 end
		ix = mod(IC-mseed,M)
		for j=1,97 do 
			ix = mod(IA*ix+IC,M); r[j] = ix
		end
		ix,iy =  mod(IA*ix+IC,M), ix
	end
	j = floor(1+97*iy/M); iy = r[j]
	ix = mod(IA*ix+IC,M); r[j] = ix
	return iy/M
end
setfenv(ran2,{M=714025,IA=1366,IC=150889,ix,iy,r={},minit=0,mod=math.fmod,
	floor=math.floor})

ran3 = function(mseed) -- Third random number generator
	local mj
	mseed = mseed or minit
	if mseed then 
		local ii,mk
		while mseed>MSEED do mseed = mseed/2 end
		mj = MSEED - abs(mseed)
		mj = mod(mj,MBIG)
		ma[55],mk = mj,1
		for i=1,54 do
			ii = mod(21*i,55)
			ma[ii] = mk
			mk = mj - mk
			if mk<MZ then mk = mk + MBIG end
			mj = ma[ii]
		end
		for k=1,4 do
			for i=1,55 do
				ma[i] = ma[i] - ma[1+mod(i+30,55)]
				if ma[i]<MZ then ma[i] = ma[i] +MBIG end
			end
		end
		inext,inextp = 0,31
		minit = nil
	end
	inext,inextp = inext+1,inextp+1
	if inext==56 then inext=1 end
	if inextp==56 then inextp=1 end
	mj = ma[inext] - ma[inextp]
	if mj<MZ then mj = mj + MBIG end
	ma[inext] = mj
	return mj*FAC
end 
setfenv(ran3,{MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1/1000000000,inext,inextp,
	minit=161803398,ma={},mod=math.fmod,abs=math.abs})

gnormal = function() -- Normally distributed random variable, unity variance
	local fac,r,v1,v2
	if iset==0 then 
		repeat
			v1,v2 = 2.*random()-1, 2.*random()-1
			r = v1*v1 + v2*v2
		until (r<1)
		fac = sqrt(-2.*log(r)/r)
		gset,iset = v1*fac,1
		return v2*fac
	else iset = 0; return gset end
end
setfenv(gnormal,{random=math.random,sqrt=math.sqrt,log=math.log,
	iset=0,gset})

lag = function(x,lag) -- Set up Lag data tables -- Default of one lag
	lag = lag or 1
	local xx,yy = {},{}
	for i=1,#x-lag do xx[i],yy[i] = x[i],x[i+lag] end
	return xx,yy
end

autocorr = function(x,nc) -- Autocorrelation function for nc lags in data
	local xa,autc,xx,yy = {},{}
	local n = #x
	nc = nc or math.min(100,n)
	local mean,sum,sum1 = 0,0,0
	for i=1,n do mean,sum1 = mean + x[i], sum1+x[i]^2 end
	mean = mean/n
	autc[1],xa[1] = sum1-n*mean^2,0
	for i=1,nc do
		xx,yy = lag(x,i)
		sum,sum1,n = 0,0,#xx
		for i=1,n do 
			sum = sum + xx[i]*yy[i] 
			sum1 = sum1 + xx[i] + yy[i]
		end
		xa[i+1] = i
		autc[i+1] = sum - mean*sum1 + n*mean^2
	end
	sum = autc[1]
	for i=1,nc+1 do autc[i] = autc[i]/sum end
	return xa,autc
end

