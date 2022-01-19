-- /* File elemfunc.lua */ -- Some definitions of elementary functions
-- Probability function, Error function, Complementary Error function, Gamma function

elemfunc = {} -- Table for elementary functions
local pp = function(x) -- Function used by Pn,Qn
	return 0.5*(((((((.0000053830*x + .0000488906)*x + .0000380036)*x +
		.0032776263)*x + .0211410061)*x + .0498673470)*x + 1)^-16)
end
elemfunc.Pn = function(x) -- Probability function (~7 digits of accuracy)
	if x>=0 then return 1.0 - pp(x)
	else return pp(-x) end
end
elemfunc.Qn = function(x) -- 1-P(x) (~7 digits of accuracy)
	if x>=0 then return pp(x)
	else return 1.0-pp(-x) end
end

elemfunc.erf = function(x) -- Error function (~7 digits of accuracy)
	local sgn = 1; if x<0 then x,sgn = math.abs(x), -1 end
	return sgn*(1-((((((.0000430638*x + .0002765672)*x + .0001520143)*x + 
		.0092705272)*x + .0422820123)*x + .0705230784)*x + 1)^-16)
end
elemfunc.erfc = function(x) return 1-elemfunc.erf(x) end -- Complementary Error function

elemfunc.gamma = function(x) -- Gamma function (x-1)! (~7 digits of accuracy)
	local ans,sign,fac=1,1,1
	if x<0 then sign,fac,x = -1, math.pi/(x*math.sin(-math.pi*x)), -x end
	if x<1 then ans,x = ans/x, x+1 end
	x = x-1
	while x>1 do ans,x = ans*x, x-1 end
	ans = ans*((((((((.035868343*x - .193527818)*x + .482199394)*x - .756704078)*x +
		.918206857)*x - .897056937)*x + .988205891)*x - .577191652)*x + 1)
	if sign>0 then return ans*fac
	else return fac/ans end
end

local gser = function(x,a) -- Series expansion for Incomplete Gamma function
	local sum,del,ap,n
	local gln = log(gamma(a))
	ap,del,sum = a,1/a,1/a
	for i=1,ITMAX do
		n = i
		ap = ap + 1
		del = del*x/ap
		sum = sum + del
		if abs(del) < abs(sum)*EPS then break end
	end
	if n==ITMAX then print('ITMAX exceeded in gser,Pigmma ',ITMAX ) end
	return sum*exp(-x + a*log(x) - gln),gln 
end
setfenv(gser,{abs=math.abs,exp=math.exp,log=math.log,gamma = elemfunc.gamma,
	ITMAX=200,	EPS=3.e-7,print=print})

local gcf = function(x,a) -- Continued fraction for Incomplete Gamma function
	local gold,fac,b1,g,n = 0,1,1
	local b0,a0,a1,anf,ana,an,ai = 0,1
	local gln = log(gamma(a))
	a1 = x
	for i=1,ITMAX do
		ana, n = i - a, i
		a0, b0, anf = (a1+a0*ana)*fac, (b1+b0*ana)*fac, i*fac
		a1, b1 = x*a0 + anf*a1, x*b0 + anf*b1
		if a1 then 
			fac = 1/a1
			g = b1*fac
			if abs((g-gold)/g)<EPS then break end
			gold = g
		end
	end
	if n==ITMAX then print('ITMAX exceeded in gcf,Pigmma ',ITMAX) end
	return 1-exp(-x + a*log(x) - gln)*g,gln
end
setfenv(gcf,getfenv(gser))

elemfunc.Pigamma = function(x,a) -- Incomplete P Gamma function P(a,x)
	if (x<0.0 or a<=0.0) then return 0 end --print('Invalid arguments in Pingam'); return 0 end
	if x<(a+1) then return gser(x,a) 
	else return gcf(x,a) end
end

elemfunc.Qigamma = function(a,x) -- 1 - Incomplete P Gamma function
	local x,g = elemfunc.Pigamma(a,x)
	return 1-x,g
end

elemfunc.igamma = function(a,x) -- Incomplete lower case gamma function g(a,x)
	local x,g = elemfunc.Pigamma(a,x)
	return x*math.exp(g)
end

elemfunc.beta = function(z,w) --Beta function
	return elemfunc.gamma(z)*elemfunc.gamma(w)/elemfunc.gamma(z+w)
end

local betacf = function(x,a,b) -- Continued fraction expansion for ibeta
	local qap,qam,qab,em,tem,d
	local bm,bz,bp,bpp = 1
	local az,am,ap,app,aold = 1,1
	qab,qap,qam = a+b,a+1,a-1
	bz = 1 - qab*x/qap
	for i=1,ITMAX do
		em,tem = i,i+i
		d = em*(b-em)*x/((qam+tem)*(a+tem))
		ap,bp = az + d*am, bz + d*bm
		d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
		app,bpp,aold = ap + d*az, bp+d*bz, az
		am,bm,az,bz = ap/bpp, bp/bpp, app/bpp, 1
		if abs(az-aold)<EPS*abs(az) then break end
	end
	if em==ITMAX then print('ITMAX exceeded in betacf,ibeta') end
	return az
end
setfenv(betacf,getfenv(gser))


elemfunc.ibeta = function(x,a,b) -- Incomplete beta function
	local bt
	if x<0 then return 0 end
	if x>1 then return 1 end
	--if x<0 or x>1 then print('Bad argument to ibeta') end
	if x<=0 or x>=1 then bt = 0
	else bt = elemfumc.gamma(a+b)*x^a*(1-x)^b/(elemfunc.gamma(a)*elemfunc.gamma(b)) end
	if x<(a+1)/(a+b+2) then return bt*betacf(x,a,b)/a
	else return 1-bt*betacf(1-x,b,a)/b end
end

elemfunc.E1 = function(x) -- Exponential integral (~7 digits of accuracy)
	if x<1 then return ((((((.00107857*x-.00976004)*x+.05519968)*x-.24991055)*x+
		.99999193)*x-.57721566) - math.log(x))
	else return (((((x+8.5733287401)*x+18.059016973)*x+8.6347608925)*x+.2677737343)/
		((((x+9.5733223454)*x+25.6329561486)*x+21.0996530827)*x+3.9584969228))*
		math.exp(-x)/x end
end	
elemfunc.En = function(x,n) -- General Exponential integral of order n (~7 digits of accuracy)
	local ans,n1 = E1(x),1
	if x==0 and n>1 then ans = 0 end
	while n1<n do
		ans = (math.exp(-x) - x*ans)/n1
		n1 = n1+1
	end
	return ans
end

local f = function(x) -- Function used by Si, Ci
	local x2 = x*x
	return ((((x2+38.027264)*x2+265.187033)*x2+335.67732)*x2+38.102495)/
		((((x2+40.021433)*x2+322.624911)*x2+570.23628)*x2+157.105423)/x
end
local g = function(x) -- Function used by Si, Ci
	local x2 = x*x
	return ((((x2+42.242855)*x2+302.757865)*x2+352.018498)*x2+21.821899)/
		((((x2+48.196927)*x2+482.485984)*x2+1114.978885)*x2+449.690326)/x2
end
elemfunc.Si = function(x) -- Sine integral (~7 digits of accuracy)
	x = math.abs(x)
	if x<1 then return ((((1.5371909046e-3*x+1.6492137837e-4)*x-
		5.5629931428e-2)*x+1.1051345311e-5)*x+1)*x 
	else return math.pi/2 - f(x)*math.cos(x) - g(x)*math.sin(x) end
end
elemfunc.Ci = function(x) -- Cosine integral (~7 digits of accuracy)
	x = math.abs(x)
	if x<1 then return 0.5772157 + math.log(x) - elemfunc.Cin(x)
	else return f(x)*math.sin(x) - g(x)*math.cos(x) end
end
elemfunc.Cin = function(x) -- Alternative Cosine integral
	x = math.abs(x)
	if x<1 then return ((((6.5335276944e-4*x-1.1119178653e-2)*x+
		3.4884526934e-4 )*x+2.4992271202e-1)*x+5.7923592050e-6)*x
	else return math.log(x) + 0.5772157 - elemfunc.Ci(x) end
end

elemfunc.K = function(mx) -- Complete Eliptic Integral of First Kind (~8 digits of accuracy)
	local m=1-mx
	return ((((.01451196212*m+.03742563713)*m+.03590092383)*m+
	.09666344259)*m+1.38629436112) + ((((.00441787012*m+.03328355346)*
	m+.06880248576)*m+.12498593597)*m+.5)*math.log(1/m)
end
elemfunc.E = function(mx) -- Complete Eliptic Integral of Second Kind (~8 digits of accuracy)
	local m=1-mx
	if m<1.e-14 then return 1 end
	return ((((.01736506451*m+.04757383546)*m+.0626060122)*m+
	.44325141463)*m+1) + ((((.00526449639*m+.04069697526)*
	m+.09200180037)*m+.24998368310)*m)*math.log(1/m)
end
local JYfth = function(x) -- Function used by J0,Y0
	local x3,fo,th = 3./x
	fo = (((((.00014476*x3-.00072805)*x3+.00137237)*x3-.00009512)*x3-
	.0055274)*x3-.00000077)*x3+.79788456
	th = (((((.00013558*x3-.00029333)*x3-.00054125)*x3+.00262573)*x3-
	.00003954)*x3-.04166397)*x3-.78539816+x
	return fo,th
end
elemfunc.J0 = function(x) -- Bessel function J0 (~8 digits of accuracy)
	local x3 = math.abs(x/3); x = math.abs(x)
	if x3<1 then
		x3 = x3*x3
		return (((((.00021*x3-.00394)*x3+.0444478)*x3-.3163866)*x3+
		1.2656208)*x3-2.2499997)*x3+1.0
	else 
		local fo,th = JYfth(x)
		return fo*math.cos(th)/math.sqrt(x)
	end
end
elemfunc.Y0 = function(x) -- Bessel function Y0 (~8 digits of accuracy)
	local x3 = math.abs(x/3); x = math.abs(x)
	if x3<1 then 
		x3 = x3*x3
		return (((((-.00024846*x3+.00427916)*x3-.04261214)*x3+
		.25300117)*x3-.74350384)*x3+.60559366)*x3+.36746691+
		(2/math.pi)*math.log(.5*x)*elemfunc.J0(x)
	else
		local fo,th = JYfth(x)
		return fo*math.sin(th)/math.sqrt(x)
	end
end
local JY1fth = function(x) -- Function used by J1,Y1
	local x3,fo,th = 3./x
	fo = (((((-.00020033*x3+.00113653)*x3-.00249511)*x3+.00017105)*x3+
	.01659667)*x3+.00000156)*x3+.79788456
	th = (((((-.00029166*x3+.00079824)*x3+.00074348)*x3-.00637879)*x3+
	.0000565)*x3+.12499612)*x3-2.35619449+x
	return fo,th
end
elemfunc.J1 = function(x) -- Bessel function J1 (~8 digits of accuracy)
	local x3,sgn = math.abs(x/3),1
	if x3<1 then -- x < 3
		x3 = x3*x3
		return x*((((((.00001109*x3-.00031761)*x3+.00443319)*x3-.03954289)*x3+
		.21093573)*x3-.56249985)*x3+.5)
	else -- x >= 3
		if x<0 then x,sgn = math.abs(x),-1 end
		local fo,th = JY1fth(x)
		return sgn*fo*math.cos(th)/math.sqrt(x)
	end
end
elemfunc.Y1 = function(x) -- Bessel function Y1 (~8 digits of accuracy)
	local x3 = math.abs(x/3); x = math.abs(x)
	if x3<1 then -- x < 3
		x3 = x3*x3
		return ((((((.0027873*x3-.0400976)*x3+.3123951)*x3-
		1.3164827)*x3+2.1682709)*x3+.2212091)*x3-.6366198)/x+
		(2/math.pi)*math.log(.5*x)*elemfunc.J1(x)
	else -- x >= 3
		local fo,th = JY1fth(x)
		return fo*math.sin(th)/math.sqrt(math.abs(x))
	end
end
elemfunc.Jn = function(x,n) -- Nth order Jn Bessel function,
	if n==0 then return elemfunc.J0(x)
	elseif n==1 then return elemfunc.J1(x)
	elseif x==0.0 then  return 0.0
	elseif x>n/2 then -- Forward recurrence
		local j0,j1 = elemfunc.J0(x),elemfunc.J1(x)
		for k=3,n+1 do j0,j1 = j1, 2*(k-2)*j1/x-j0 end
		return j1
	else -- Reverce recurrence
		local ja,j0,n2={},0,math.max(n+n,10); ja[n2],ja[n2-1] = 0,1
		for k=n2-2,0,-1 do ja[k] = 2*(k+1)*ja[k+1]/x-ja[k+2] end
		j0 = ja[0]; for k=2,n2,2 do j0 = j0+2*ja[k] end
		return ja[n]/j0
	end
end
elemfunc.Yn = function(x,n) -- Nth order Yn Bessel function
	local y0,y1 = elemfunc.Y0(x), elemfunc.Y1(x)
	for k=3,n+1 do y0,y1 = y1, 2*(k-2)*y1/x-y0 end -- Recursion
	return y1
end
