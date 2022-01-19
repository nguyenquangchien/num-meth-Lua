-- /* File cltwodim.lua */ -- Some basic probability and statistical functions

requires("DataFit","prob")

cltwodim = function(cx,cy,prb)
	prb = prb or 90
	local ix,c1,c2,fm,nd = 1,{},{},{},#cx
	local m1,s1,m2,s2,ul
	local ca,cb,csq = {},{},{}
	local xb,yb,NC = {},{},400
	local c1l,c1u,c2l,c2u = FH,FL,FH,FL
		f1 = datafit({cy,cx},3) -- Fit function to data
	for i=1,nd do
		c1[i],fm[i] = cx[i], f1(cx[i]); c2[i] = cy[i] - fm[i]
	end
	m1,s1 = stats(c1); m2,s2 = stats(c2)
	for i=1,nd do 
		csq[i] = ((c1[i]-m1)/s1)^2 + ((c2[i]-m2)/s2)^2  
	end
	_,ul = climits(csq,(2*prb-100))
	for i=1,nd do 
		if csq[i]>ul then 
			ca[ix],cb[ix] = cx[i],cy[i]; ix = ix+1
		end
	end
	ul = math.sqrt(ul)
	for i=1,NC+1 do
		xb[i] = ul*s1*math.cos(2*math.pi*(i-1)/NC) + m1
		yb[i] = ul*s2*math.sin(2*math.pi*(i-1)/NC) + m2 + f1(xb[i])
	end
	for i=1,NC+1 do
		if xb[i]<c1l then c1l = xb[i] end
		if xb[i]>c1u then c1u = xb[i] end
		if yb[i]<c2l then c2l = yb[i] end
		if yb[i]>c2u then c2u = yb[i] end
	end
	return {{c1l,c1u},{c2l,c2u}},{xb,yb},{ca,cb}
end
setfenv(cltwodim,{FH=math.huge,FL=-math.huge,math=math,datafit=datafit,
	climits=climits,stats=stats})




