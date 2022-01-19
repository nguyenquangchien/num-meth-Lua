--/* File list7_2.lua */  Example of solar data for linear data fit
require'prob'

jk,yyd,ssd,_ = {},{},{},{}
read_data('monthssn.dat',jk,yyd,ssd) -- read data
nmths = #jk -- Number of months

prd = 60 -- Approximate period in months
pn,pks,mns,ssp = {},{},{},{}
i,k,ixx = 1,1,1
while i<nmths do
	mx = 0
	for j=ixx,prd+ixx do -- search for peak
		if i>=nmths then break end
		if ssd[i]>mx then mx,mtn,ixm = ssd[i],yyd[i],i end
		i = i+1
	end -- Now have peak
	i = ixm; prd = 120
	pks[k] = mtn; ssp[k] = mx
	mx,ixx = 1.e10,ixx+ixm
	for j=ixx,prd+ixx do -- search for minimum
		if i>=nmths then break end
		if ssd[i]<mx then mx,mtn,ixm = ssd[i],yyd[i],i end
		i = i+1
	end -- Now have minimum
	pn[k],mns[k],k = k,mtn,k+1
	ixx = ixx+ixm
	if i>=nmths then break end
	i = ixm
end
print('Number peaks and minima =',#pks,#mns)
a1,a2 = clinear(pn,pks), clinear(pn,mns)
print('solar cycle in years =',a1,a2)
write_data('list7_2.dat',pks,mns,ssp)

