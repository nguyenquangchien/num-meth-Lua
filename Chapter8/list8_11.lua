-- /* File list8_11.lua */ -- KS test

require"prob"; require"elemfunc"
local Pn = elemfunc.Pn

paper = {}; read_data('paper.txt',paper) -- Data set

Pnms = function(x,mean,std) -- Gaussian definition for tests
	return Pn((x-mean)/std)
end
KS_test = function(data1,data2,...) -- Includes possible arguments to function
	local arg,nd1,nd2,imax = {...}
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
	print('Null Hypothesis HO:		Distribution Fits the Data')
	print('Alternate Hypothesis HA:	Distribution does not Fit the Data')
	print('Number of observations = ',engr(nd))
	print('\nK-S test statistic: D, sqrt(N)*D = '..engr(D)..','..engr(kstest)..'(Test value)')
	print('\nConfidence Level\tCutoff\t\tConclusion')
	prob = {.1,.05,.01,.001}
	for j=1,#prob do
		lv = iCDF(Qks,prob[j],1)
		if lv>kstest then print('',engr((1-prob[j])*100,'%'),'',engr(lv),'\t Accept HO')
		else print('',engr((1-prob[j])*100,'%'),'',engr(lv),'\t Reject HO') end
	end
	print('Accept Above Confidence Level of ',engr((1-Qks(kstest))*100,'%\n'))
end	

mean,std = stats(paper)
KS_test(paper,Pnms,mean,std) -- Call with comparison equation
KS_test(paper,Pnms,1.05*mean,std)
KS_test(paper,Pnms,mean,1.5*std)
