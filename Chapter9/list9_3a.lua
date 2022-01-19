-- /* File list9_3a.lua */ -- MC analysis for Figure 9.1 data

requires("nlstsq","mcpar")

xd,yd = {},{} -- Data arrays
read_data('list7_1.dat',xd,yd) -- Read in data set to be fitted

ft = function(x,c) -- Define function to fit data
	return c[1]*x[2]/(x[2]^c[2]+c[3])^(1/c[2]) - x[1]
end

c,actv = {1,2.7,100},{1,0,1}; nc = #c -- Initial guess at coefficients
del,err,nmax = nlstsq({yd,xd},fw,ft,c,actv,step) -- Call fiting 
print(err,nmax); io.flush() -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
-- Now perform MC simulations
getfenv(MCpar).nprint=1 -- Follow MC development
step = {1.1,1.1,1.2} -- Need some limits for convergence
cvar,nm = MCpar({yd,xd},fw,ft,c,actv,step,1000) -- 1000 MC loops
print('Maximum Newton iterations =',nm)
-- Analyze results
xcdf,ycdf = {},{}
for i=1,nc do xcdf[i],ycdf[i] = makeCDF(cvar[i]) end
write_data('list9_3a.dat',cvar[1],cvar[2],cvar[3],xcdf[1],ycdf[1],xcdf[2],ycdf[2],xcdf[3],ycdf[3])
for i=1,nc do plot(xcdf[i],ycdf[i]) end
for i=1,nc do stem(hist(cvar[i],40)) end
for i=1,nc do print(stats(cvar[i])) end
-- Correlation plots
scatterplot(cvar[1],cvar[2]); scatterplot(cvar[1],cvar[3]); scatterplot(cvar[2],cvar[3])
prob = {68.27,95.45,90,95,99,99.9}
np = #prob
for j=1,nc do
	print('') -- spacer line
	for i=1,np do printf('Limits for c[%d] at %s = %12.4e  to  %12.4e\n',
		j,engr(prob[i],'%'),climits(xcdf[j],prob[i])) end
end

