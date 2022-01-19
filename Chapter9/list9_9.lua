-- /* File list9_9.lua */ -- General file for data fitting, Parameter estimation 

requires("nlstsq","mcpar","DataFit","cltwodim")

-- Section 1. Input data file and define fitting function
-- Change lines 6 through 14 for different data sets
infile = 'Rat42'
xd,yd = {},{}; read_data(infile..'.txt',yd,xd) -- y stored first
nd = #xd

ft = function(x,c) -- Define function to fit data
	return c[1]/(1 + math.exp(c[2] - c[3]*x[2])) - x[1]
end
c = {100,1,.1}; nc = #c; step ={.5,.5,.5}

-- Section 2. Perform data fit and print fitted parameters
del,err,nm = nlstsq({yd,xd},fw,ft,c,actv,step) -- Call fiting, print max iterations 
print('Max iteration number =',nm)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end

-- Section 3. Generate 6-Plot Data, Only Plots 1,2,3 and 6 are used here
ycalc,res = {},{} -- Calculate fitted values and residuals
for i=1,nd do ycalc[i] = ft({0,xd[i]},c); res[i] = yd[i] - ycalc[i] end
fres1,fres2 = datafit({res,xd}), datafit({res,yd}) -- Fit residuals
yres1,yres2,xyp,xypn = {},{},{},{}
for i=1,nd do yres1[i],yres2[i] = fres1(xd[i]),fres2((ycalc[i])) end
xyp[1],xyp[2] = makeCDF(res); xypn[1],xypn[2] = normalCDF(res) -- CDF data for residuals
plot(xd,yd,ycalc,{"Data plot"}); 
scatterplot(xd,yres2,res,{"Residual vs X","X","Res"})
scatterplot(ycalc,yres2,res,{"Residual vs Y","Y","Res"}); plot(xyp,xypn,{"CDF plot","Res","CDF"}) 
write_data(infile..'a.dat',xd,yd,ycalc,res,yres1,yres2,xyp,xypn)

-- Section 4. Generate MC Data
cvar = MCpar({yd,xd},fw,ft,c) -- 1000 MC loops -- default value
xcdf,ycdf = {},{}
for i=1,nc do xcdf[i],ycdf[i] = makeCDF(cvar[i]) end
scatterplot(cvar[1],cvar[2],{"c[2] vs c[1]","c[1]","c[2]"}); 
if nc>2 then scatterplot(cvar[1],cvar[3],{"c[3] vs c[1]","c[1]","c[3]"})
	scatterplot(cvar[2],cvar[3],{"c[3] vs c[2]","c[2]","c[3]"}) end
prob = {68.27,95.45,90,95,99,99.9}; np = #prob
print('Individual confidence limits')
for j=1,nc do
	for i=1,np do printf('Limits for c[%d] at %s = %12.4e  to  %12.4e\n',j,engr(prob[i],'%'),
		climits(xcdf[j],prob[i])) end
end

-- Section 5. Generate data for joint confidence bounds
cl,xyb,cx = cltwodim(cvar[1],cvar[2],90); cl2,xyb2,cx2 = cltwodim(cvar[1],cvar[2],95)
print('Joint confidence limits')
for i=1,2 do
	for j=1,2 do print('c['..i..'], 95% limits = ',cl2[i][j],'90% limits = ',cl[i][j]) end
end
scatterplot(cvar[1],cvar[2],cx,cx2,'with lines',xyb,xyb2,{"c[2] vs c[1] with limits","c[1]","c[2]"}) -- Joint confidence limits plot
write_data(infile..'b.dat',cvar,xyb,cx,xyb2,cx2)
