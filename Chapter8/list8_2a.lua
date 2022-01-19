-- /* File list8.1.lua */ Test of random number generators

require"prob"

rand = math.random -- or ran0, ran1, ran2, ran3 or gnormal
nmb = 40000 -- or 1000 or 10000
x = {}
for i=1,nmb do x[i] = rand() end
print('mean,std,variance,skew = ',stats(x))
xx,yy = makeODF(x)
write_data('list8_2a.dat',xx,yy)
xxd,xd = {},{}
for i=1,nmb-1 do xxd[i],xd[i] = i,xx[i+1]-xx[i] end
write_data('list8.tmp.dat',xxd,xd)
scatterplot(xxd,xd)
plot(xx,yy)
xl,yl = lag(x)
write_data('list8_2b.dat',xl,yl)
scatterplot(xl,yl)
xd,yd = makeCDF(x)
write_data('list8_2c.dat',xd,yd)
plot(xd,yd)

