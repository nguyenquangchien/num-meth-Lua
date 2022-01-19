-- /* File list8_2.lua */ Test of random number generators

require"prob"

rand = ran0 -- or math.random, ran1, ran2, ran3 or gnormal
nmb = 400 -- or 1000 or 10000

x = {}; for i=1,nmb do x[i] = rand() end
print('mean,std,variance,skew = ',stats(x))

xx,yy = makeODF(x)
write_data('list8_2a.dat',xx,yy); plot(xx,yy)

xx,yy = lag(x)
write_data('list8_2b.dat',xx,yy); scatterplot(xx,yy)

xx,yy = makeCDF(x)
write_data('list8_2c.dat',xx,yy); plot(xx,yy)

