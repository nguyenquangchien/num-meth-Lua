-- /* File list8_3.lua */ Test of gaussian random number generator

require"prob"

rand = gnormal -- Normal density function
nmb = 10000 -- or 1000 
x = {}
for i=1,nmb do x[i] = rand() end -- Generate random numbers
print('mean,std,variance,skew = ',stats(x))

xx,yy = makeCDF(x) -- Make CDF from random numbers
xn,yn = normalCDF(x) -- Make normal CDF for comparison
write_data(100,'list8_3a.dat',xx,yy,yn) -- Save every 100 data point
plot(xx,yy,yn)

xl,yl = lag(x) -- Generate lag data
write_data('list8_3b.dat',xl,yl) -- Save lag data
scatterplot(xl,yl) -- Plot lag data

xd,yd = hist(x,50,-4,4) -- Generate histogram, 50 bins from -4 to +4
xp,yp = histnorm(x,50,-4,4) -- Normal histogram for comparison
write_data('list8_3c.dat',xd,yd,yp) -- Save histograms
plot(xd,yd,yp) -- Plot histograms

