-- /* File list8_5.lua */  Analysis of data for paper thickness

require"prob"

x,y = {},{}
read_data('paper.txt',y)
nd = #y
for i=1,nd do x[i] = i end
plot(x,y)
m,std,var,skew = stats(y)
print('mean,std,variance,skew = ',m,std,var,skew)
xx,yy = makeCDF(y)
xn,yn = normalCDF(y)
plot(xx,yy,yn)
xl,yl = lag(y)
scatterplot(xl,yl)
xh,yh = hist(y,20)
xph,yph = histnorm(y,20)
plot(xh,yh,yph)
write_data('list8_5a.dat',x,y,xx,yy,yn)
write_data('list8_5b.dat',xl,yl)
write_data('list8_5c.dat',xh,yh,yph)

