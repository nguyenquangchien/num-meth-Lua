-- /* File list4_21.lua */
-- Examples of Roots of 20th order Polynomial

require "Polynomial"

p= Polynomial.new{-1,1}*{-2,1}*{-3,1}*{-4,1}*{-5,1}*{-6,1}
p = p*{-7,1}*{-8,1}*{-9,1}*{-10,1}*{-11,1}*{-12,1}*{-13,1}
p = p*{-14,1}*{-15,1,}*{-16,1}*{-17,1}*{-18,1}*{-19,1}*{-20,1}
print(p)

rts = p:roots()
for i=1,#rts do
	rtt = math.floor(rts[i]+.5)
	print(i,rtt,rts[i],(rts[i]-rtt)/rtt)
end
