 --/* File list3_6.lua */ -- Diode equation example
 
require "newton_1"
 
Is,vt = 1.0e-14,0.026 -- Diode parameters
R = 2.0e3 -- Resistance value
 
function f1(vd) return Is*R*(math.exp(vd/vt)-1) + vd - vs end
 
function f2(vd) return vd - vt*math.log((vs - vd)/(Is*R) + 1) end

vs = 15 -- Positive voltage source
print("Positive voltage solutions")
print(newton(f1,0))
print(newton(f1,.4))
print(newton(f1,.8))
print(newton(f2,0))
print(newton(f2,.4))
print(newton(f2,.8))

vs = -15 -- Negative voltage source
print("Negative voltage solutions")
print(newton(f1,0))
print(newton(f1,-1))
print(newton(f1,.8))
print(newton(f2,0))
print(newton(f2,-1))
print(newton(f2,-14.999))
