--/* File list3_9.lua */  -- Diode equation with iteration limits

require "newton"
 
Is,vt = 1.0e-14,0.026 -- Diode parameters
R = 2.0e3 -- Resistance value

function f1(vd) return Is*R*(math.exp(vd/vt)-1) + vd - vs end

vs = 15 -- Positive source voltage
print("Positive voltage solutions")
print(newton(f1,0,-1))
print(newton(f1,.4,2))
print(newton(f1,.8,-1))

vs = -15 -- Negative voltage source
print("Negative voltage solutions")
print(newton(f1,0,-1))
print(newton(f1,-1,2))
print(newton(f1,.8,-2))
