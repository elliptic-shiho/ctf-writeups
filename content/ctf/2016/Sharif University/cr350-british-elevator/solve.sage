from sage.all import *

"""
p = 16857450949524777441941817393974784044780411511252189319
A = 16857450949524777441941817393974784044780411507861094535
B = 77986137112576 

E = EllipticCurve(GF(p), [A, B])

print E.order() == p
print E.order()

g = E(5732560139258194764535999929325388041568732716579308775, 14532336890195013837874850588152996214121327870156054248)
v = E(2609506039090139098835068603396546214836589143940493046, 8637771092812212464887027788957801177574860926032421582)

"""

p = 853
A = 108
B = 4

E = EllipticCurve(GF(p), [A, B])

g = E(0, 2)
v = E(563, 755)

def hensel_lift(curve, p, point):
    A, B = map(long, (E.a4(), E.a6()))
    x, y = map(long, point.xy())

    #fr = y**2 - (x**3 + A*x + B)
    fr =  ((x**3 + A*x + B) - y**2) / p # = k
    t = fr % p
    print t
    t *= inverse_mod(2 * y, p)  # (y**2)' = 2 * y
    t %= p
    print t
    new_y = y + p * t
    print new_y
    return x, new_y

# lift points
x1, y1 = hensel_lift(E, p, g)
x2, y2 = hensel_lift(E, p, v)

print "P :",x1, y1
print "Q :",x2, y2

# calculate new A, B (actually, they will be the same here)
mod = p ** 2

A2 = y2**2 - y1**2 - (x2**3 - x1**3)
print x2-x1
A2 = A2 * inverse_mod(x2 - x1, mod)
print "modinv:",repr(inverse_mod(x2 - x1, mod))
A2 %= mod

B2 = y1**2 - x1**3 - A2 * x1
B2 %= mod

print "A, B :", A2, B2

# new curve
E2 = EllipticCurve(IntegerModRing(p**2), [A2, B2])

# calculate dlog
g2s = (p - 1) * E2(x1, y1)
v2s = (p - 1) * E2(x2, y2)

x1s, y1s = map(long, g2s.xy())
x2s, y2s = map(long, v2s.xy())

dx1 = (x1s - x1) / p % p
dx2 = (y1s - y1) / p
dy1 = (x2s - x2)
dy2 = (y2s - y2) % p

print "%d, %d, %d, %d, %d" % (dx1, dy1, dx2, dy2, p)

m = dy1 * inverse_mod(dx1, p) * dx2 * inverse_mod(dy2, p)
print dy1, inverse_mod(dx1, p), dx2, inverse_mod(dy2, p)
m %= p

print m
