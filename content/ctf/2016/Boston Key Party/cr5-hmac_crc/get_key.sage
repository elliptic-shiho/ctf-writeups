from sage.all import *

def ntopoly(npoly):
    return sum(c*X**e for e, c in enumerate(Integer(npoly).bits()))

def polyton(poly):
    return sum(int(poly[i])*(1 << i) for i in xrange(poly.degree() + 1))

X = GF(2).polynomial_ring().gen()

N = 64
poly = (2**64) + 0xeff67c77d13835f7
poly = ntopoly(poly)

CONST = 0xabaddeadbeef1dea
CONST = ntopoly(CONST)

I = ntopoly(int(("36"*8), 16))
O = ntopoly(int(("5C"*8), 16))

def crc(m, l):
  return polyton((ntopoly(m) * X^64 + CONST) % poly)

def hmac_crc64(m, k, l):
  m = ntopoly(m)
  k = ntopoly(k)
  a = crc(polyton((k + I) * X^l + m), l+64)
  b = crc(polyton((k + O) * X^64 + ntopoly(a)), 128)
  return b

msg = "zupe zecret"
M = 88
plain = int(msg.encode("hex"), 16)
plain = ntopoly(plain)
mesg = 0xa57d43a032feb286
mesg = ntopoly(mesg)

part = (CONST * (X^64 + 1) + X^128 * (plain + I * X^M + O) + mesg) % poly

print part
k = polyton(part * inverse_mod((X^128 + X^(128 + M)), poly) % poly)
print "[+] Key = 0x%016x" % k
print hex(hmac_crc64(polyton(plain), k, 88))
flag = "BKPCTF"
l = len(flag)*8
flag = int(flag.encode("hex"), 16)
print "[+] Flag is BKPCTF{0x%s}" % hex(hmac_crc64(flag, k, l))
