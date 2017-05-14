hmac\_crc (Crypto 5pts)
==================

Challenge file(s) is here: [0c7433675c3c555afb77271d6a549bf5d941d2ab](0c7433675c3c555afb77271d6a549bf5d941d2ab)

$$
\def\Ipad{\mathbb{I}}
\def\Opad{\mathbb{O}}
$$

## Writeup
This is HMAC Hash using CRC64. We need to get HMAC Key from a hash/plaintext pair.

### Notation:
* $P(x)$ : Generator polynomial of $GF(2^64)$
* $x$    : Generator of $GF(2^64)$
* $\Ipad$ : Inner Padding
* $\Opad$ : Outer Padding
* $C$    : Some Constant

All Operation is operation over GF(264). i.e. Addition means XOR operation.

---

$CRC\_{64}$ function is defined by $CRC _ {64}(m) = mx^{64} + C\mod P(X)$. HMAC using CRC64 is Defined as follows:

$$
\begin{aligned}
HMAC\\\_CRC \_ {64}(k, m) &= CRC \_ {64}((k + \Opad) x^{64} + CRC \_ 64(k + \Ipad) x^M + m)\mod P(x)\\\\
                 &= CRC \_ {64}(k x^{64} + \Opad x^{64} + CRC \_ {64}(k x^M + \Ipad x^M + m))\mod P(x)
\end{aligned}
$$

where $M$ is bit length of $m$.

We need only $k$. so we transform that equation to solve for $k$.

$$
\begin{aligned}
HMAC\\\_CRC\_{64}(k, m) &= CRC\_{64}(kx^{64} + \Opad x^{64} + k x^{M+64} + \Ipad x^{M+64} + m x^{64} + C)\mod P(x)\\\\
                 &= k x^{128} + \Opad x^{128} + k x^{M+128} + \Ipad x^{M+128} + m x^{128} + C x^{64} + C\mod P(x)
\end{aligned}
$$

$$
HMAC\\\_CRC\_{64}(k, m) = k (x^{128} + x^{M+128}) + \Opad x^{128} + \Ipad x^{M+128} + m x^{128} + C (x^{64} + 1)\mod P(x)
$$

$$
\begin{aligned}
k (x^{128} + x^{M+128}) &= \Opad x^{128} + \Ipad x^{M+128} + m x^{128} + C (x^{64} + 1) + HMAC\\\_CRC\_{64}(k, m)\mod P(x)\\\\
k &= (x^{128} (\Opad + \Ipad x^M + m) + C (x^{64} + 1) + HMAC\\\_CRC\_{64}(k, m)) / x^{128} (x^M + 1)\mod P(x)
\end{aligned}
$$

We got a equation for $k$. In addition, we can calculate right-side of equation. Therefore, We got $k$. Furthermore, we need only k to got flag. so, we got flag.


Flag: `BKPCTF{0xd2db2b8b9002841f}`

[get_key.sage](get_key.sage)

```python
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
```

## References
* [MMA CTF 2015 - Motto Mijikai Address (Crypto/Web 100+300) | More Smoked Leet Chicken](http://mslc.ctf.su/wp/mma-ctf-2015-motto-mijikai-address-cryptoweb-100300/)
