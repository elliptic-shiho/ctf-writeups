RSA? (Crypto 2pts)
===================================

Challenge file(s) is here: [rsa.zip](rsa.zip)

## Writeup
We have RSA public key and encrypted flag. we know public key is here:

$$
\begin{aligned}
e &= 3\\\\
n &= 23292710978670380403641273270002884747060006568046290011918413375473934024039715180540887338067
\end{aligned}
$$

We noticed n is sufficiently small (only 314-bits long). so, we try factoring it.

Factoring result is here:

$$
\begin{aligned}
p &= 26440615366395242196516853423447\\\\
q &= 27038194053540661979045656526063\\\\
r &= 32581479300404876772405716877547\\\\
n &= p\times q\times r
\end{aligned}
$$

Wow, This is **Multi-Prime RSA**. but we know all private factor of $n$, we can compute $\phi(n) = \phi(p)\phi(q)\phi(r) = (p-1)(q-1)(r-1)$.

$$
\begin{aligned}
\phi(n) &= 23292710978670380403641273270000427421848709005360280557445800298810723014218767619832560713992\\\\
&= 2^3\times 3^3 \times 1468923075910846788695380745747 \times 5430246550067479462067619479591\times 13519097026770330989522828263031
\end{aligned}
$$

BTW, In RSA cryptosystem needs to compute Private Exponent $d$. It is defined by $d\equiv e^{-1}\mod \phi(n)$.

However, In this case $e$ is _not_ invertible (because $\gcd(e, \phi(n)) \ne 1$). so, we can't compute private exponent!

well, we use **Chinese-Remainder-Theorem** (CRT) to solve challenge. CRT statement is here:

_**Chinese Remainder Theorem**_ (with 3-equations)
> Let $N\_1$, $N\_2$, and $N\_3$ be coprime, and let $x\_1$, $x\_2$, and $x\_3$ be integers. Then the system of equations
> $$\begin{aligned} x &\equiv x\_1 \mod N\_1\\\\ x &\equiv x\_2 \mod N\_2\\\\x &\equiv x\_3 \mod N\_3\end{aligned}$$
> has a unique solution for $x$ modulo $N\_1N\_2N\_3$.

Let $m$ be a flag and let ciphertext $c$ be a RSA ciphertext of $m$. Indeed $c = m^3\mod n$. We know $n = pqr$ and $p$, $q$ and $r$ is coprime. and we can compute modular cubic root over finite field. 

We use following procedure to solve this:

1. Compute $x\_1 \equiv \sqrt[3]{c}\mod p$
2. Compute $x\_2 \equiv \sqrt[3]{c}\mod q$
3. Compute $x\_3 \equiv \sqrt[3]{c}\mod r$
4. Apply CRT for $(x\_1, p), (x\_2, q), (x\_3, r)$
5. we got $m$. this is flag.

Implementation is here:

```python
from scryptos import *

p1 = 32581479300404876772405716877547
p2 = 27038194053540661979045656526063
p3 = 26440615366395242196516853423447
n = p1*p2*p3
e = 3

c = int(open("flag.enc", "rb").read().encode("hex"), 16)

# from User's Guide to PARI/GP, nth_root function
sqrtnall = 'sqrtnall(x,n)={my(V,r,z,r2);r=sqrtn(x,n,&z);if(!z,error("Impossible case in sqrtn"));if(type(x)=="t_INTMOD"||type(x)=="t_PADIC",r2 = r*z;n=1;while(r2!=r,r2*=z;n++));V=vector(n);V[1]=r;for(i=2,n,V[i]=V[i-1]*z);V}'

c1 = eval(parigp([sqrtnall, "Vec(liftall(sqrtnall(Mod(%d, %d), 3)))" % (c, p1)]))
c2 = eval(parigp([sqrtnall, "Vec(liftall(sqrtnall(Mod(%d, %d), 3)))" % (c, p2)]))
c3 = eval(parigp([sqrtnall, "Vec(liftall(sqrtnall(Mod(%d, %d), 3)))" % (c, p3)]))

"""
c1 = [6149264605288583791069539134541, 13404203109409336045283549715377, 13028011585706956936052628027629]
c2 = [19616973567618515464515107624812]
c3 = [13374868592866626517389128266735, 7379361747422713811654086477766, 5686385026105901867473638678946]
"""

for x in c1:
  for y in c2:
    for z in c3:
      crt = chinese_remainder_theorem([(x, p1), (y, p2), (z, p3)])
      d = hex(crt, 2)[2:].decode("hex")
      if "0ctf" in d:
print d[d.find("0ctf"):].strip()
```

Execute log:

```
Mon Mar 14 09:14:36 JST 2016 ~/ctf/0ctf-2016/crypto2 Battery 0: Full, 100%
> python solve.py 
0ctf{HahA!Thi5_1s_n0T_rSa~}
```

Flag: `0ctf{HahA!Thi5_1s_n0T_rSa~}`


# References