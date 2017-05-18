Alpha Complex 1 (Crypto 300pts)
===================================

Challenge file(s) is here: [AlphaComplex1.zip](AlphaComplex1.zip)

## Writeup
We have RSA-like cryptosystem. in fact, this is RSA over the **Gaussian Integer**.

First, We construct a RSA over the Gaussian Integer. Gaussian Integer is extended field of $\mathbb{Z}$:

$$
\mathbb{Z}[i]\quad\mathrm{where}\ i^2 = -1
$$

It is similar to Complex Field. indeed, some properties of complex number is similarly holds (norm, inverse calculation, ...). However, this is useless for cryptosystem. so, we use $\mathbb{Z} / n\mathbb{Z}$ as base field!

$$
\mathbb{Z}/n\mathbb{Z}[i]\quad\mathrm{where}\ i^2 = -1\\\\
\iff (\mathbb{Z}/n\mathbb{Z}[x])/(x^2 + 1)
$$

Suppose $n = pq$ where $p$, $q$ is prime, we know $\\#(\mathbb{Z}/n\mathbb{Z}[i])^\* = (p^2 - 1)(q^2-1)$ and for all $x \in (\mathbb{Z}/n\mathbb{Z}[i])^\*$ holds $x^{\\#(\mathbb{Z}/n\mathbb{Z}[i])} = 1$. so we pick $e$ which has e is relatively prime to $\\#(\mathbb{Z}/n\mathbb{Z}[i])$, we can construct RSA-like cryptosystem using $d\equiv e^{-1} \mod \\#(\mathbb{Z}/n\mathbb{Z}[i])$.

----

In this challenge, we know ciphertext for plaintext $m = x + yi$ with known $y$. The public key is generate randomly at every connection, but $e$ was same number at all public-key. we use broadcast-attack, maybe? - No, we don't use broadcast-attack. because $e$ is big($0\mathrm{x}\\!1337 = 4919$).

We know ciphertext $c \equiv m^e \mod n$ is represented as follows:

$$
\begin{aligned}
c&\equiv m^e \mod n\\\\
&= (x + yi)^e \mod n\\\\
&= x' + y'i
\end{aligned}
$$

We know $y$. Moreover $x'$ and $y'$ is able to represent by polynomial of $x$ and $y$  for each. so, Let $x$ be a indeterminate variable, we have **two** univariate polynomials over $\mathbb{Z}/n\mathbb{Z}$. we denote those polynomials by $x'(x)$ and $y'(x)$. then $x'(x)$ and $y'(x)$ has **same roots**.

Finally, we compute $g(x) = \gcd(x'(x), y'(x))$ and compute roots of $g(x)$. From actual results, $g(x)$ is linear polynomial. so, we can compute roots of that easily.

[solve.py](solve.py)

```python
from roputils import Proc
from scryptos import *
from itertools import *
from scryptos.wrapper.parigp import set_gp_memalloc_size
import hashlib
import string

"""
p = Proc(host="ac1.pwn.seccon.jp", port=31337)
p.read_until("prefix: ")
prefix = p.readline().strip()

for x in product(string.printable, repeat=16):
  if hashlib.sha1(prefix + "".join(x)).hexdigest()[:5] == "00000":
    print repr("".join(x))
    p.write("".join(x))
    break

p.interact()
"""

n=85913081289093219314279414199911465975826993645255099108711586608621408324014021787874679342549096244529235932586709419001980975828323391172145792836579199724931013820275272212225060955572024694062338567977737291578560186799120393834106834908399529875670248954839313057047385142493049783648656051266200679539
e=4919
y = 109823
c = [7343836202785630351137266838524322944196218302395251619708621278307028537163796438017957382700469113231065259304446768938126604072292056771961938036191098188575796583664333657843736555276833510818568441987946108063187232471647366362615025383861649632471364729561625292996669722943156467517886902288663376530, 78805899264837359127934361704813657728028495417225274572152492070945573113436642696718282244964253253778011656522230634179033236200136511174534988947704854580923553360972284483058586894235508191718244903134393200667608836101550862777775562092269758188940713875300829889369313825258345746287364875750361695453]

# sage: z = expand((x+y*I)^4919)
# z.real() -> out2_1.txt
# z.imag() -> out2_2.txt
d1 = open("out2_1.txt").read()
d2 = open("out2_2.txt").read()
# allocate 512MB for calculation
set_gp_memalloc_size(1024 * 1024 * 512)
expr = []
expr += ["n = %d" % n]
expr += ["y = %d" % y]
expr += ["re = Mod(%s, n)" % d1]
expr += ["im = Mod(%s, n)" % d2]
expr += ["g = gcd(Mod(re-%d, n), Mod(im-%d, n))" % (c[0], c[1])]
expr += ["liftall(-Vec(g / Vec(g)[1])[2])"]
x = eval(parigp(expr))
print repr(long_to_bytes(x))
```

Execute Log:

```
Tue Dec 13 11:25:28 JST 2016 ~/ctf/2016/seccon-quals-2016/cr/alpha_complex_1 100%
> time python solve.py 
  ***   Warning: new stack size = 536870912 (512.000 Mbytes).
'yK#\xaf\x1c\xfc\xae\x02\x92\xd7\xa8ES\x8fp\x92\xe0\xaf\xf7\xe6\xcc\x9fR%\xe6T\x90\x9cB4\xbcX;\xe4\x11?\xb3\x8d\xb2)q\x1eM\\\xfav\x86\xe0\x0c\xe1G\x05yo"8/\xa7\xbc\xb1\xefSC\x96Y\xe1\x14\x8f\x11\xb4:\xe3\xc3G\xb4\xe5\xb1\x00SECCON{3907c554839a6462b920f55274c159eea1949086}'

real    0m25.407s
user    0m25.212s
sys     0m0.164s
```

Flag: `SECCON{3907c554839a6462b920f55274c159eea1949086}`


# References