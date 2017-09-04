Authentication&Secrecy (Crypto 1000pts)
===================================

Challenge file(s) is here: [rsa_66080e40bb4a92fd4358d7236cda7c91.zip](rsa_66080e40bb4a92fd4358d7236cda7c91.zip)

## Writeup
We have a *weired* RSA-Key Generator, some communication program, a ciphertext, and transmitter/receiver's RSA Public-Key.

First, We need to known **Dual RSA**. Dual RSA public key has a Public Exponent $e$ and *Different* Modulus $N_1$, $N_2$ satisfy $ed\equiv 1\mod\phi(N_1)$ and $ed\equiv 1\mod\phi(N_2)$ where $d$ is a corresponding private exponent of $e$. This RSA scheme invented by Sun, Wu, Ting, and Hinek [1], this is called **Dual RSA**.

We searching paper and found [2]. So we implemented it. However, that paper is maybe not perfect...

* Actual results does not following to paper computation. (Especially determinant of LLL-reduced Lattice determinant)
* We have a small question in this paragraph([2] p.10):

> Moreover, the more integer equations corresponding to the vectors we choose, the less time calculating Gröbner basis takes. For the previous example, when we chose 50 of these integer equations, the calculation of Gröbner basis took 1933.664 seconds. However, when we chose all of them, it took only 16.700 s.

Maybe, **all of them** means **all of got polynomials**. However, "The vector is basis of LLL-reduced lattice. therefore It satisfy Howgrave-graham's bound." is **not** always true. So we can't compute a solution using all of them! (Is this wrong point of paper?)

Therefore, we select maybe shorter vectors by heuristic. Indeed, we use a 19 vectors from the second one.

[solve.sage](solve.sage)

```python
from sage.all import *
import math
import itertools
import alice_public_key
import bob_public_key

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    print a

def dual_rsa_liqiang_et_al(e, n1, n2, delta, mm, tt):
  '''
  Attack to Dual RSA: Liqiang et al.'s attack implementation

  References:
    [1] Liqiang Peng, Lei Hu, Yao Lu, Jun Xu and Zhangjie Huang. 2016. "Cryptanalysis of Dual RSA"
  '''
  N = (n1+n2)/2
  A = ZZ(floor(N^0.5))

  _XX = ZZ(floor(N^delta))
  _YY = ZZ(floor(N^0.5))
  _ZZ = ZZ(floor(N^(delta - 1./4)))
  _UU = _XX * _YY + 1

  # Find a "good" basis satisfying d = a1 * l'11 + a2 * l'21
  M = Matrix(ZZ, [[A, e], [0, n1]])
  B = M.LLL()
  l11, l12 = B[0]
  l21, l22 = B[1]
  l_11 = ZZ(l11 / A)
  l_21 = ZZ(l21 / A)

  modulo = e * l_21
  F = Zmod(modulo)

  PR = PolynomialRing(F, 'u, x, y, z')
  u, x, y, z = PR.gens()

  PK = PolynomialRing(ZZ, 'uk, xk, yk, zk')
  uk, xk, yk, zk = PK.gens()

  # For transform xy to u-1 (unravelled linearlization)
  PQ = PK.quo(xk*yk+1-uk)

  f = PK(x*(n2 + y) - e*l_11*z + 1)

  fbar = PQ(f).lift()

  # Polynomial construction
  gijk = {}
  for k in xrange(0, mm + 1):
    for i in xrange(0, mm-k + 1):
      for j in xrange(0, mm-k-i + 1):
        gijk[i, j, k] = PQ(xk^i * zk^j * PK(fbar) ^ k * modulo^(mm-k)).lift()

  hjkl = {}
  for j in xrange(1, tt + 1):
    for k in xrange(floor(mm / tt) * j, mm + 1):
      for l in xrange(0, k + 1):
        hjkl[j, k, l] = PQ(yk^j * zk^(k-l) * PK(fbar) ^ l * modulo^(mm-l)).lift()

  monomials = []
  for k in gijk.keys():
    monomials += gijk[k].monomials()
  for k in hjkl.keys():
    monomials += hjkl[k].monomials()

  monomials = sorted(set(monomials))[::-1]
  assert len(monomials) == len(gijk) + len(hjkl) # square matrix?
  dim = len(monomials)

  # Create lattice from polynmial g_{ijk} and h_{jkl}
  M = Matrix(ZZ, dim)
  row = 0
  for k in gijk.keys():
    for i, monomial in enumerate(monomials):
      M[row, i] = gijk[k].monomial_coefficient(monomial) * monomial.subs(uk=_UU, xk=_XX, yk=_YY, zk=_ZZ)
    row += 1
  for k in hjkl.keys():
    for i, monomial in enumerate(monomials):
      M[row, i] = hjkl[k].monomial_coefficient(monomial) * monomial.subs(uk=_UU, xk=_XX, yk=_YY, zk=_ZZ)
    row += 1

  matrix_overview(M)
  print '=' * 128

  # LLL
  B = M.LLL()

  matrix_overview(B)

  # Construct polynomials from reduced lattices
  H = [(i, 0) for i in xrange(dim)]
  H = dict(H)
  for j in xrange(dim):
    for i in xrange(dim):
      H[i] += PK((monomials[j] * B[i, j]) / monomials[j].subs(uk=_UU, xk=_XX, yk=_YY, zk=_ZZ))
  H = H.values()

  PQ = PolynomialRing(QQ, 'uq, xq, yq, zq')
  uq, xq, yq, zq = PQ.gens()

  # Inversion of unravelled linearlization
  for i in xrange(dim):
    H[i] = PQ(H[i].subs(uk=xk*yk+1))

  # Calculate Groebner basis for solve system of equations
  '''
  Actually, These polynomials selection (H[1:20]) is heuristic selection.
  Because they are "short" vectors. We need a short vector less than
  Howgrave-Graham bound. So we trying test parameter(at [1]) and decided it.
  '''
  I = Ideal(*H[1:20])
  g = I.groebner_basis('giac')[::-1]
  mon = map(lambda t: t.monomials(), g)

  PX = PolynomialRing(ZZ, 'xs')
  xs = PX.gen()

  x_pol = y_pol = z_pol = None

  for i in xrange(len(g)):
    if mon[i] == [xq, 1]:
      print g[i] / g[i].lc()
      x_pol = g[i] / g[i].lc()
    elif mon[i] == [yq, 1]:
      print g[i] / g[i].lc()
      y_pol = g[i] / g[i].lc()
    elif mon[i] == [zq, 1]:
      print g[i] / g[i].lc()
      z_pol = g[i] / g[i].lc()

  if x_pol is None or y_pol is None or z_pol is None:
    print '[-] Failed: we cannot get a solution...'
    return

  x0 = x_pol.subs(xq=xs).roots()[0][0]
  y0 = y_pol.subs(yq=xs).roots()[0][0]
  z0 = z_pol.subs(zq=xs).roots()[0][0]

  # solution check
  assert f(x0*y0+1, x0, y0, z0) % modulo == 0

  a0 = z0
  a1 = (x0 * (n2 + y0) + 1 - e*l_11*z0) / (e*l_21)

  d = a0 * l_11 + a1 * l_21
  return d

if __name__ == '__main__':
  delta = 0.334
  mm = 4
  tt = 2

  n1 = alice_public_key.N1
  n2 = alice_public_key.N2
  e = alice_public_key.e

  d1 = dual_rsa_liqiang_et_al(e, n1, n2, delta, mm, tt)
  print '[+] d for alice = %d' % d1

  n1 = bob_public_key.N1
  n2 = bob_public_key.N2
  e = bob_public_key.e

  d2 = dual_rsa_liqiang_et_al(e, n1, n2, delta, mm, tt)
  print '[+] d for bob = %d' % d2

  d = ZZ(open('send.bak', 'rb').read().encode('hex'), 16)
  m = ZZ(Mod(ZZ(Mod(d, bob_public_key.N1)^d2), alice_public_key.N2)^alice_public_key.e)
  m_ = hex(m)
  if len(m_) % 2 == 1:
    m_ = '0' + m_
  print repr(m_.decode('hex'))
```

```
Sun Jun  4 17:25:20 JST 2017 ~/Downloads/rsa 57%
> sage solve.sage
00                                                                                           X X X         
01                       X X X X     X X X     X X     X                                                   
02                                                                                                   X     
03                                                                                       X X   X           
04                                                     X                                                   
05                                                                             X X   X                     
06                                               X                                                         
07                                   X X X     X X     X                                                   
08                                                                                               X         
09                             X                                                                           
10                                           X X X   X X   X                                               
11                                                                     X X X   X X   X                     
12                                                               X X       X                               
13                     X X X X X   X X X X   X X X   X X   X                                               
14                                                                               X                         
15                                                                                         X               
16                                                           X X X X   X X X   X X   X                     
17                                             X X     X                                                   
18                         X X X       X X       X                                                         
19                                                                                     X X X X X X         
20                                                                                   X                     
21                                       X                                                                 
22                                                                                                 X X X   
23                                                   X X   X                                               
24                                                                                                       X 
25                                                                       X X     X                         
26                                 X X X X   X X X   X X   X                                               
27                                                             X X X     X X     X                         
28                                                                                             X           
29                                                                         X                               
30                                     X X       X                                                         
31                                                         X                                               
32                                                                 X                                       
33                           X X         X                                                                 
34                                                                                                     X   
35                                                                           X     X         X         X   
36   X X X X   X X X               X X     X X     X     X             X X     X             X X X         
37                   X                                                                                     
38                                         X       X     X             X X     X             X X X         
39         X                                                                                               
40                                                                   X       X     X   X X   X     X X X   
41                 X X                               X                               X                     
42                                                 X     X                     X                 X         
43 X X X X X X X X X   X X X     X X X     X X     X     X   X X X     X X     X       X X X X X X         
44       X X       X                                     X                                                 
45     X X X     X X                         X     X     X                     X                 X         
46           X X X X X X X X X     X X X     X X     X       X X X X   X X X   X X   X                     
47                                                       X                                                 
48               X X X                       X X     X                         X X   X                     
49                                                                                 X                       
50             X X X X             X X X     X X     X                 X X X   X X   X                     
51                               X         X       X     X   X X X     X X     X       X X X X X X         
================================================================================================================================
00           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
01           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
02           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
03           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
04           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
05           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
06           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
07           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
08           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
09           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
10           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
11           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
12           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
13           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
14           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
15           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
16           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
17           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
18           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
19           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
20           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
21           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
22           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
23           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
24           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
25           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
26           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
27           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
28           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
29           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
30                                                                                     X X X X X X         
31                                                                                                 X X X   
32                                                                                     X X X X X X         
33                                                                                                       X 
34                                                                                                 X X X   
35                                                                                     X X X X X X         
36           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
37                                                                                     X X X X X X         
38                                                                                                 X X X   
39           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
40           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
41           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
42           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
43           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
44           X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   
45           X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   
46                                                                   X       X     X   X X X X X X X X X   
47 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
48 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
49 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
50 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
51 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
// Giac share root-directory:/usr/lib/sagemath//local/share/giac/
// Giac share root-directory:/usr/lib/sagemath//local/share/giac/
Help file /usr/lib/sagemath//local/share/giac/doc/fr/aide_cas not found
Added 0 synonyms

// Groebner basis computation time 0.078343 Memory 0.203556M
zq - 61140973128283104870537838
yq + 8407381940047047520115757582815846503329543922728691632771955750973825326885843225149837173391777915655641415687874176390003516776943256011646863368652171
xq - 7919404215158064938288843511454164188553173880252141278917970265662701060110142326401620725399524328585
[+] d for alice = 5175600911919085867132554056898403461564674635722636201137442473058967531761922627791714515652901422461
00                                                                                           X X X         
01                       X X X X     X X X     X X     X                                                   
02                                                                                                   X     
03                                                                                       X X   X           
04                                                     X                                                   
05                                                                             X X   X                     
06                                               X                                                         
07                                   X X X     X X     X                                                   
08                                                                                               X         
09                             X                                                                           
10                                           X X X   X X   X                                               
11                                                                     X X X   X X   X                     
12                                                               X X       X                               
13                     X X X X X   X X X X   X X X   X X   X                                               
14                                                                               X                         
15                                                                                         X               
16                                                           X X X X   X X X   X X   X                     
17                                             X X     X                                                   
18                         X X X       X X       X                                                         
19                                                                                     X X X X X X         
20                                                                                   X                     
21                                       X                                                                 
22                                                                                                 X X X   
23                                                   X X   X                                               
24                                                                                                       X 
25                                                                       X X     X                         
26                                 X X X X   X X X   X X   X                                               
27                                                             X X X     X X     X                         
28                                                                                             X           
29                                                                         X                               
30                                     X X       X                                                         
31                                                         X                                               
32                                                                 X                                       
33                           X X         X                                                                 
34                                                                                                     X   
35                                                                           X     X         X         X   
36   X X X X   X X X               X X     X X     X     X             X X     X             X X X         
37                   X                                                                                     
38                                         X       X     X             X X     X             X X X         
39         X                                                                                               
40                                                                   X       X     X   X X   X     X X X   
41                 X X                               X                               X                     
42                                                 X     X                     X                 X         
43 X X X X X X X X X   X X X     X X X     X X     X     X   X X X     X X     X       X X X X X X         
44       X X       X                                     X                                                 
45     X X X     X X                         X     X     X                     X                 X         
46           X X X X X X X X X     X X X     X X     X       X X X X   X X X   X X   X                     
47                                                       X                                                 
48               X X X                       X X     X                         X X   X                     
49                                                                                 X                       
50             X X X X             X X X     X X     X                 X X X   X X   X                     
51                               X         X       X     X   X X X     X X     X       X X X X X X         
================================================================================================================================
00           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
01           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
02           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
03           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
04           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
05           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
06           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
07           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
08           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
09           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
10           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
11           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
12           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
13           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
14           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
15           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
16           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
17           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
18           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
19           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
20           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
21           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
22           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
23           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
24           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
25           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
26           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
27           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
28           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
29           X X X X X X X X X X   X X X X   X X X   X X   X X X X X   X X X   X X   X                     
30                                                                                                       X 
31                                                                                                 X X X   
32                                                                                     X X X X X X         
33                                                                                                 X X X   
34                                                                                     X X X X X X         
35                                                                                                 X X X   
36                                                                                     X X X X X X         
37                                                                                     X X X X X X         
38                                                                                     X X X X X X         
39           X X X X X X X X X X X X X X X   X X X   X X   X X X X X   X X X   X X   X X X X X X X         
40                                                                                     X X X X X X         
41           X X X X X X X X X X X X X X X X X X X   X X   X X X X X   X X X   X X   X X X X X X X         
42           X X X X X X X X X X X X X X X X X X X X X X   X X X X X   X X X   X X   X X X X X X X         
43           X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
44                                                                   X                 X X X X X X X X X   
45                                                                   X       X         X X X X X X X X X   
46                                                                   X       X     X   X X X X X X X X X   
47 X         X X X X X X X X X X X X X X X   X X X   X X   X X X X X   X X X   X X   X X X X X X X         
48 X X       X X X X X X X X X X X X X X X X X X X   X X   X X X X X   X X X   X X   X X X X X X X         
49 X X X     X X X X X X X X X X X X X X X X X X X X X X   X X X X X   X X X   X X   X X X X X X X         
50 X X X X   X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         
51 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   X X X   X X   X X X X X X X         

// Groebner basis computation time 0.072379 Memory 0.207636M
zq - 143145407796820419691921116
yq + 15352514463512545429409965865844475808601686641370002128014114724769375278708730306689584590428850472250709065231080455974049923847598927444538942404357569
xq - 6046141695019859939882486410408884140504805414384153622531566562183840881867683963179563847287957629047
[+] d for bob = 6238257262527822874210763706228720299138117835182790315423266965504012841647775705014056496108579279573
"\x7f\xc6\x0b/\x89\xd5\xf4u\x8d\xcf\nH&\x04V\xff\xbd:\xe8\x02\x9c\xd5\xca\xe8\x9bFo'?\x1d\x195Ci\xa3\xb8\xa0w\x02W\x9d*@$\xc5\xaa=\xba\x00Well, OK... Here is what you what: flag{eNj0y_th3_CrYpt4na1ysi5_0F_Dual_RSA~~}"
```

Flag: `flag{eNj0y_th3_CrYpt4na1ysi5_0F_Dual_RSA~~}`

I solved this challenge, only one of the participated teams!

# References
* [1] Hung-Min Sun, Mu-En Wu, Wei-Chi Ting, and M. Jason Hinek. 2007. [_Dual RSA and its security analysis_](https://www.researchgate.net/publication/3086435_Dual_RSA_and_its_security_analysis)
* [2] Liqiang Peng, Lei Hu, Yao Lu, Jun Xu, and Zhangjie Huang. 2016. [_Cryptanalysis of Dual RSA_](https://link.springer.com/article/10.1007/s10623-016-0196-5)
