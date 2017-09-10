Gracias (Crypto 287pts)
===================================

Challenge file(s) is here: [gracias_5cf9b9e3f11e527b1b4cd99d67e19bee34d1aab4](gracias_5cf9b9e3f11e527b1b4cd99d67e19bee34d1aab4)

## Writeup
This Challenge is Small Secret Exponent Attack against Multi-Prime RSA. 

**Key-Generation Algorithm of Multi-Prime RSA** (in this challenge): Generate prime $p$, $q$ and $r$ and Compute $n = pqr$ (Note : $\phi(n) = (p-1)(q-1)(r-1)$). 
Select $d\approx 2^{256}$ and Compute $e\equiv d^{-1}\pmod{\phi(n)}$. Public-Key is $(n, e)$ and Private-Key is $d$.

In this generation algorithm, $d$ is quite small relatively to $n$. In fact, we know $\log\_2(d)/\log\_2(n) \approx 0.167$. In normal RSA, We know some attacks for small private key. Especially Wiener's Attack [1]\($d\lt n^{0.25}$) and Boneh-Durfee's Attack [2]\($d\lt n^{0.292}$) is most famous attacks. However, We didn't know any tool/attacks for Multi-Prime RSA.

Now, We consider to breaking Multi-Prime RSA which has small private key. First, We know following equation:

$$
ed \equiv 1 \pmod{\phi(n)}.
$$

Using $\phi(n) = (p-1)(q-1)(r-1) = n - (pq+pr+qr) + p+q+r - 1$, $A = n-1$ and $s = -(pq+pr+qr)+p+q+r$, We know below that:

$$
\begin{aligned}
ed &= 1 + k\phi(n)\\\
&= 1 + k(A + s)\\
\end{aligned}
$$

$$
\therefore k(A + s) + 1\equiv 0\pmod e. \tag{1}
$$

where $k$ is a integer.

Recall Boneh-Durfee's Attack, That are use the equation like $x(A + y) - 1 \equiv 0\mod e$. It is very similar to our equation $(1)$!

We consider upper-bound of the solution. Let $(x, y) = (k, s)$ then we have $|x|\leq X = n^{0.167}$ and $|y|\leq Y = n^{2/3}$. 

Finally, We write solver, based [David Wong's Boneh-Durfee's Attack implementation](https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage), using consideration so far.

[boneh_durfee.sage](boneh_durfee.sage)
```python
from sage.all import *
# Original: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage

dimension_min = 7

def remove_unhelpful(BB, monomials, bound, current):
  if current == -1 or BB.dimensions()[0] <= dimension_min:
    return BB
  for ii in range(current, -1, -1):
    if BB[ii, ii] >= bound:
      affected_vectors = 0
      affected_vector_index = 0
      for jj in range(ii + 1, BB.dimensions()[0]):
        if BB[jj, ii] != 0:
          affected_vectors += 1
          affected_vector_index = jj
      if affected_vectors == 0:
        #print "* removing unhelpful vector", ii
        BB = BB.delete_columns([ii])
        BB = BB.delete_rows([ii])
        monomials.pop(ii)
        BB = remove_unhelpful(BB, monomials, bound, ii-1)
        return BB
      elif affected_vectors == 1:
        affected_deeper = True
        for kk in range(affected_vector_index + 1, BB.dimensions()[0]):
          if BB[kk, affected_vector_index] != 0:
            affected_deeper = False
        if affected_deeper and abs(bound - BB[affected_vector_index, affected_vector_index]) < abs(bound - BB[ii, ii]):
          #print "* removing unhelpful vectors", ii, "and", affected_vector_index
          BB = BB.delete_columns([affected_vector_index, ii])
          BB = BB.delete_rows([affected_vector_index, ii])
          monomials.pop(affected_vector_index)
          monomials.pop(ii)
          BB = remove_unhelpful(BB, monomials, bound, ii-1)
          return BB
  return BB

def boneh_durfee_small_roots(pol, modulus, mm, tt, XX, YY):
    PR.<u, x, y> = PolynomialRing(ZZ)
    Q = PR.quotient(x*y + 1 - u) # u = xy + 1
    polZ = Q(pol).lift()
    UU = XX*YY + 1
    gg = []
    for kk in range(mm + 1):
      for ii in range(mm - kk + 1):
        xshift = x^ii * modulus^(mm - kk) * polZ(u, x, y)^kk
        gg.append(xshift)
    gg.sort()
    monomials = []
    for polynomial in gg:
      for monomial in polynomial.monomials():
        if monomial not in monomials:
          monomials.append(monomial)
    monomials.sort()
    for jj in range(1, tt + 1):
      for kk in range(floor(mm/tt) * jj, mm + 1):
        yshift = y^jj * polZ(u, x, y)^kk * modulus^(mm - kk)
        yshift = Q(yshift).lift()
        gg.append(yshift)
        monomials.append(u^kk * y^jj)
    nn = len(monomials)
    BB = Matrix(ZZ, nn)
    for ii in range(nn):
      BB[ii, 0] = gg[ii](0, 0, 0)
      for jj in range(1, ii + 1):
        if monomials[jj] in gg[ii].monomials():
          BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU,XX,YY)
    BB = remove_unhelpful(BB, monomials, modulus^mm, nn-1)
    nn = BB.dimensions()[0]
    if nn == 0:
      print "failure"
      return 0,0
    BB = BB.LLL()
    PR.<w,z> = PolynomialRing(ZZ)
    pol1 = pol2 = 0
    for jj in range(nn):
      pol1 += monomials[jj](w*z+1,w,z) * BB[0, jj] / monomials[jj](UU,XX,YY)
      pol2 += monomials[jj](w*z+1,w,z) * BB[1, jj] / monomials[jj](UU,XX,YY)
    PR.<q> = PolynomialRing(ZZ)
    rr = pol1.resultant(pol2)
    if rr.is_zero() or rr.monomials() == [1]:
      print "the two first vectors are not independant"
      return 0, 0
    rr = rr(q, q)
    soly = rr.roots()
    if len(soly) == 0:
      print "Your prediction (delta) is too small"
      return 0, 0
    soly = soly[0][0]
    ss = pol1(q, soly)
    solx = ss.roots()[0][0]
    return solx, soly

def boneh_durfee(n, e):
  delta = RR(0.167) # d ~ n^0.167
  m = 5
  t = round((1-2*delta) * m)
  X = ZZ(2*floor(n^delta))
  # we have n = p^2q. so `phi(n) = n + {-(pq+pr+qr) + p+q+r)} - 1`.
  # we reconsidered boneh-durfee's attack then we have `x(A+y) + 1 = 0 mod e` where `A = (n-1)`
  # and (x, y) = (k, -(pq+pr+qr)+p+q+r). 
  Y = ZZ(floor(n^(2/3)))
  P.<x,y> = PolynomialRing(ZZ)
  A = ZZ((n-1)/2)
  pol = 1 + x * (A + y)
  solx, soly = boneh_durfee_small_roots(pol, e, m, t, X, Y)
  print solx, soly
  if solx > 0:
    return int(pol(solx, soly) / e)
  return 0

if __name__ == "__main__":
  N = 1696852658826990842058316561963467335977986730245296081842693913454799128341723605666024757923000936875008280288574503060506225324560725525210728761064310034604441130912702077320696660565727540525259413564999213382434231194132697630244074950529107794905761549606578049632101483460345878198682237227139704889943489709170676301481918176902970896183163611197618458670928730764124354693594769219086662173889094843054787693685403229558143793832013288487194871165461567
  e = 814161885590044357190593282132583612817366020133424034468187008267919006610450334193936389251944312061685926620628676079561886595567219325737685515818965422518820810326234612624290774570873983198113409686391355443155606621049101005048872030700143084978689888823664771959905075795440800042648923901406744546140059930315752131296763893979780940230041254506456283030727953969468933552050776243515721233426119581636614777596169466339421956338478341355508343072697451
  print boneh_durfee(N, e)
```


[solve.py](solve.py)
```python
from scryptos import *

c = (1569733526826523065259704222721381245770313117205865099913421859731162526943498524936251685846967970606251353344665893442015804015671457823645874503670136308040791285744658847419176471348768113798503897694020110157476679833746227801224812046930570487233225157924912272791212802495997329083436189937249314855532400635293522270501567950040825794060896420481676398789310029592608176167251882124182145471818654414925639589921023176070657483148482403065241178276749773L, 139537660044872985880471632333334179976891152860359271230202507995985566816703080930428310461057387079799847266510420206696052591677854190150642820963140050439023069266243433278700748622126726137374130247097863526461696642750021196138340072411724739383716017406022211953417323065831672315854246554523225039827L)
pubkey = (1696852658826990842058316561963467335977986730245296081842693913454799128341723605666024757923000936875008280288574503060506225324560725525210728761064310034604441130912702077320696660565727540525259413564999213382434231194132697630244074950529107794905761549606578049632101483460345878198682237227139704889943489709170676301481918176902970896183163611197618458670928730764124354693594769219086662173889094843054787693685403229558143793832013288487194871165461567L, 814161885590044357190593282132583612817366020133424034468187008267919006610450334193936389251944312061685926620628676079561886595567219325737685515818965422518820810326234612624290774570873983198113409686391355443155606621049101005048872030700143084978689888823664771959905075795440800042648923901406744546140059930315752131296763893979780940230041254506456283030727953969468933552050776243515721233426119581636614777596169466339421956338478341355508343072697451L, 171012227587318507773834753911094468358648971527111097308935888531930900156798659257578479378777764146070352809723708236353390208094909385240006920137781562826981091183813955039359863361624869703055918575613667858215532572602435432258750639197322091887713402631456113333645709142822182724397962837201266977523L, 96969753191136466007366303619618019752521508403657426306543836447071659732926802256183021740376016065813234292694535879838415771865207311953800116203362150588941093508091412441933752168889516206420588410478242229762908362637083786338280959547015086176046206126019992386890758970740552952647510652431386064722L)

def main():
  n, e, a, g = pubkey
  c1, c2 = c
  d = 100556095937036905102538523179832446199526507742826168666218687736467897968451
  print pow(pow(2, e, n), d, n) == 2
  k = pow(c1, d, n)
  K = pow(g, k, a)
  print long_to_bytes((c2 * modinv(K, a)) % a)


if __name__ == '__main__':
  main()
```

```
Mon Sep 11 02:35:17 JST 2017 ~/Downloads/gracias 100%
> sage boneh_durfee.sage
96495049525709737646237784043989949922984947124527440290534285859764556084120 -213794805403874041582980303133544562301251683427479249295482840787267897208520588155734823678603567894373005996132986582548573443974113270682593599108980563590099083306897483648896074029256401888893482421279582317475577712512180341863233891668553849278499685187092538645494649114612591891435800686745132070463
100556095937036905102538523179832446199526507742826168666218687736467897968451
Mon Sep 11 02:35:22 JST 2017 ~/Downloads/gracias 100%
> python solve.py 
True
ASIS{Wiener_at7ack_iN_mUlt1_Prim3_RSA_iZ_f34sible_t0O!}
```

Flag: `ASIS{Wiener_at7ack_iN_mUlt1_Prim3_RSA_iZ_f34sible_t0O!}`

# References
[1] M.Wiener. 1990. *Cryptanalysis of short RSA secret exponents*.
[2] D.Boneh and G.Durfee. 1999. *Cryptanalysis of RSA with Private Key $d$ Less Than $N^{0.292}$*.
