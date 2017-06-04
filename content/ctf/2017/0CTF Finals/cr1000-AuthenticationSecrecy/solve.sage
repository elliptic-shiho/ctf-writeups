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
