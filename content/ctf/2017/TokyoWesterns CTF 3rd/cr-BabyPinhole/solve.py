from collections import namedtuple
from roputils import *
from scryptos import *
import hashlib
import random
import sys

PaillierPublicKey = namedtuple('PaillierPublicKey', ['n', 'n2', 'g'])
PaillierSecretKey = namedtuple('PaillierSecretKey', ['sk1', 'sk2'])

def encrypt(m, pubkey):
  assert isinstance(pubkey, PaillierPublicKey)
  r = random.randint(1, pubkey.n2)
  return (pow(pubkey.g, m, pubkey.n2) * pow(r, pubkey.n, pubkey.n2)) % pubkey.n2

def decrypt(c, pubkey, secretkey):
  def L(x, n):
    return (x - 1) / n
  return (L(pow(c, secretkey.sk1, pubkey.n2), pubkey.n) * secretkey.sk2) % pubkey.n

def add(c1, c2, pubkey):
  return (c1 * c2) % pubkey.n2

def scalarmult(c, k, pubkey):
  if k == 0:
    return 1
  if k < 0:
    c = modinv(c, pubkey.n2)
    k = -k
  return pow(c, k, pubkey.n2)

def bit_on(c, k, pubkey):
  return add(c, encrypt(1<<k, pubkey), pubkey)

def bit_off(c, k, pubkey):
  return add(c, scalarmult(encrypt(1<<k, pubkey), -1, pubkey), pubkey)

p = Proc(host='ppc2.chal.ctf.westerns.tokyo', port=38264)
# p = Proc(host='localhost', port=38264)#, display=True)

def oracle(c):
  c = hex(c)[2:].rstrip('L')
  if len(c) % 2 == 1:
    c = '0' + c
  p.writeline(c)
  return eval(p.readline())

pubkey = PaillierPublicKey(*map(lambda x: int(x, 16), open('publickey').read().split()))

c = int(open('ciphertext').read(), 16)

cur = c
mbits = 1024
cbits = mbits * 2
b = mbits // 2
b0 = oracle(cur)
bits = {b: b0}

if b0 == 1:
  cur = bit_off(cur, b, pubkey)

print '[+] Leak lower-side...'
i = 1
while b - i >= 0:
  t = bit_off(cur, b-i, pubkey)
  if oracle(t) == 1:
    bits[b-i] = 0
  else:
    bits[b-i] = 1
    cur = t
  i += 1
  print '\r[+] i = %d' % i,
  sys.stdout.flush()

print '\n[+] Done.'
known_bits = [bits[x] for x in xrange(b + 1)]
s = ''.join(map(str, known_bits))[::-1]

print '[+] Leak upper-side...'
z = modinv(2, pubkey.n)
ginv = modinv(pubkey.g, pubkey.n2)
res = []
while True:
  x = oracle(cur)
  res = [x] + res
  cur = pow(cur, z, pubkey.n2)
  if len(res) > 510:
    break
  print '\r[+] l = %d' % len(res),
  sys.stdout.flush()

print '\n[+] Done.'

data = (int(''.join(map(str, res))[:-1] + s, 2))
print 'TWCTF{%s}' % hashlib.sha1(str(data)).hexdigest()
print long_to_bytes(data).encode("hex")