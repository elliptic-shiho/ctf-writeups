BabyPinhole (Crypto)
===================================

Challenge file(s) is here: [baby-pinhole.7z-876abd71fe4430fa31f03a890d25e66a8f132cc6e44f81ff7544f8eac98d5fc3](baby-pinhole.7z-876abd71fe4430fa31f03a890d25e66a8f132cc6e44f81ff7544f8eac98d5fc3)

## Writeup
In this challenge, We have a decryption oracle for Paillier Cryptosystem [1]. but we can see **only 1-bit of middle of plaintext**.
Indeed, we can see following code in given source code:

```python
  while True:
    line = fin.readline()[:4+cbits//4]
    ciphertext = int(line, 16) # Note: input is HEX
    m = decrypt(ciphertext, sk1, sk2, n, n2)
    fout.write(str((m >> b) & 1) + "\n")
    fout.flush()
```

where `b = mbits / 2` and `mbits = 1024`. 

Okay, we consider to use this "pinhole" for reveal all plaintext bits.
First, We use *binary borrow* for observe lowest-bits.
We use this fact: "If 0b10000....0 minus one then all bits are one" (For example, we can observe `0b10001 - 0b10 = 0b1000 - 0b1 = 0b111` easily).

We denote $i$-th bit of $x$ by $x[i]$, ciphertext by $c$ and plaintext corresponding to $c$ by $m$ (Note : In this challenge, we have only $c$ and Paillier cryptosystem public key).

---

**Algorithm 1:**

1. Set $m' \leftarrow m$.
2. Observe $b$-th bit and set $m'[b] \leftarrow 0$ using homomorphic property of Paillier Cryptosystem.
3. Set $i = 1$.
4. Compute $Enc(x) \leftarrow Enc(m' - 2^{b-i})$ and Observe $x[b]$.
5. If $x[b] = 1$ then guess $m[b-i] = 1$ and $m'[b - i] \leftarrow 0$. Otherwise, guess $m[b-i] = 0$.
6. Set $i \leftarrow i + 1$.
7. If $b - i \geq 0$ then Back to Step 3. Otherwise, Terminate.

---

Use this algorithm, we have $m[0], \ldots, m[b]$, and $m'$. Moreover, we have $m'[0] = m'[1] = \cdots = m'[b] = 0$.

We want to know $m[b+1], \cdots, m[mbits]$. In this case, we use this lemma:

---

> **Lemma 1 ([2] p.7 Lemma 1):**
> 
> Let $N$ be a random $n$-bit RSA modulus, $y\in\mathbb{Z}\_N^{ \* }$, $c$ an even element of $\mathbb{Z} _ N$ and $g$ an element in $\mathcal{B}$. Then, denoting $z = 2^{-1}\mod N$,
> $$
> (g^cy^N)^z = g^\frac{c}{2}y'^N\mod N^2
> $$
> for some $y'\in\mathbb{Z}\_N^{ \* }$.

---

**Proof**: See [2] p.7 Lemma 1.

---

We construct following algorithm using $m'$ and Lemma 1:

---

**Algorithm 2:**

1. Set $i \leftarrow 1$.
2. Use Lemma 1, "Bit-Shift" $m'$. i.e. $m' \leftarrow m'^z \mod N^2$ where $z = 2^{-1}\mod N$.
3. Observe $m'[b]$ and set $m[b + i] \leftarrow m'[b]$.
4. Set $i \leftarrow i + 1$.
5. If $b + i < mbits$ then back to Step 2. Otherwise, Terminate.

---

Now we have $m[0], \ldots, m[b]$ and $m[b + 1], \ldots, m[mbits]$. So we can calculate flag!

[solve.py](solve.py)

```python
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
```

```
Sun Sep  3 01:42:52 JST 2017 ~/Downloads/baby-pinhole/orig 100%
> python solve.py 
[+] Leak lower-side...
[+] i = 513 
[+] Done.
[+] Leak upper-side...
[+] l = 510 
[+] Done.
TWCTF{ccb71c01f350cf0bc844e87d161f84b9b479b439}
232d73f8df61b0d993648fb442a91458ab9b88f6d77414d7fb766441b308fc3d85a26271444afd1437e69dd884377e5b56872d2fa4f784f18ec8bea72a43b23d8343fc9f46bbaeaeee573908afe1a6c3d293deda3d24e29f16202dffaf47f8e26ad7296f8264782136a40237c5792c5d078b388e04c6dc34fe0f2001e4d11000
```

Flag: `TWCTF{ccb71c01f350cf0bc844e87d161f84b9b479b439}`

## References
[1] Pascal Paillier. 1999. [*Public-Key Cryptosystems Based on Composite Degree Residuosity Classes*](https://link.springer.com/chapter/10.1007/3-540-48910-X_16)
[2] Dario Catalano, Rosario Gennaro and Nick Howgrave-Graham. 2001. [*The Bit Security of Paillier’s Encryption Scheme and its Applications*](https://link.springer.com/chapter/10.1007/3-540-44987-6_15).
