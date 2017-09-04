The Worst (Crypto/Misc)
===================================

Challenge file(s) is here: [the_worst.7z-272c73136b20d9433378ecd8d7ca2d424e1aa0228aa3da7a42347b1a5c80eba3](the_worst.7z-272c73136b20d9433378ecd8d7ca2d424e1aa0228aa3da7a42347b1a5c80eba3)

## Writeup
We have encryption program and execute log of that. The Encryption program using Low-exponent ($e = 3$!) RSA but the program are padded message enough.

We focused on `rand` function. The program doesn't call `srand`, so we can predict `rand()` deterministically. (Indeed, `sprintf` has a bug. `flag` variable is `signed char`, so if `rand() % 256 >= 128` then `"%02x"` is ignored and sprint `"ffffffXX"`. )

Now we can predict all padding bytes. so we know following equation:
$$
(2^\beta x + padding)^3 - c \equiv 0 \pmod N
$$
where $\beta$ is padding bit length, $c$ is known ciphertext, and $N$ is RSA modulus. 

However, we don't know flag length. Equivalently, we can't know padding length. but we have $2^\beta < n$ and flag length are less than 32-characters.
therefore, we can bruteforce $\beta$, and solve equation using Coppersmith's Method.

[solve.sage](solve.sage)

```python
from sage.all import *

c = 107797257717677608972937851535768450709511610558576880946974945352552014654578789957451627927062856875276247258537360204856065541745711767779202724883839942699975062610942686285595750185802199113561904929512045813037669798287030278015252176567317799073700505323362947735036768009089368962950593202307632448303
n = 147391055337422035628254503121086681936558549846700068002467521885608242008107714483022482590086265945147038822970959474588506763687846415276055867381614762806814336658075914791763287728087033132858626992976930135645795181618355388123999915086700928323326053433449469842426405802845975581734834128440695492649
e = 3

# PR := (Z/nZ)[x]
ZmodnZ = Zmod(n)
PR = PolynomialRing(ZmodnZ, 'x')
x = PR.gen()

# Padding (fixed param)
K = '67ff697351ff4aff29ffffffffffff467cff54ff1bffffff765a2e6333ffffff66320dff3158ff5a255d051758ff5effffffffffffff54110eff7441213dffff70ff3eff41ffff673e017effffff6bffff385c2affff3bff32ff3c54ff18ff5c021aff43ffffff3aff29ffff053c7cff75ffff61ffff5cffffff0f'

for i in xrange(1, 27):
  k = ZZ(K[:-2*i], 16)
  pol = (2^(len(hex(k))*4) * x + k)^e - c
  pol = pol.monic()
  roots = pol.small_roots(X=2^256+1)
  if len(roots) > 0:
    print '[+] Found Root! i = %d' % i
    print '[+] roots = %r' % roots
    x0 = roots[0]
    print '[+] long_to_bytes(roots[0]) = %r' % hex(ZZ(x0)).decode('hex')
    break
```

```
Sun Sep  3 15:21:12 JST 2017 ~/Downloads/the_worst 100%
> sage solve.sage
[+] Found Root! i = 23
[+] roots = [8882130418134271787060459933578365880192674082998543906681274137866]
[+] long_to_bytes(roots[0]) = 'TWCTF{wmPib5vUl5bly4TlJlqf}\n'
```

Flag: `TWCTF{wmPib5vUl5bly4TlJlqf}`
