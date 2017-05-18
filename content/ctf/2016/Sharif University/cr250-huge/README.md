Huge (Crypto 250pts)
==================

Challenge file(s) is here: [huge2.tar.gz](huge2.tar.gz)

and here: [before_replace/huge.tar.gz](before_replace/huge.tar.gz) (Challenge was replaced once. This is challenge file before replace.)

## Writeup
This is crazy and huge RSA Modulus (we were waiting 3 minutes to load rsa modulus!).

My teammate(@\_193s) said: "This RSA Modulus has very long bits. So I think it is $m^e$ < n."

we tried to calculate $e$-th root... Yeah, we got flag!

[solve.py](solve.py)

```python
from scryptos import *
import gmpy

d = int(open("enc2.raw", "rb").read().encode("hex"), 16)

d = gmpy.root(d, 65537)
print d

open("out_e_root.txt", "w").write(repr(d))
```

Flag: `d604a27c3abcbe3e077b98`

## References
