equation (Crypto 2pts)
===================================

Challenge file(s) is here: [equation.zip](equation.zip)

## Writeup
We have partially masked private key and encrypted flag. so we need to "OCR" image first.

OCR-ed private key is here:

```
-----BEGIN RSA PRIVATE KEY-----
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    ~~~some~~~lines~~~
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Os9mhOQRdqW2cwVrnNI72DLcAXpXUJ1HGwJBANWiJcDUGxZpnERxVw7s0913WXNt
V4GqdxCzG0pG5EHThtoTRbyX0aqRP4U/hQ9tRoSoDmBn+3HPITsnbCy67VkCQBM4
xZPTtUKM6Xi+16VTUnFVs9E4rqwIQCDAxn9UuVMBXlX2Cl0xOGUF4C5hItrX2woF
7LVS5EizR63CyRcPovMCQQDVyNbcWD7N88MhZjujKuSrHJot7WcCaRmTGEIJ6TkU
8NWt9BVjR4jVkZ2EqNd0KZWdQPukeynPcLlDEkIXyaQx
-----END RSA PRIVATE KEY-----
```

We decode this masked key according to PEM structure then we can know three parameters $d\mod (p-1)$, $d\mod (q-1)$ and $q^{-1}\mod p$ where $d$ is RSA private exponent and $p$, $q$ is prime (but, we don't know $d$, $p$, and $q$).

For simple, we denote $d\mod (p-1)$ as $d\_p$, $d\mod (q-1)$ as $d\_q$ and $q^{-1}\mod p$ as $q\_{\mathrm{inv}}$.

Well, we have a part of RSA-CRT parameters. they are satisfied following equation:

$$
\begin{aligned}
\begin{eqnarray\*}
d&\equiv d\_p\mod (p-1)\\\\
\iff d &= d\_p + k\_p(p-1) & \exists k\_p\in\mathbb{Z}\tag{1}\\\\
\end{eqnarray\*}
\end{aligned}
$$

$$
\begin{aligned}
d&\equiv d\_q\mod (q-1)\\\\
\iff d &= d\_q + k\_q(q-1) & \exists k\_q\in\mathbb{Z}\\\\
\end{aligned}
$$

In addition, $d$ is satisfy $ed\equiv 1 \mod \phi(n) \iff ed - 1 = k\phi(n)$ where $k$ is some integer. so, we know $d\_p e\equiv 1\mod(p-1)$ and $d\_q e\equiv 1\mod (q-1)$. and we use $(1)$, we know $ed\_p - 1 = k\_{ep}(p-1)$ where $k\_{ep}$ is some integer.

We guess $e = 65537$, and $k\_{ep}$ is sufficiently small(equivalently, we guess $k\_{ep}$ is able to bruteforce.). so, we bruteforce $k\_{ep}$.

If we guessed correct $k\_{ep}$, we know $\displaystyle\frac{ed\_p - 1}{k\_{ep}}-1 = p$. furthermore, $p$ is a prime. thus we construct algorithm as follows:

1. Guess $k\_{ep}$ (we use a count-up variable)
2. If $ed\_p - 1\mod k\_{ep} \ne 0$, back to 1.
3. Compute $\displaystyle p' = \frac{ed\_p - 1}{k\_{ep}}-1$
4. If $p'$ is a prime, terminate. otherwise, back to 1.

Finally, We implement this algorithm and factor $n$ and decrypt a flag.

```python
from scryptos import *
import gmpy

dp = 0xD5A225C0D41B16699C4471570EECD3DD7759736D5781AA7710B31B4A46E441D386DA1345BC97D1AA913F853F850F6D4684A80E6067FB71CF213B276C2CBAED59
dq = 0x1338C593D3B5428CE978BED7A553527155B3D138AEAC084020C0C67F54B953015E55F60A5D31386505E02E6122DAD7DB0A05ECB552E448B347ADC2C9170FA2F3
iq = 0xD5C8D6DC583ECDF3C321663BA32AE4AB1C9A2DED6702691993184209E93914F0D5ADF415634788D5919D84A8D77429959D40FBA47B29CF70B943124217C9A431

e = 65537
c = int(open("flag.enc", "rb").read().encode("hex"), 16)

kp_p_1 = e*dp-1

cand = 1
while True:
  if kp_p_1%cand == 0:
    p = kp_p_1/cand + 1
    if gmpy.is_prime(p):
      break
  cand += 1


q = modinv(iq, p)

rsa = RSA(e, p*q, p=p, q=q)
print rsa
m = hex(rsa.decrypt(c), 2)[2:].decode("hex")
print m[m.find("0ctf"):].strip()
```

Execute log:

```
Mon Mar 14 09:45:16 JST 2016 ~/ctf/0ctf-2016/crypto3 100%
> python solve.py 
RSA Private Key(e=65537, n=161080154188292201430717335450301702574211890587423028785946588452513709903864566907797711402814280216429284407010865117658741411399738837015270166197792615276511302372234182990420185803542388458087342116253675425489502589540709488892694405415013333511961708962693793627275736479090319881934245022826824347203, p=12883429939639100479003058518523248493821688207697138417834631218638027564562306620214863988447681300666538212918572472128732943784711527013224777474072569, q=12502893634923161599824465146407069882228513776947707295476805997311776855879024002289593598657949783937041929668443115224477369136089557911464046118127387, d=12441639692099655517376308833932392670257420848582256919212988552216677594845086557017745931627109670194928630671056032860651983223301005431608062335676428430110171020554477490159485308455680772826276447201841772149055876380443034602731403064627486237285806604612267999273183028007861118868108999965277036321)
0ctf{Keep_ca1m_and_s01ve_the_RSA_Eeeequati0n!!!}
```

Flag: `0ctf{Keep_ca1m_and_s01ve_the_RSA_Eeeequati0n!!!}`

# References