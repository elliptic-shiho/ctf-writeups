from scryptos import *
import gmpy

d = int(open("enc2.raw", "rb").read().encode("hex"), 16)

d = gmpy.root(d, 65537)
print d

open("out_e_root.txt", "w").write(repr(d))
