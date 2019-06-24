from sage.all import *
import itertools

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      if BB[ii,jj] == 0:
        a += ' ' 
      elif BB[ii,jj] == 1:
        a += '1'
      else:
        a += 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    print a

RRo = RR
RR = RealField(4096)
### YQLAC5DCJR57PYVUBQ4PXMH47IO5IETPUI7EDFUR7JWTNIHNTEAA====
x1  = RR("1.1523683371067123349656382131317521483572497142298337343260320993424398216178709236539938709110452634807784715262138995755342817993685509357712770340229750497263666806472679182677446813307514752989998709851733038910120353481027193255623282195303273679305610332766192092795306395336470967297906833935279894252616095703554546236457175492266398641941143086328952666283556115648410157333967927769488481273065721795561095918412150012421542250430997607592179132116258276065709318979773496974099022073099502")
y1 = RR("785701515561796928936294283336388546237.76198513405518052745128365573717412032431857838308609736805706641397119523106356273567386950051527641654288830128296680222705035406388428925103738738880089143710871996591308092142413410462799119104211628940723436162204485785271731710389128475080795001712865906678930222773707391836436633044651508284465733937444061563710525469368853621204597936989528982407800283786961077055617880725523871857845636253399548919606316914163179078711224606329197843485595974741371")
x2 = RR("1.9788681780306785745156989575455836997628595976195551967457035744531437329843012720941996363829938236558497904581798051243772329161582321879360660717582006051876040760971474767304024929464189882580767726915416095996327004491382227802281015240630789402853666506660443130218440489869153208264839264361191403711079005451092879896580775565516043311227204323256062060636302634873383609155444647328923902763251787493352130750547736994397271532409511747826853294706030568422435743878427864639456394057079862")
y2 = RR("2800780348429075081385422902824585625983.9718926020338437020050029020823624261231050171252891745710847437841708047345973255413221720133608100491056729443460854596666227899379299259825534978829077355033242582097430638128055602755533412987750066891446221229772669000110400372735908811477249717115822683690394865669693744895646210678759522710089141476167873126736495304431404870994048252140031420103893970444269041047766479070093625398030832513319343802074702488746482787664792595305879288767850463661089")
x3 = RR("1.1670934471504288389286699040597145166483859406607818722290870762005830879563520820444583769798094486100304620471097338034767093997643829487907758687772023239707260734472162763033737562629817749955408516815773124364141649872257509424163110108092639724918657438264373626094282680210615254900773574389965613279176813046462954963793860609569983518700922195346251798337750431074023542866696934221388642181566976684008412054371564255018302161417018138623829465811427229918451188130718662890832098982222125")
y3 = RR("802707509246523948167720804790529379981.65200627237446296237041791338739143005671462858728927775450751967070365345683907995507470923968298285874283562502374335568470739347447815600897365404905225832360747786484781020758259038661605921813716501226097097691811131220851765551439580788138654797589031081162956293274990820524344279481532872541315373551218364140305887507030965924866440561850791854542913794394621606495561223518482006711916111610968280784647353332520349002673666762620118053665941771430482")

L = 10^200
L2 = 1
K = 1

x1, x2, x3 = map(lambda x: x.exact_rational(), [x1, x2, x3])
y1, y2, y3 = map(lambda x: x.exact_rational(), [y1, y2, y3])

xs = [x1, x2, x3]
ys = [y1, y2, y3]

# x1.denominator()
# x1.numerator()

M = []

for i in xrange(6):
  t = []
  t += [0] * i + [L2] + [0] * (5 - i)
  if i < 5:
    for j in xrange(3):
      t += [L*xs[j]^(4-i)]
    print "i: deg {}".format(4-i)
  elif i == 5:
    for j in xrange(3):
      t += [L*-ys[j]]
    print "i: y"
  M += [t]

print

M = Matrix(QQ, M) * K
matrix_overview(M)

print

B = M.LLL()

matrix_overview(B)

XX = var("x")

for x in B:
  x /= K
  a, b = x[:6], x[6:]
  x = vector(map(RR, list(a / L2) + list(b / L)))
  A, B, C, D, E, k = a
  if k == 1:
    print "Norm: {}".format(RR(x.norm()).str(no_sci=2, skip_zeroes=True).rstrip("."))
    print "(" + ', '.join(map(lambda t: t.str(no_sci=2, skip_zeroes=True).rstrip("."), x)) + ")"
    # f = A * XX^4 + B * XX^3 + C * XX^2 + D * XX^1 + E
    print E