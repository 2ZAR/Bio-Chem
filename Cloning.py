#Practice#
#Cloning - Using Pydna#

import pydna.parsers as parsers
from pydna.amplify import pcr

C3 = pydna.parsers.parse('CEACAM3.ape',ds=True)
C3_F = pydna.parsers.parse_primers('C3_F.ape')
C3_R = pydna.parsers.parse_primers('C3_R.ape')

C3_pcr_prod = pcr(C3_F,C3_R,C3)

%matplotlib inline
from pydna.gel import weight_standard_sample, Gel
st = weight_standard_sample('1kb+_GeneRuler')
Gel([ st,[C3_pcr_prod]] , gel_len=16).run()

print(C3_pcr_prod.program())

print(C3_pcr_prod.dbd_program())

print(C3_pcr_prod.figure())

C5 = parsers.parse('CEACAM5.ape',ds=True)
C5_F = parsers.parse_primers('Cyno_C5-F.ape')
C5_R = parsers.parse_primers('Cyno_C5-R.ape')
C5_pcr_prod = pcr(C5_F,C5_R,C5)