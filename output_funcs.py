import numpy as np

# translate an arbitrary sequence of 'integrified bases' into ACTG. Notice the input is an integer, not a string of integers.
def translate_int_seq(int_seq):
    strint_seq=''.join(list(str(int_seq)))
    trans_dict={'1':'A','2':'C','3':'T','4':'G'}
    return(''.join(list(map((lambda x: trans_dict.get(x)),list(strint_seq)))))

