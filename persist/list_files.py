# LIST FILES CREATING
# ***
# FOR ASTRONIRCAM
# ***

import glob

def toDig(s):
        return int(s[s.rfind('-')+1:s.rfind('.')])
     
def sort_files(list):
        d = dict()
        dig = []
        for el in list:
            d[toDig(el)]=el
            dig.append(toDig(el))
        dig.sort()
        k = 0
        for i in dig:
            list[k] = d[i]
            k = k + 1
        return list   

#l=sort_files(glob.glob('*-1-*.fits'))
#print l