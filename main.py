from scipy.stats import spearmanr
## Given an input string of nucleotides (DNA template strand), we would like to check if the nearby region of nucleotides after are C-rich or not. Output: binary string, 0 for not double C-rich (meaning not favorable for an r-loop) or 1 for double G-rich. We then use another function to determine if we have "enough" double-G-rich values in the neighborhood to call the region favorable for R-loop (1) or unfavorable (0). So it is a two-step process: First decide if each substring is double_C_rich, then decide if we have enough double_c_rich regions close enough together to start an R-loop.

## There are many helper functions that appear as part of the first-step function (which is double_c_rich), and several auxiliary functions that are used for determining the appropriate thresholds to use as parameters in double_c_rich.  

##double_c_rich(input_string, m, j, threshold1, threshold2)
# input_string: should be a template DNA sequence
# m: how often we would like to check for double_c_richness. m = 1 means we look at a substring starting at every nucleotide. m=5 checks only every fifth nucleotide. 
# j: how long a substring we would like to examine for double_c_richness. j = 10 means a moving window of size 10. 
#threshold1: this is the threshold in order to consider a substring as C-skewed. The calculation for C-skewness is (# of C's - # of G's)/(# of C's + # of G's), so a threshold1 of 0.3 means the C's to G's is at least 65% to 35%, or essentially a 2 to 1 ratio. 
#threshold2: this is the threshold in order to consider a substring as C-rich. This calculation is # of C's / # of nt. Note that C-skewness is not a purely stronger measure, since it doesn't include information about the proportion of A's and T's. A threshold of 0.6 means that at least 60% of the j nucleotides are C's. 

##The second-step function is findjonesoutoften(input_string, j).
# input-string: a binary string indicating the double_c_richness of a template DNA strand
# j: how many times we would need to see a 1 in a substring of length 10 before we consider it R-loop favorable (1) or unfavorable (0). 

#I could certainly combine those into a single function, but I prefer the two-step process at this point since we might still tune our thresholds as we get more data. 


sample_string = 'AACTAGCGTAAGGGCGGAGGGGGTAACAAATTAAA'
## for an input string of nucleotides, the function gratio gives the ratio of C's to the length of the string.
def cratio(input_string):
  number_of_cs=0
  for k in range(len(input_string)):
    if input_string[k]=='C':
      number_of_cs +=1
  return number_of_cs/len(input_string)


## apply gratio every mth nucleotide, looking ahead j nucleotides
def list_cratio(input_string, m, j):
  output_list=[]
  for k in range(len(input_string)//m):
    a = cratio(input_string[m*k:m*k+j])
    output_list.append(a)
  return output_list

def crich(ratio_list, threshold):
  output_string = ''
  for a in ratio_list:
    if a > threshold:
      output_string+='1'
    else:
      output_string+='0'
  return output_string

def crichstring(input_string, m, j, threshold):
  return crich(list_cratio(input_string, m, j), threshold)

def cskewratio(input_string):
  number_of_gs=0
  number_of_cs=0
  for k in range(len(input_string)):
    if input_string[k]=='G':
      number_of_gs +=1
    elif input_string[k]=='C':
      number_of_cs +=1
  if number_of_gs+number_of_cs == 0:
    a = 0
  else: a = (number_of_cs - number_of_gs)/(number_of_gs+number_of_cs)
  return a




def list_cskewratio(input_string, m, j):
  output_list=[]
  for k in range(len(input_string)//m):
    a = cskewratio(input_string[m*k:m*k+j])
    output_list.append(a)
  return output_list


def cskew(ratio_list, threshold):
  output_string = ''
  for a in ratio_list:
    if a > threshold:
      output_string+='1'
    else:
      output_string+='0'
  return output_string

def cskewstring(input_string, m, j, threshold):
  return cskew(list_cskewratio(input_string, m, j), threshold)


def double_c_rich(input_string, m, j, threshold1, threshold2):
  output_string=''
  str1=cskewstring(input_string, m, j, threshold1)
  str2=crichstring(input_string, m, j, threshold2)
  for k in range(len(input_string)//m):
    a=int(str1[k])
    b=int(str2[k])
    c=str(a*b)
    output_string+=c
  return output_string



##In pfc53, cluster 1 is 65-425bp, so 13-85 halftwists
##cluster 2 is 435 - 790bp, so 87 - 158
##cluster 3 is 805 - 1305 bp, so 161 - 261
##cluster 4 is the rest



pfc53long ='CCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCCATTCGCCATTCAGGCTACGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGGGATCCTGCAGTAATACGACTCACTATAGGGCGAATTGGAGCTCCACATCAAGCTTTTGAGCTTGCCTCTCTTGCAACGTGGCACTTTTGAGTTCATCTCTCCTGTAACATGGCACTTTTGAGCGTTGGCCCTGATCACAGAACCCTTCGAATCCTCCCCTTGTGCAGTTTGCACCCTCAGGATACCTCGGAATCCTCCGAGCCGTCTCTTCCCCTCCCCTCGCGCAACTTGGCATAACCAGAATCACAGCCCAACCGGAATCGCATTAAAACCCTCCGAACCTTTGGGCAGCGCGGCACCGGGGCTCACGCAGTGCCGCGGAACCCTCGGAACCCTGCCCTTGTGCAGCTTCGCACCCTCAGGATAGCTCGGAACCCTCTGAGCTTTCCCTTCCCTTTCCCTCACCGCAACTCAGCACAACCAAGGATCACGGCACAACCGGAATTGCTCCGAACCCTCGGGCAGCACTGCGGCGCAGTGCCGCGGAACCCTCTAAACCCTCCGAATCCCCTCCTTGTGCAGCTTTTCACCCTTAGCCTCTCTCGCAACTCGGCACAACCAGAATCGCAGCACAATTGGAGTTGTGCTAAAATCCTCTGAACGGGGGGGGGGGGGGGGGCAGCGCACACGCCAGCACCCTGGCTCGCGCAGGGCCGCGGAACCCACCGAACCTTCCCCCTGCCCCGCAAACCGCGGATCAGGTCGTGTAGAAATCCTCCGATCTCCACCCCCCCCCCAACTCTCCAGCAGCGTGGTGCCCTAGCGCGAACCCCTGCTCCCCTTGCGCAGCGTGGCACCCTCGTGCTCTAAATCGCCCGTAAACCCCCCCCCCCCAGTGCAGCAGTCGGCCCGCGCTGGACGCCCCCCTCCCTTCTCCTCTTGCTGACGCGGCATGGCGGCAGCGCGGCTTTCGGCTTTCTACCGGCCCTATCGGCCCTCGTGTAGTTCAGAACACTGGTGAGCAGTGGGGCACCCTAGTGCTGAGTTTTCTTGTAGCCCAGAAATCTTCACCCTAGCGCTGAATCTCGCGTAGGGGAACCTTTGAGCGGCCTGGAACTCTTGGTCGGAGCCCTCGAGCTGCCGGATCAGTGGAACCCTCACCTCGCGTAGAGCCCTTCCCTCCTGTGGAACCCTTCCTTTGCGGAATCCTCTAACGCGTGGAATCCTCCCCTTGCAACGGAATCCTACCCTCATCTGCAGAATCCTCCCCTCGCATTGCCGCGCTTCACGCGGGGAACCCTAGCTTACCTCAGTTCCCCAGGATGTCTGCGTGGTAACTGGCGGTGCTGTGCTTCTTGCTGCCCGCTAGAGCAAGGAGGGATAGATATAATGTTGAAGCCTCGGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTCTAGAGCATGCCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCGGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTCCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTGCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGAAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGCTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGAAATTGTAAACGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTT'

pfc53 = pfc53long[:1750]


testones = double_c_rich(pfc53long, 1, 10, 0.3, 0.62)

test = double_c_rich(pfc53, 5, 10, 0.3, 0.62)

testpresent = double_c_rich(pfc53long, 1, 10, 0.3, 0.6)

testten = double_c_rich(pfc53, 10, 10, 0.3, 0.62)

tester = crichstring(pfc53, 5, 10, 0.3)

pfc19 = 'CCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCCATTCGCCATTCAGGCTACGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGGGATCCTGCAGTAATACGACTCACTATAGGGCGAATTGGAGCTCCACATCAAGCTTATCGATACCGTCGGATTCCACGCATTTCTGCAGCCCCTGCGGCAGTGTAGGCATTGCGCAGTTTTAATAAAAAAGCACCACCACCACAGTAGGCAAACCAGATGACCATCGCAGGTCACAGGAAAATTAAAGGCTGCGGACTGTGCTACTGCCCCTTCTGATGCCCCCTCCTCTACACAGCAATCATTCAGCGTCCCTTAGTCACTCCGGACAGCGACAGGCCCCGCGGCCGCCATGCCCACCGCCTCCATGCCATGCCCACCGCCGCCATGCCTACCGCCGCCAAAGTCCACCACCGCCATGCCTACCCGCTGCCAATGCCCACCGCCGCCAATACCCACTGTCGCCGCCTTCCCCCTACCTCCCAGCCACTTCCTACGGACTCTCCCCGCGCCGCGACCACCAACACAACCCCCACCACTGTCACACCGACTCATCCCCCTGGTCCACTGCCATAGCCTCCTCGCCTCGGTCACTGCGACGAATTCCCCCCCCAGTCGCCCCACGTACCCTGCTCCACCACGCAGTGGTCACTATTATACACCTACCTGCGCTCAACACCCCCTAAATACCGATCACTTCACGTACCTTCGCCCCGCCACAATCACTCCAATATACCTACCTCCGCCTAAAATCCCTATGCACTGGTCCCCCCACGTACCCTCGCCACACGGAACTGCAATCACCCTGATGTACCCACCTCCACCCATGTCCCTTGCCCACTGCGGTTACCCCGCATGCTCCCAGTCACCACCGCCCTTCCCACCGCAGACACCCGCAATAGGACCTGTCGCGACACCACAGTTGGGGGCGGATGGGGGACGCGCCCCAATGCGAGCGGACAGGATACCATCGGGGCAGAACGGCACAACAGCAAGCCTCTGAACATTCCGGATCTGGTTCTCCAGAACAAAGGACTTTAGGGCCCAAATTCCGTTTATTCAGTACTCCAAGACCTCGAGGGGGGGCCCGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTCTAGAGCATGCCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCGGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTCCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTGCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGAAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGCTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGAAATTGTAAACGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGG'

string1='010011100101011111111000000'

paperexample = 'AGCCCCCGATCCCGA'

secondtest = double_c_rich(pfc19, 1, 10, 0.3, 0.6)
def findjonesoutoften(input_string,j):
  where_at=[]
  for k in range(len(input_string)-10):
    numb = 0
    for i in range(10):
      numb+=int(input_string[k+i])
    if numb >= j:
      where_at.append(k)
  return where_at



    
def findones(input_string, j):
  where_at=[]
  for k in range(len(input_string)):
    if input_string[k]==1:
      where_at.append(j*k)
  return where_at

def findjonesoutoffive(input_string,j):
  where_at=[]
  for k in range(len(input_string)-5):
    numb = 0
    for i in range(5):
      numb+=int(input_string[k+i])
    if numb >= j:
      where_at.append(k)
  return where_at

testagain = test = double_c_rich(pfc19, 1, 10, 0.3, 0.62)

## The values we used in our presentation were 6 out of 10 double-c-rich values, and thresholds 0.3 and 0.6. 


##Let's test spearmanr. 

def pfc19_list(threshold1, threshold2):
  p = double_c_rich(pfc19, 1, 10, threshold1, threshold2)
  listp = list(p)
  floated = [float(char) for char in listp]
  return(floated)

def pfc53_list(threshold1, threshold2):
  p = double_c_rich(pfc53long, 1, 10, threshold1, threshold2)
  listp = list(p)
  floated = [float(char) for char in listp]
  return(floated)


def accumlist53(j, t1, t2):
  where_at=[]
  p = pfc53_list(t1, t2)
  for i in range(10):
    where_at.append(0)
  for k in range(len(p)-10):
    numb = 0
    for i in range(10):
      numb+=int(p[k+i])
    if numb >= j:
      where_at.append(1)
    else: where_at.append(0)
  
  return where_at

##print(pfc19_list(0.3, 0.6))

x = [0,0,0,1,1,1,0,0,0]

x_corr = [0, 0.1, 0.1, 0.1, 0.2, 0.3, 0.7, 0.1, 0]
corr, p_value = spearmanr(x, x_corr)
#print(corr)
##print(p_value)
pfc53_full = open('pfc53_full_output.txt','r')
test_floats = open('testfloats.txt','r')

list53 = pfc53_full.readlines()
list53[0] ='0'
#print(list53)
float53 = [float(char) for char in list53]
#print(float53)

pfc19_full = open('pfc19_full_output.txt','r')
list19 = pfc19_full.readlines()
list19[0]='0'
float19 = [float(char) for char in list19]
#print(float19)

pfc530 = open('pfc53zeros.txt','r')

list530 = pfc530.readlines()
list530[0] ='0'
#print(list53)
float53 = [float(char) for char in list53]

float530 = [float(char) for char in list530[:3906]]
#corr53, p_value53 = spearmanr(float53, pfc53_list(0.3, 0.6))
#print(corr53)

pfc190 = open('pfc190.txt', 'r')

list190 = pfc190.readlines()
list190[0]='0'
float190 = [float(char) for char in list190[:3669]]

#corr19, p_value19 = spearmanr(float19, pfc19_list(0.3, 0.6))
#print(corr19)
#print(float53)

#print(float('34e-12'))
def tune53(h):
  best = [0,0,0,0]
  for t1 in range(0,1):
    for t2 in range(0,h,1):
      for j in range(0, 11):
        corr, p_value = spearmanr(float53, accumlist53(j, 0.1*t1, (1/h)*t2))
        if corr >= best[0]:
          best[0] = corr
          best[1] = j
          best[2] = 0.1*t1
          best[3] = (1/h)*t2
  return(best)

pd = pfc53_list(0,0.4)

import simplejson
f = open('output.txt', 'w')
simplejson.dump(pd, f)
f.close()

pd_bad = pfc53_list(0.3, 0.6)



pd_19 = pfc19_list(0, 0.8)

f19 = open('output19.txt','w')
simplejson.dump(pd_19, f19)
f19.close()




def tune53old(h):
  best = [0,0,0]
  for t1 in range(0, 1):
    for t2 in range(0,h,1):
      corr, p_value = spearmanr(float530, pfc53_list(0.1*t1, (1/h)*t2))
      if corr >= best[0]:
        best[0] = corr
        best[1] = 0.1*t1
        best[2] = (1/h)*t2
  return(best)

def tune19(h):
  best = [0,0,0]
  for t1 in range(0, 1):
    for t2 in range(0,h,1):
      corr, p_value = spearmanr(float190, pfc19_list(0.1*t1, (1/h)*t2))
      if corr >= best[0]:
        best[0] = corr
        best[1] = (1/h)*t1
        best[2] = (1/h)*t2
  return(best)

def countones(inputlist):
  bins = []
  for i in range(len(inputlist)-50):
    started=0
    for j in range(50):
      if inputlist[i+j]==1:
        started+=1
    bins.append(started)
  return(bins)

count53 = countones(pd)

f3 = open('outputcum.txt', 'w')
simplejson.dump(count53, f3)
f3.close()

count19 = countones(pd_19)

f4 = open('outputcum19.txt', 'w')
simplejson.dump(count19, f4)
f4.close()

pd19less = pfc19_list(0, 0.6)

count19less = countones(pd19less)

f5 = open('outputcum19less.txt','w')
simplejson.dump(count19less, f5)
f5.close()

count53bad = countones(pd_bad)

f2 = open('outputbadcum.txt', 'w')
simplejson.dump(count53bad, f2)
f2.close()


def expandlistbyj(inputlist, j):
  newlist = []
  for char in inputlist:
    for i in range(j):
      newlist.append(char)
  return(newlist)


def tune53ones(h):
  best = [0,0,0]
  for t1 in range(0,1):
    for t2 in range(0, h, 1):
      corr, p_value = spearmanr(float530[:3800], expandlistbyj(countones(pfc53_list(0.1*t1, (1/h)*t2)),10))
      if corr >= best[0]:
        best[0]=corr
        best[1]=0.1*t119
        best[2]=(1/h)*t2

  return(best)

float53short = float53[0:3800]
float53shorter = float53short[0::10]

float530short = float530[0:3800]
float530shorter = float530short[0::10]

def tune53onesbetter(h):
  best = [0,0,0]
  for t1 in range(0,1):
    for t2 in range(0, h, 1):
      corr, p_value = spearmanr(float530[:3856], countones(pfc53_list(0.1*t1, (1/h)*t2)))
      if corr >= best[0]:
        best[0]=corr
        best[1]=0.1*t1
        best[2]=(1/h)*t2

  return(best)

def tune19onesbetter(h):
  best = [0,0,0]
  for t1 in range(0,1):
    for t2 in range(0, h, 1):
      corr, p_value = spearmanr(float190[:3619], countones(pfc19_list(0.1*t1, (1/h)*t2)))
      if corr >= best[0]:
        best[0]=corr
        best[1]=0.1*t1
        best[2]=(1/h)*t2

  return(best)


