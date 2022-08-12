#!/usr/bin/env python
# coding: utf-8

# # PySam/Bamnostic: Navigating through a SAM/BAM file

# In[87]:


"""
# Install a pip package in the current Jupyter kernel
import sys
!{sys.executable} -m pip install bamnostic
"""


# In[88]:


import bamnostic as bs 
import pysam
import numpy as np
import pandas as pd
import re
from collections import defaultdict 


# In[89]:


# Sorted BAM files

p1mat = '1mat_sort.bam'
p2mat = '2mat_sort.bam'
p1pat = '1pat_sort.bam'
p2pat = '2pat_sort.bam'

hap1mat = bs.AlignmentFile(p1mat, 'rb')
hap2mat = bs.AlignmentFile(p2mat, 'rb')
hap1pat = bs.AlignmentFile(p1pat, 'rb')
hap2pat = bs.AlignmentFile(p2pat, 'rb')

# Read Count

p1mrc = int(pysam.view('-c', p1mat))
p2mrc = int(pysam.view('-c', p2mat))
p1prc = int(pysam.view('-c', p1pat))
p2prc = int(pysam.view('-c', p2pat))

# Maternal Comparisons 

m1 = []
rm1 = []

for i in range(p1mrc):
    hap1matreads = next(hap1mat)
    
    a_seq = hap1matreads.query_sequence
    
    if a_seq is not None:
        q_len = len(a_seq)
    else:
        q_len = 0
    m1.append([hap1matreads.read_name, hap1matreads.cigarstring, q_len])

# Removes all reads that have an alignment length of "0"

for i in range(len(m1)):
    if m1[i][2] != 0:
        rm1.append(m1[i])

      
m2= []
rm2 = []

for i in range(p2mrc):
    
    hap2matreads = next(hap2mat)
    
    a_seq = hap2matreads.query_sequence
    
    if a_seq is not None:
        q_len = len(a_seq)
    else:
        q_len = 0    

    m2.append([hap2matreads.read_name, hap2matreads.cigarstring, q_len])
    
for i in range(len(m2)):
    if m2[i][2] != 0:
        rm2.append(m2[i])
        
# Paternal Comparisons 
p1 = []
rp1 = []

for i in range(p1prc):
    hap1patreads = next(hap1pat)
    
    a_seq = hap1patreads.query_sequence
    
    if a_seq is not None:
        q_len = len(a_seq)
    else:
        q_len = 0
        #print(bs.utils.flag_decode(hap1patreads.flag))
        
    p1.append([hap1patreads.read_name, hap1patreads.cigarstring, q_len])
    
for i in range(len(p1)):
    if p1[i][2] != 0:
        rp1.append(p1[i])

p2 = []
rp2 = []

for i in range(p2prc):
    hap2patreads = next(hap2pat)
    
    a_seq = hap2patreads.query_sequence
    if a_seq is not None:
        q_len = len(a_seq)
    else:
        q_len = 0
        #print(bs.utils.flag_decode(hap2patreads.flag))    
    p2.append([hap2patreads.read_name, hap2patreads.cigarstring, q_len])

for i in range(len(p2)):
    if p2[i][2] != 0:
        rp2.append(p2[i])


# In[90]:


def cigsort(cigstring):
    matches = re.findall(r'(\d+)([=SMIDHX]{1})', cigstring)
    cigar_op = ([{'type':m[1], 'length':int(m[0])} for m in matches])
    return cigar_op


# In[91]:


def cigturn(cigdict, operation):
    
    #Segment of the query sequence that does not appear in the alignment. 
    seg = []
    #Match (alignment column containing two letters). 
    match = []
    #Insertion (gap in the query sequence). 
    ins = []
    #Deletion (gap in the target sequence).
    dele = []
    #Segment of the query sequence that does not appear in the alignment. 
    nota = []
    #Alignment column containing a mismatch, i.e. two different letters.
    xma = []
    # Alignment column containing two identical letters.
    eq = []

    for i in cigdict:
        if i['type'] == 'S':
            seg.append(i['length'])
        elif i['type'] == 'M':
            match.append(i['length'])
        elif i['type'] == 'I':
            ins.append(i['length'])
        elif i['type'] == 'D':
            dele.append(i['length'])
        elif i['type'] == 'H':
            nota.append(i['length'])
        elif i['type'] == 'X':
            xma.append(i['length'])
        elif i['type'] == '=':
            eq.append(i['length'])
        

    if operation == 'S':
        return seg
    elif operation == 'M':
        return match
    elif operation == 'I':
        return ins
    elif operation == 'D':
        return dele
    elif operation == 'H':
        return nota
    elif operation == 'X':
        return xma
    elif operation == '=':
        return eq
    else:
        return ["Inv"]


# In[92]:


def matcount(rlist):
    r_match = []
    for i in range(len(rlist)):
        rx = sum(cigturn(cigsort(rlist[i][1]),'X'))
        rsclip = sum(cigturn(cigsort(rlist[i][1]),'S'))
        r_match.append([rlist[i][0], rx/(rlist[i][2] - rsclip)])
    return r_match


# In[93]:


def pmmatch(nmatch1, nmatch2):
    nmatch = []
    for i in range(len(nmatch1)):
        for j in range(len(nmatch2)):
            if nmatch1[i][0] == nmatch2[j][0]:
                nmatch.append([nmatch1[i][0], nmatch1[i][1], nmatch2[j][1]])
    return nmatch 


# In[94]:


# Counting Matches Between Reads

# Maternal Reads

m1_match = matcount(rm1)
m2_match = matcount(rm2)

# Maternal Graph Parameters

m_match = pmmatch(m1_match , m2_match)

mp1 = []
mp2 = []
mat_label = []
            
for i in range(len(m_match)):
    mat_label.append(m_match[i][0])
    mp1.append(m_match[i][1])
    mp2.append(m_match[i][2])


# In[95]:


# Counting Matches Between Reads

# Paternal Reads

p1_match = matcount(rp1)
p2_match = matcount(rp2)
    
# Paternal Graph Parameters

p_match = pmmatch(p1_match , p2_match)

pp1 = []
pp2 = []
pat_label = []

for i in range(len(p_match)):
    pat_label.append(p_match[i][0])
    pp1.append(p_match[i][1])
    pp2.append(p_match[i][2])


# In[96]:


# Read Graphs Comparing Maternal and Paternal

mp_label = []
mp_len = []
patp1 = []
patp2 = []
matp1 = []
matp2 = []

for j in range(len(mat_label)):
    for i in range(len(pat_label)):
        if pat_label[i] == mat_label[j]:
            mp_label.append(mat_label[j])
            patp1.append(pp1[i])
            patp2.append(pp2[i])
            matp1.append(mp1[j])
            matp2.append(mp2[j])


# In[97]:


# Maternal Graph
list_dict = {'Read Name': mat_label, 'Phased 1': mp1, 'Phased 2': mp2} 
df = pd.DataFrame(list_dict) 
df.to_csv('matp1p2.csv', index=False) 


# In[98]:


# Paternal Graph 
list_dict = {'Read Name': pat_label, 'Phased 1': pp1, 'Phased 2': pp2} 
df = pd.DataFrame(list_dict) 
df.to_csv('patp1p2.csv', index=False) 


# In[99]:


# Comparing Maternal Paternal Graph
list_dict = {'Read Name': mp_label, 'Maternal Phased 1': matp1, 'Maternal Phased 2': matp2, 'Paternal Phased 1': patp1, 'Paternal Phased 2': patp2}
df = pd.DataFrame(list_dict) 
df.to_csv('matpat.csv', index=False) 


# In[101]:


"""
p1pat = '1pat_sort.bam'
hap1pat = bs.AlignmentFile(p1pat, 'rb')
p1prc = int(pysam.view('-c', p1pat))

for i in range(p1prc):
    hap1patreads = next(hap1pat)
    if hap1patreads.read_name == '0_phaseblock_93':
        print("Mismatch",sum(cigturn(cigsort(hap1patreads.cigarstring),'X')))
        print("Soft clip",sum(cigturn(cigsort(hap1patreads.cigarstring),'S')))
        print("Hard clip",sum(cigturn(cigsort(hap1patreads.cigarstring),'H')))
        print("Q align length",len(hap1patreads.query_sequence))
        print(hap1patreads.flag)
        print(hap1patreads.cigarstring)
"""


# In[ ]:





# In[ ]:




