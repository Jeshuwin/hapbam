#!/usr/bin/env python
# coding: utf-8

# # PySam/Bamnostic: Navigating through a SAM/BAM file

# In[3]:


"""
# Install a pip package in the current Jupyter kernel
import sys
!{sys.executable} -m pip install bamnostic
"""


# In[4]:


import bamnostic as bs 
import pysam
import numpy as np
import pandas as pd
import re
from collections import defaultdict 


# In[15]:


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
        #print(bs.utils.flag_decode(hap1matreads.flag) )
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


# In[6]:


def cigsort(cigstring):
    matches = re.findall(r'(\d+)([=SMIDHX]{1})', cigstring)
    cigar_op = ([{'type':m[1], 'length':int(m[0])} for m in matches])
    return cigar_op


# In[7]:


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


# In[8]:


# Counting Matches Between Reads

# Maternal Reads

m1_match = []
for i in range(len(rm1)):
    m1mx = sum(cigturn(cigsort(rm1[i][1]),'X'))
    m1sclip = sum(cigturn(cigsort(rm1[i][1]),'S'))
    m1_match.append([rm1[i][0], m1mx/(rm1[i][2] - m1sclip)])
    


m2_match = []
for i in range(len(rm2)):
    m2mx = sum(cigturn(cigsort(rm2[i][1]),'X'))
    m2sclip = sum(cigturn(cigsort(rm2[i][1]),'S'))
    m2_match.append([rm2[i][0], m2mx/(rm2[i][2] - m2sclip)])
    
# Maternal Graph Parameters

m_match = []
mp1 = []
mp2 = []
mat_label = []

for i in range(len(m1_match)):
    for j in range(len(m2_match)):
        if m1_match[i][0] == m2_match[j][0]:
            m_match.append([m1_match[i][0], m1_match[i][1], m2_match[j][1]])
            
for i in range(len(m_match)):
    mat_label.append(m_match[i][0])
    mp1.append(m_match[i][1])
    mp2.append(m_match[i][2])
    


# In[9]:


# Counting Matches Between Reads

# Paternal Reads

p1_match = []
for i in range(len(rp1)):
    p1x = sum(cigturn(cigsort(rp1[i][1]),'X'))
    p1sclip = sum(cigturn(cigsort(rp1[i][1]),'S'))
    p1_match.append([rp1[i][0], p1x/(rp1[i][2] - p1sclip)])
    
p2_match = []

for i in range(len(rp2)):
    p2x = sum(cigturn(cigsort(rp2[i][1]),'X'))
    p2sclip = sum(cigturn(cigsort(rp2[i][1]),'S'))
    p2_match.append([rp2[i][0], p2x/(rp2[i][2] - p2sclip)])
    
# Paternal Graph Parameters

p_match = []
pp1 = []
pp2 = []
pat_label = []


for i in range(len(p1_match)):
    for j in range(len(p2_match)):
        if p1_match[i][0] == p2_match[j][0]:
            p_match.append([p1_match[i][0], p1_match[i][1], p2_match[j][1]])
            #p_match.append([p1_match[i][0], p1_match[i][1], p2_match[j][1]])

for i in range(len(p_match)):
    pat_label.append(p_match[i][0])
    pp1.append(p_match[i][1])
    pp2.append(p_match[i][2])
 


# In[10]:


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


# In[11]:


# Maternal Graph
list_dict = {'Read Name': mat_label, 'Phased 1': mp1, 'Phased 2': mp2} 
df = pd.DataFrame(list_dict) 
df.to_csv('matp1p2.csv', index=False) 


# In[12]:


# Paternal Graph 
list_dict = {'Read Name': pat_label, 'Phased 1': pp1, 'Phased 2': pp2} 
df = pd.DataFrame(list_dict) 
df.to_csv('patp1p2.csv', index=False) 


# In[13]:


# Comparing Maternal Paternal Graph
list_dict = {'Read Name': mp_label, 'Maternal Phased 1': matp1, 'Maternal Phased 2': matp2, 'Paternal Phased 1': patp1, 'Paternal Phased 2': patp2}
df = pd.DataFrame(list_dict) 
df.to_csv('matpat.csv', index=False) 


# In[14]:


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


# In[ ]:




