#!/usr/bin/env python
# coding: utf-8

# # HapBamG: PySam/Bamnostic: Navigating through a SAM/BAM file

# In[1]:


"""
# Install a pip package in the current Jupyter kernel
import sys
!{sys.executable} -m pip install bamnostic
"""


# In[2]:


import bamnostic as bs 
import sys
import os
import pysam
import numpy as np
import pandas as pd
import re
from collections import defaultdict 


# In[3]:


def startup(xinput):
    # All .bam files in directory
    
    bamlist = []
    for x in os.listdir():
        if x.endswith(".bam"):
            bamlist.append(x)
            
    # Error Messages
    erms = 'INPUT FORMAT:|Maternal Phased 1.bam||Maternal Phased 2.bam||Paternal Phased 1.bam||Paternal Phased 2|.bam |-g| |-flag|'
    helpc = '-g\t generate a .csv graph comparison\n-mm\t calculate mismatches frequency\n-ma\t calculate match frequency'
    
    uput = xinput.split()
    
    if len(uput) < 5:
        print(erms)
        print(helpc)
        startup(input())
        
    if len(uput) > 6:
        print('Too many arguments:')
        print(erms)
        print(helpc)
        startup(input())
    
        
    for i in range(0, 3):
        if uput[i][-4:] != '.bam':
            print("Invalid file type:")
            print(erms) 
            print(helpc)
            startup(input())
        elif uput[i] not in bamlist:
            print("Invalid .bam files:")
            print(erms)
            startup(input())
            
    else:
        flags = []
    
        for i in range(len(uput)):
            p1mat = uput[0]
            p2mat = uput[1]
            p1pat = uput[2]
            p2pat = uput[3]
        
        for i in range(3, len(uput)):
            
            flag = uput[i]
            
            if flag == '-g':
                flags.append('g')
            elif flag == '-mm':
                flags.append('X')
            elif flag == '-ma':
                flags.append('M')
                
        if len(flags) < 1:
            print("Invalid flags:")
            print(erms)
            startup(input())
            
        if len(flags) > 2:
            print("Too many flags:")
            print(erms)
            startup(input())


        print('Generating fastas for:', p1mat, p2mat, p1pat, p2pat,'...')
        return([p1mat, p2mat, p1pat, p2pat], flags)


# In[4]:


user_in = input()
startup_list = startup(user_in)


# In[6]:


# Sorted BAM files

p1mat = startup_list[0][0]
p2mat = startup_list[0][1]
p1pat = startup_list[0][2]
p2pat = startup_list[0][3]

hap1mat = bs.AlignmentFile(p1mat, 'rb')
hap2mat = bs.AlignmentFile(p2mat, 'rb')
hap1pat = bs.AlignmentFile(p1pat, 'rb')
hap2pat = bs.AlignmentFile(p2pat, 'rb')

# Read Count

p1mrc = int(pysam.view('-c', p1mat))
p2mrc = int(pysam.view('-c', p2mat))
p1prc = int(pysam.view('-c', p1pat))
p2prc = int(pysam.view('-c', p2pat))

print(p1mrc, p2mrc, p1prc, p2prc)

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


# In[ ]:


def cigsort(cigstring):
    matches = re.findall(r'(\d+)([=SMIDHX]{1})', cigstring)
    cigar_op = ([{'type':m[1], 'length':int(m[0])} for m in matches])
    return cigar_op


# In[ ]:


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


# In[ ]:


def matcount(rlist, flag):
    oper = []
    for i in range(len(flag)):
        if flag[i] != 'g':
            oper.append(flag[i])
    r_match = []
    for i in range(len(rlist)):
        rx = sum(cigturn(cigsort(rlist[i][1]), oper[0]))
        rsclip = sum(cigturn(cigsort(rlist[i][1]),'S'))
        r_match.append([rlist[i][0], rx/(rlist[i][2] - rsclip), rx, rlist[i][2]])
        
    return r_match


# In[ ]:


def pmmatch(nmatch1, nmatch2):
    nmatch = []
    for i in range(len(nmatch1)):
        for j in range(len(nmatch2)):
            if nmatch1[i][0] == nmatch2[j][0]:
                nmatch.append([nmatch1[i][0], nmatch1[i][1], nmatch2[j][1], nmatch1[i][2], nmatch2[j][2], nmatch1[i][3], nmatch2[j][3] ])
    return nmatch 


# In[ ]:


# Counting Matches Between Reads

# Maternal Reads

m1_match = matcount(rm1, startup_list[1])
m2_match = matcount(rm2, startup_list[1])

    
# Maternal Graph Parameters

m_match = pmmatch(m1_match , m2_match)

mp1 = []
mp2 = []
mat_label = []

            
for i in range(len(m_match)):
    mat_label.append(m_match[i][0])
    mp1.append(m_match[i][1])
    mp2.append(m_match[i][2])
    


# In[ ]:


# Counting Matches Between Reads

# Paternal Reads

p1_match = matcount(rp1, startup_list[1])
p2_match = matcount(rp2, startup_list[1])
    
# Paternal Graph Parameters

p_match = pmmatch(p1_match , p2_match)

pp1 = []
pp2 = []
pat_label = []

for i in range(len(p_match)):
    pat_label.append(p_match[i][0])
    pp1.append(p_match[i][1])
    pp2.append(p_match[i][2])


# In[ ]:


# Read Graphs Comparing Maternal and Paternal

mp_label = []
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


# In[ ]:


def graphn(label, p1, p2, gname):
    list_dict = {'Read Name': label, 'Phased 1': p1, 'Phased 2': p2} 
    df = pd.DataFrame(list_dict) 
    df.to_csv(gname, index=False) 


# In[ ]:


def cgraphn(label, p1, p2, pp1, pp2, gname):
    list_dict = {'Read Name': label, 'Maternal Phased 1': p1, 'Maternal Phased 2': p2, 'Paternal Phased 1': pp1, 'Paternal Phased 2': pp2}
    df = pd.DataFrame(list_dict) 
    df.to_csv(gname, index=False) 


# In[ ]:


if 'g' in startup_list[1]:
    graphn(mat_label, mp1, mp2, 'Maternalp1p2.csv')
    graphn(pat_label, pp1, pp2, 'Paternalp1p2.csv')
    cgraphn(mp_label, matp1, matp2, patp1, patp2,'MatPat.csv' )


# In[ ]:


hapdup1_pat_hapdup2_mat = []
hapdup1_mat_hapdup2_pat = []
patp1_inval = []
patp2_inval = []
matp1_inval = []
matp2_inval = []

for i in range(len(matp1)):
    # Invalid Cases
    if patp1[i] > 0.1:
        patp1_inval.append([mp_label[i]])
    elif patp2[i] > 0.1:
        patp2_inval.append([mp_label[i]])
    elif matp1[i] > 0.1:
        matp1_inval.append([mp_label[i]])
    elif matp2[i] > 0.1:
        matp2_inval.append([mp_label[i]])
    elif patp1[i] > patp2[i] and matp1[i] > matp2[i]:
        patp1_inval.append([mp_label[i]])
        patp2_inval.append([mp_label[i]])
        matp1_inval.append([mp_label[i]])
        matp2_inval.append([mp_label[i]])
    elif patp1[i] < patp2[i] and matp1[i] < matp2[i]:
        patp1_inval.append([mp_label[i]])
        patp2_inval.append([mp_label[i]])
        matp1_inval.append([mp_label[i]])
        matp2_inval.append([mp_label[i]])     
    # Valid Cases
    elif patp1[i] < patp2[i] and matp1[i] > matp2[i]:
        hapdup1_pat_hapdup2_mat.append([mp_label[i]])
    elif patp2[i] < patp1[i] and matp2[i] > matp1[i]:
        hapdup1_mat_hapdup2_pat.append([mp_label[i]])


# In[ ]:


def hapfast(hapt, hapm, fname):
    ht = bs.AlignmentFile(hapt, 'rb')
    pj = int(pysam.view('-c', hapt))
    for i in range(pj):
        hapre = next(ht)
        for j in range(len(hapm)):
            if hapre.read_name == hapm[j][0]:
                hapm[j].append(hapre.query_sequence)

    ofile = open(fname+".txt", "w")
    for i in range(len(hapm)): 
        ofile.write(">" + hapm[i][0] + "\n" + hapm[i][1] + "\n")
    ofile.close()


# In[ ]:


def invalfast(mp1, mp2, pp1, pp2, fname):
    p1mat = startup_list[0][0]
    p2mat = startup_list[0][1]
    p1pat = startup_list[0][2]
    p2pat = startup_list[0][3]
    
    hap1mat = bs.AlignmentFile(p1mat, 'rb')
    hap2mat = bs.AlignmentFile(p2mat, 'rb')
    hap1pat = bs.AlignmentFile(p1pat, 'rb')
    hap2pat = bs.AlignmentFile(p2pat, 'rb')
    
    p1mrc = int(pysam.view('-c', p1mat))
    for i in range(p1mrc):
        hap1matreads = next(hap1mat)
        for j in range(len(mp1)):
            if hap1matreads.read_name == mp1[j][0]:
                mp1[j].append(hap1matreads.query_sequence)

    p2mrc = int(pysam.view('-c', p2mat))
    for i in range(p2mrc):
        hap2matreads = next(hap2mat)
        for j in range(len(mp2)):
            if hap2matreads.read_name == mp2[j][0]:
                mp2[j].append(hap2matreads.query_sequence)

    p1prc = int(pysam.view('-c', p1pat))
    for i in range(p1prc):
        hap1patreads = next(hap1pat)
        for j in range(len(pp1)):
            if hap1patreads.read_name == pp1[j][0]:
                pp1[j].append(hap1patreads.query_sequence)
    
    p2prc = int(pysam.view('-c', p2pat))
    for i in range(p2prc):
        hap2patreads = next(hap2pat)
        for j in range(len(pp2)):
            if hap2patreads.read_name == pp2[j][0]:
                pp2[j].append(hap1patreads.query_sequence)
                
    inval_hap = mp1 + mp2+ pp1 + pp2 
    ofile = open(fname+".txt", "w")
    for i in range(len(inval_hap)): 
        ofile.write(">" + inval_hap[i][0] + "\n" + inval_hap[i][1] + "\n")
    ofile.close()    
    


# In[ ]:


hapfast(startup_list[0][3], hapdup1_pat_hapdup2_mat, 'Hapdup1_Paternal_Hapdup2_Maternal')
hapfast(startup_list[0][0], hapdup1_mat_hapdup2_pat, 'Hapdup1_Maternal_Hapdup2_Paternal')
invalfast(matp1_inval, matp2_inval, patp1_inval, patp2_inval, 'Inval_p1p2')
print("Completed.")


# In[ ]:





# In[ ]:




