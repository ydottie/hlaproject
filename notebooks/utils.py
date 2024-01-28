def fix(s):
    firstcolon = s.find(":")
    
    # if there is only one field prediction
    if firstcolon == -1:
        star = s.find("*")
        s_new = s[star+1:]
        s_new.ljust(4, '0')
    else:
        s_new = s[firstcolon-2:firstcolon] + s[firstcolon+1:firstcolon+3]
        
    return s_new

def compute_resolution(gs_val,pre_val):
    gs_all = gs_val.split("/")
    pre_fixed = fix(pre_val)
    flag = False #false = inaccurate, turn to true = 2 digits accurate
    
    for val in gs_all:
        gs_fixed = fix(val)

        if (gs_fixed[0:2] == pre_fixed[0:2]):
            if (gs_fixed[2:4] == pre_fixed[2:4]):
                return 4
            flag = True # use flag, rather than "return 2", because there can be multiple val to compare
        
    return 2 if flag==True else 0

# formats val to exact X*XX:XX format
def reformat(gs_val):
    numbers = fix(gs_val)
    if (gs_val[0] == 'D'): 
        return gs_val[0:4]+'*'+numbers[0:2]+':'+numbers[2:4]
    else:
        return gs_val[0]+'*'+numbers[0:2]+':'+numbers[2:4]
    

# updates accuracy for paired allele comparison (2 pre - 2 gs)
def pair_allele_tally(ans1,ans2,ans3,ans4,zerodig,twodig,fourdig,gs_locus):
    
    # parallel comparison
    if (ans1+ans2 > ans3+ans4):
        
        # assign accuracy for allele 1
        if (ans1 == 0):
            if gs_locus != 'D':
                zerodig[0]  = zerodig[0] + 1 
            else:
                zerodig[1] = zerodig[1] + 1
        if (ans1 == 2):
            if gs_locus != 'D':
                twodig[0] = twodig[0] + 1 
            else:
                twodig[1] = twodig[1] + 1
        if (ans1 == 4):
            if gs_locus != 'D':
                fourdig[0] = fourdig[0] + 1 
            else:
                fourdig[1] = fourdig[1] + 1
        
        # assign accuracy for allele 2
        if (ans2 == 0):
            if gs_locus != 'D':
                zerodig[0] = zerodig[0] + 1 
            else:
                zerodig[1] = zerodig[1] + 1
        if (ans2 == 2):
            if gs_locus != 'D':
                twodig[0] = twodig[0] + 1 
            else:
                twodig[1] = twodig[1] + 1
        if (ans2 == 4):
            if gs_locus != 'D':
                fourdig[0] = fourdig[0] + 1 
            else:
                fourdig[1] = fourdig[1] + 1
                
    # crosswise comparison
    else:
        # assign accuracy for allele 1
        if (ans3 == 0):
            if gs_locus != 'D':
                zerodig[0] = zerodig[0] + 1 
            else:
                zerodig[1] = zerodig[1] + 1
        if (ans3 == 2):
            if gs_locus != 'D':
                twodig[0] = twodig[0] + 1 
            else:
                twodig[1] = twodig[1] + 1
        if (ans3 == 4):
            if gs_locus != 'D':
                fourdig[0] = fourdig[0] + 1 
            else:
                fourdig[1] = fourdig[1] + 1
                
        # assign accuracy for allele 2
        if (ans4 == 0):
            if gs_locus != 'D':
                zerodig[0] = zerodig[0] + 1 
            else:
                zerodig[1] = zerodig[1] + 1
        if (ans4 == 2):
            if gs_locus != 'D':
                twodig[0] = twodig[0] + 1 
            else:
                twodig[1] = twodig[1] + 1
        if (ans4 == 4):
            if gs_locus != 'D':
                fourdig[0] = fourdig[0] + 1 
            else:
                fourdig[1] = fourdig[1] + 1
    
    return zerodig,twodig,fourdig

# updates accuracy for single allele comparison (2 pre - 1 gs, or 1 gs - 2 pre)
def single_allele_tally(ans1,ans2,zerodig,twodig,fourdig,gs_locus):
    ans = max(ans1,ans2)
    
    if (ans == 0):
        if gs_locus != 'D':
            zerodig[0] = zerodig[0] + 1 
        else:
            zerodig[1] = zerodig[1] + 1
    if (ans == 2):
        if gs_locus != 'D':
            twodig[0] = twodig[0] + 1 
        else:
            twodig[1] = twodig[1] + 1
    if (ans == 4):
        if gs_locus != 'D':
            fourdig[0] = fourdig[0] + 1 
        else:
            fourdig[1] = fourdig[1] + 1
    return zerodig,twodig,fourdig

    
# DESCRIPTION
# for general accuracy calculation. Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Loci in GS, but not in PRE, are tallied in the "nocall" variable 
# REQUIREMENTS: 
# 1. gs accession numbers are under a column labeled "Run" 
# 2. pre accession numbers are under a column labeled "ERR" 
# 3. accession numbers/column titles are labeled identically between gs and pre csv's
# INPUTS:
# Pandas dataframes for gs and pre
# OUTPUTS:
# count number of zerodig, twodig, fourdig, and nocall. 
# Each is a 2-item list, the first being Class I and second being Class II

def compute_matches(pre,gs):

    # Initializing variables
    fourdig = [0,0]
    twodig = [0,0]
    zerodig = [0,0]
    nocall = [0,0]

    d5 = ['Run', 'A', 'B']
    d6 = ['Run', 'C']

    accession_numbers = gs["Run"].values.tolist()
    genes = gs.columns.values.tolist()

    # iterate through accesssion numbers to check accuracy one sample at a time
    for number in accession_numbers:
        
        # get the gs and pre alleles for the sample
        pre_row = pre.loc[pre['ERR'] == number]
        gs_row = gs.loc[gs['Run'] == number]

        # if we are working with d5 or d6, the monoallelic datasets
        if (gs.columns.tolist() == d5 or gs.columns.tolist() == d6):
            
            # use single_allele_tally to tally accuracy of only 1 allele
            for i in range(1,len(genes)):
                gs_val = gs_row[genes[i]].astype(str).values[0]
                pre_val1 = pre_row[genes[i]].astype(str).values[0]
                pre_val2 = pre_row[genes[i]+".1"].astype(str).values[0]
                gs_locus = gs_val[0] # stores first digit of gs_locus to differentiate class I vs II
                
                # addresses edge case in D6 only where some gs_val are na. 
                if gs_val == 'nan':
                    # we ignore these alleles entirely. 'continue' to avoid tallying in any category
                    continue

                ans1 = compute_resolution(gs_val,pre_val1)
                ans2 = compute_resolution(gs_val,pre_val2)
                    
                zerodig,twodig,fourdig=single_allele_tally(ans1, ans2, zerodig, twodig, fourdig, gs_locus)

        # if we are working with d1-d4, the biallelic datasets
        else:
            
            # check accuracy, one locus (pair of alleles) at a time
            for i in range(1,len(genes),2):
                gs_val1 = gs_row[genes[i]].astype(str).values[0]
                gs_val2 = gs_row[genes[i+1]].astype(str).values[0]
                gs_locus = reformat( gs_val1.split("/")[0])[0] # stores first digit of gs_locus to differentiate class I vs II
                
                try:
                    # try to get the 2 pre alleles at the current iteration loci
                    pre_val1 = pre_row[genes[i]].astype(str).values[0]
                    pre_val2 = pre_row[genes[i+1]].astype(str).values[0]
                    
                    
                except: # exception occurs if pre alleles do not exist, in which case we will tally as nocall
                    # if allele in class I (ie, gs_locus is not D), increment nocall[1]
                    if gs_locus != "D": 
                        nocall[0] += 2
                        
                    # else (if allele is Class I), increment nocall[0]
                    else: 
                        nocall[1] += 2
                        
                    # skip rest of loop
                    continue
               
                    
                # handling cases where one or both alleles are "no call"
                no_val1 = False
                no_val2 = False

                if pre_val1 == None or pre_val1 == 'nan':
                    no_val1=True
                    if gs_locus != 'D':
                        nocall[0] += 1 
                    else:
                        nocall[1] += 1

                if pre_val2 == None or pre_val2 == 'nan':
                    no_val2=True
                    if gs_locus != 'D':
                        nocall[0] += 1 
                    else:
                        nocall[1] += 1

                #if both alleles are "no call" simply end this iteration
                if no_val1 and no_val2:
                    continue

                # if one allele is "no call", calculate accuracy as if a mono allele
                elif no_val1:
                    ans1 = compute_resolution(gs_val1,pre_val2)
                    ans2 = compute_resolution(gs_val1,pre_val2)
                    zerodig,twodig,fourdig=single_allele_tally(ans1,ans2,zerodig,twodig,fourdig, gs_locus)


                elif no_val2:
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val1,pre_val1)
                    zerodig,twodig,fourdig=single_allele_tally(ans1,ans2,zerodig,twodig,fourdig, gs_locus)

                else: # most typical case -- both calls are valid alleles

                    # assuming no swapping 
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val2,pre_val2)

                    # assuming swapping
                    ans3 = compute_resolution(gs_val1,pre_val2)
                    ans4 = compute_resolution(gs_val2,pre_val1)

                    zerodig,twodig,fourdig=pair_allele_tally(ans1,ans2,ans3,ans4,zerodig,twodig,fourdig, gs_locus)
 

    return zerodig,twodig,fourdig,nocall 

####################
#
# # BELOW FUNCTION COMMENTED AS A REMINDER TO REWRITE IT MIRRORING THE BUG FIXES TO COMPUTE_MATCHES
#
####################
#
# # requirements: gs accession numbers are under a column labeled "Run" 
# #pre accession numbers are under a column labeled "ERR" 
# # accession numbers/column titles are labeled identically between gold standard and results csv
# # Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Samples in GS, but not in PRE, are tallied in the "failed" variable 
# def get_inaccurate_and_all_alleles(pre,gs):

#     zerodig = []
#     all_alleles = [] # holds all alleles in gold standard
#     fail = 0

#     accession_numbers = gs['Run'].values.tolist()
#     genes = gs.columns.values.tolist()

#     for number in accession_numbers:
#         pre_row = pre.loc[pre['ERR'] == number]
#         gs_row = gs.loc[gs['Run'] == number]
        
        
#         # if we are working with d5 or d6, the monoallelic datasets
#         if (gs.columns.tolist() == ['Run', 'A', 'B'] or gs.columns.tolist() == ['Run', 'C']):
#             for i in range(1,len(genes)):
#                 gs_val = gs_row[genes[i]].astype(str).values[0]
#                 pre_val1 = pre_row[genes[i]].astype(str).values[0]
#                 pre_val2 = pre_row[genes[i]+".1"].astype(str).values[0]
                
#                 # if the gold standard contains many allele possibilities, and the caller is incorrect,
#                 # we will return only the first value in the gs
#                 gs_primary = reformat( gs_val.split("/")[0])
#                 all_alleles.append( gs_primary ) 

#                 ans1 = compute_resolution(gs_val,pre_val1)
#                 ans2 = compute_resolution(gs_val,pre_val2)
#                 if (max(ans1,ans2) == 0):
#                     zerodig.append(gs_primary)

#         # if we are working with d1-d4, the biallelic datasets
#         else:
#             for i in range(1,len(genes),2):
#                 try:
#                     gs_val1 = gs_row[genes[i]].astype(str).values[0]
#                     pre_val1 = pre_row[genes[i]].astype(str).values[0]
#                     gs_val2 = gs_row[genes[i+1]].astype(str).values[0]
#                     pre_val2 = pre_row[genes[i+1]].astype(str).values[0]

#                     if (gs_val1 == None) or (pre_val1 == None) or (gs_val2 == None) or (pre_val2 == None):
#                         fail = fail+1
#                         continue
                        
#                     # if the gold standard contains many allele possibilities, and the caller is incorrect,
#                     # we will return only the first value in the gs
#                     gs_primary1 = reformat( gs_val1.split("/")[0] ) 
#                     all_alleles.append(gs_primary1)
#                     gs_primary2 = reformat( gs_val2.split("/")[0] ) 
#                     all_alleles.append(gs_primary2)

#                     # assuming no swapping 
#                     ans1 = compute_resolution(gs_val1,pre_val1)
#                     ans2 = compute_resolution(gs_val2,pre_val2)

#                     # assuming swapping
#                     ans3 = compute_resolution(gs_val1,pre_val2)
#                     ans4 = compute_resolution(gs_val2,pre_val1)

#                     if (ans1+ans2 > ans3+ans4):
#                         if (ans1 == 0):
#                             zerodig.append(gs_primary1)
#                         if (ans2 == 0):
#                             zerodig.append(gs_primary2)
#                     else:
#                         if (ans3 == 0):
#                             zerodig.append(gs_primary1)
#                         if (ans4 == 0):
#                             zerodig.append(gs_primary2)
#                 except:
#                     fail = fail+1

#     return zerodig, all_alleles #,fail #onzero fail indicates exception occurred

# requirements: gs accession numbers are under a column labeled "Run" 
#pre accession numbers are under a column labeled "ERR" 
# accession numbers/column titles are labeled identically between gold standard and results csv
# Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Samples in GS, but not in PRE, are tallied in the "failed" variable 
def get_inaccurate_alleles(pre,gs):
    ret = get_inaccurate_and_all_alleles(pre,gs)
    return ret[0] 

# TODO: would this be better in the notebook ? is it really a global function
def sum_euro_groups(data,s=False):
    if s:
        ret = []
        for group in data:
            for allele in group:
                ret.append(allele)
        return ret
    else:
        ret = [0,0,0]
        for group in data:
            ret[0] += group[0]
            ret[1] += group[1]
            ret[2] += group[2]
        return ret