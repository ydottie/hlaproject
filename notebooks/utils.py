def fix(s):
    firstcolon = s.find(":")
    if firstcolon == -1:
        star = s.find("*")
        s_new = s[star+1:star+3]
    else:
        s_new = s[firstcolon-2:firstcolon] + s[firstcolon+1:firstcolon+3]
        
    s_new = s_new.replace("*", "0")
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
            flag = True
        
    return 2 if flag==True else 0

# formats val to the proper X*XX:XX format
def reformat(gs_val):
    numbers = fix(gs_val)
    if (gs_val[0] == 'D'):
        return gs_val[0:4]+'*'+numbers[0:2]+':'+numbers[2:4]
    else:
        return gs_val[0]+'*'+numbers[0:2]+':'+numbers[2:4]

# requirements: gs accession numbers are under a column labeled "Run" 
#pre accession numbers are under a column labeled "ERR" 
# accession numbers/column titles are labeled identically between gold standard and results csv
# Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Samples in GS, but not in PRE, are tallied in the "failed" variable 
def compute_matches(pre,gs):

    # index 0 is class I, index 1 is class II
    fourdig = [0,0]
    twodig = [0,0]
    zerodig = [0,0]
    fail = 0

    d5 = ['Run', 'A', 'B']
    d6 = ['Run', 'C']
    classI = ['A','B','C']

    accession_numbers = gs["Run"].values.tolist()
    genes = gs.columns.values.tolist()

    for number in accession_numbers:
        pre_row = pre.loc[pre['ERR'] == number]
        gs_row = gs.loc[gs['Run'] == number]

        # if we are working with d5 or d6, the monoallelic datasets
        if (gs.columns.tolist() == d5 or gs.columns.tolist() == d6):
            for i in range(1,len(genes)):
                gs_val = gs_row[genes[i]].astype(str).values[0]
                pre_val1 = pre_row[genes[i]].astype(str).values[0]
                pre_val2 = pre_row[genes[i]+".1"].astype(str).values[0]

                ans = max(compute_resolution(gs_val,pre_val1), compute_resolution(gs_val,pre_val2))
                if (ans == 0):
                    zerodig[0] = zerodig[0] + 1
                if (ans == 2):
                    twodig[0] = twodig[0] + 1
                if (ans == 4):
                    fourdig[0] = fourdig[0] + 1

        # if we are working with d1-d4, the biallelic datasets
        else:
            for i in range(1,len(genes),2):
                try:
                    gs_val1 = gs_row[genes[i]].astype(str).values[0]
                    pre_val1 = pre_row[genes[i]].astype(str).values[0]
                    gs_val2 = gs_row[genes[i+1]].astype(str).values[0]
                    pre_val2 = pre_row[genes[i+1]].astype(str).values[0]
                    
                    if gs_val1 == None or pre_val1 == None or gs_val2 == None or pre_val2 == None:
                        fail = fail+1
                        continue
                    
                    if gs_val1 == 'nan' or pre_val1 == 'nan':
                        fail = fail+1
                        continue
                    
                    gs_primary = reformat( gs_val1.split("/")[0])

                    # assuming no swapping 
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val2,pre_val2)

                    # assuming swapping
                    ans3 = compute_resolution(gs_val1,pre_val2)
                    ans4 = compute_resolution(gs_val2,pre_val1)

                    if (ans1+ans2 > ans3+ans4):
                        if (ans1 == 0):
                            if classI.contains(gs_primary[0]):
                                zerodig[0] += 1 
                            else:
                                zerodig[1] += 1
                        if (ans1 == 2):
                            if classI.contains(gs_primary[0]):
                                twodig[0] = twodig[0] + 1 
                            else:
                                twodig[1] = twodig[1] + 1
                        if (ans1 == 4):
                            if gs_primary[0] in classI:
                                fourdig[0] = fourdig[0] + 1 
                            else:
                                fourdig[1] = fourdig[1] + 1
                        if (ans2 == 0):
                            if gs_primary[0] in classI:
                                zerodig[0] = zerodig[0] + 1 
                            else:
                                zerodig[1] = zerodig[1] + 1
                        if (ans2 == 2):
                            if gs_primary[0] in classI:
                                twodig[0] = twodig[0] + 1 
                            else:
                                twodig[1] = twodig[1] + 1
                        if (ans2 == 4):
                            if gs_primary[0] in classI:
                                fourdig[0] = fourdig[0] + 1 
                            else:
                                fourdig[1] = fourdig[1] + 1
                    else:
                        if (ans3 == 0):
                            if gs_primary[0] in classI:
                                zerodig[0] = zerodig[0] + 1 
                            else:
                                zerodig[1] = zerodig[1] + 1
                        if (ans3 == 2):
                            if gs_primary[0] in classI:
                                twodig[0] = twodig[0] + 1 
                            else:
                                twodig[1] = twodig[1] + 1
                        if (ans3 == 4):
                            if gs_primary[0] in classI:
                                fourdig[0] = fourdig[0] + 1 
                            else:
                                fourdig[1] = fourdig[1] + 1
                        if (ans4 == 0):
                            if gs_primary[0] in classI:
                                zerodig[0] = zerodig[0] + 1 
                            else:
                                zerodig[1] = zerodig[1] + 1
                        if (ans4 == 2):
                            if gs_primary[0] in classI:
                                twodig[0] = twodig[0] + 1 
                            else:
                                twodig[1] = twodig[1] + 1
                        if (ans4 == 4):
                            if gs_primary[0] in classI:
                                fourdig[0] = fourdig[0] + 1 
                            else:
                                fourdig[1] = fourdig[1] + 1
                except:
                    fail = fail+1

    return zerodig,twodig,fourdig #,fail #onzero fail indicates exception occurred

# requirements: gs accession numbers are under a column labeled "Run" 
#pre accession numbers are under a column labeled "ERR" 
# accession numbers/column titles are labeled identically between gold standard and results csv
# Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Samples in GS, but not in PRE, are tallied in the "failed" variable 
def get_inaccurate_and_all_alleles(pre,gs):

    zerodig = []
    all_alleles = [] # holds all alleles in gold standard
    fail = 0

    accession_numbers = gs['Run'].values.tolist()
    genes = gs.columns.values.tolist()

    for number in accession_numbers:
        pre_row = pre.loc[pre['ERR'] == number]
        gs_row = gs.loc[gs['Run'] == number]
        
        
        # if we are working with d5 or d6, the monoallelic datasets
        if (gs.columns.tolist() == ['Run', 'A', 'B'] or gs.columns.tolist() == ['Run', 'C']):
            for i in range(1,len(genes)):
                gs_val = gs_row[genes[i]].astype(str).values[0]
                pre_val1 = pre_row[genes[i]].astype(str).values[0]
                pre_val2 = pre_row[genes[i]+".1"].astype(str).values[0]
                
                # if the gold standard contains many allele possibilities, and the caller is incorrect,
                # we will return only the first value in the gs
                gs_primary = reformat( gs_val.split("/")[0])
                all_alleles.append( gs_primary ) 

                ans1 = compute_resolution(gs_val,pre_val1)
                ans2 = compute_resolution(gs_val,pre_val2)
                if (max(ans1,ans2) == 0):
                    zerodig.append(gs_primary)

        # if we are working with d1-d4, the biallelic datasets
        else:
            for i in range(1,len(genes),2):
                try:
                    gs_val1 = gs_row[genes[i]].astype(str).values[0]
                    pre_val1 = pre_row[genes[i]].astype(str).values[0]
                    gs_val2 = gs_row[genes[i+1]].astype(str).values[0]
                    pre_val2 = pre_row[genes[i+1]].astype(str).values[0]

                    if (gs_val1 == None) or (pre_val1 == None) or (gs_val2 == None) or (pre_val2 == None):
                        fail = fail+1
                        continue
                        
                    # if the gold standard contains many allele possibilities, and the caller is incorrect,
                    # we will return only the first value in the gs
                    gs_primary1 = reformat( gs_val1.split("/")[0] ) 
                    all_alleles.append(gs_primary1)
                    gs_primary2 = reformat( gs_val2.split("/")[0] ) 
                    all_alleles.append(gs_primary2)

                    # assuming no swapping 
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val2,pre_val2)

                    # assuming swapping
                    ans3 = compute_resolution(gs_val1,pre_val2)
                    ans4 = compute_resolution(gs_val2,pre_val1)

                    if (ans1+ans2 > ans3+ans4):
                        if (ans1 == 0):
                            zerodig.append(gs_primary1)
                        if (ans2 == 0):
                            zerodig.append(gs_primary2)
                    else:
                        if (ans3 == 0):
                            zerodig.append(gs_primary1)
                        if (ans4 == 0):
                            zerodig.append(gs_primary2)
                except:
                    fail = fail+1

    return zerodig, all_alleles #,fail #onzero fail indicates exception occurred

# requirements: gs accession numbers are under a column labeled "Run" 
#pre accession numbers are under a column labeled "ERR" 
# accession numbers/column titles are labeled identically between gold standard and results csv
# Only accuracy for samples in both GS and PRE are calculated. Samples in PRE, but not in GS are ignored. Samples in GS, but not in PRE, are tallied in the "failed" variable 
def get_inaccurate_alleles(pre,gs):
    ret = get_inaccurate_and_all_alleles(pre,gs)
    return ret[0] 

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