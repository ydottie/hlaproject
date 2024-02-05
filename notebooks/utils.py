import pandas as pd

# formats alleles into barebones 4-char strings. for ex, input s=A*01:03 returns 0103
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

# compares gs_val to pre_val and returns 0 for miscall, 1 for one-field accurate, and 2 for two-field accurate
def compute_resolution(gs_val,pre_val):
    gs_all = gs_val.split("/")
    pre_fixed = fix(pre_val)
    flag = False #false = inaccurate, turn to true = One field accurate
    
    for val in gs_all:
        gs_fixed = fix(val)

        if (gs_fixed[0:2] == pre_fixed[0:2]):
            if (gs_fixed[2:4] == pre_fixed[2:4]):
                return 2
            flag = True # use flag, rather than "return 2", because there can be multiple val to compare
        
    return 1 if flag==True else 0

# reformats val to exact X*XX:XX format to standardize and prevent unexpected bugs/symbols
def reformat(gs_val):
    numbers = fix(gs_val)
    if (gs_val[0] == 'D'): 
        return gs_val[0:4]+'*'+numbers[0:2]+':'+numbers[2:4]
    else:
        return gs_val[0]+'*'+numbers[0:2]+':'+numbers[2:4]
    

# updates accuracy for **PAIRED** allele comparison (2 pre - 2 gs)
def pair_allele_tally(ans1,ans2,ans3,ans4,zerofield,onefield,twofield,gs_locus):
    
    # parallel comparison
    if (ans1+ans2 > ans3+ans4):
        
        # assign accuracy for allele 1
        if (ans1 == 0):
            if gs_locus != 'D':
                zerofield[0]  = zerofield[0] + 1 
            else:
                zerofield[1] = zerofield[1] + 1
        if (ans1 == 1):
            if gs_locus != 'D':
                onefield[0] = onefield[0] + 1 
            else:
                onefield[1] = onefield[1] + 1
        if (ans1 == 2):
            if gs_locus != 'D':
                twofield[0] = twofield[0] + 1 
            else:
                twofield[1] = twofield[1] + 1
        
        # assign accuracy for allele 2
        if (ans2 == 0):
            if gs_locus != 'D':
                zerofield[0] = zerofield[0] + 1 
            else:
                zerofield[1] = zerofield[1] + 1
        if (ans2 == 1):
            if gs_locus != 'D':
                onefield[0] = onefield[0] + 1 
            else:
                onefield[1] = onefield[1] + 1
        if (ans2 == 2):
            if gs_locus != 'D':
                twofield[0] = twofield[0] + 1 
            else:
                twofield[1] = twofield[1] + 1
                
    # crosswise comparison
    else:
        # assign accuracy for allele 1
        if (ans3 == 0):
            if gs_locus != 'D':
                zerofield[0] = zerofield[0] + 1 
            else:
                zerofield[1] = zerofield[1] + 1
        if (ans3 == 1):
            if gs_locus != 'D':
                onefield[0] = onefield[0] + 1 
            else:
                onefield[1] = onefield[1] + 1
        if (ans3 == 2):
            if gs_locus != 'D':
                twofield[0] = twofield[0] + 1 
            else:
                twofield[1] = twofield[1] + 1
                
        # assign accuracy for allele 2
        if (ans4 == 0):
            if gs_locus != 'D':
                zerofield[0] = zerofield[0] + 1 
            else:
                zerofield[1] = zerofield[1] + 1
        if (ans4 == 1):
            if gs_locus != 'D':
                onefield[0] = onefield[0] + 1 
            else:
                onefield[1] = onefield[1] + 1
        if (ans4 == 2):
            if gs_locus != 'D':
                twofield[0] = twofield[0] + 1 
            else:
                twofield[1] = twofield[1] + 1
    
    return zerofield,onefield,twofield

# updates accuracy for **SINGLE** allele comparison (works with either 2 pre - 1 gs, or 1 gs - 2 pre)
def single_allele_tally(ans1,ans2,zerofield,onefield,twofield,gs_locus):
    ans = max(ans1,ans2)
    
    if (ans == 0):
        if gs_locus != 'D':
            zerofield[0] = zerofield[0] + 1 
        else:
            zerofield[1] = zerofield[1] + 1
    if (ans == 1):
        if gs_locus != 'D':
            onefield[0] = onefield[0] + 1 
        else:
            onefield[1] = onefield[1] + 1
    if (ans == 2):
        if gs_locus != 'D':
            twofield[0] = twofield[0] + 1 
        else:
            twofield[1] = twofield[1] + 1
    return zerofield,onefield,twofield

###############################   
# MASTER ACCURACY FUNCTION 
###############################   
# For all purpose accuracy calculation. 
# REQUIREMENTS: 
# 1. gs accession numbers are under a column labeled "Run" 
# 2. pre accession numbers are under a column labeled "ERR" 
# 3. accession numbers/column titles are labeled identically between gs and pre csv's
# INPUTS:
# Pandas dataframes for gs and pre
# OUTPUTS:
# zerofield, onefield, twofield, and nocall: integer lists of length 2, index 0 holds class I counts, index 1 holds class II counts
# all_gs_alleles, all_pre_alleles, miscalled_gs_alleles: string list of alleles
# NOTES:
# cells in PRE, but not in GS are ignored
# cells in GS, but not in PRE, are tallied in the "nocall" variable 
###############################
def master_accuracy_function(pre,gs):
    
    # fixing empty cells in case of homozygous allele (bug fix suggested by mourisl in github issue #1)
    for gene in ["A", "B", "C", "DQB1", "DRB1"]:
        target = gene + ".1"
        if (target not in pre.columns):
            continue
        pre[target].fillna(pre[gene], inplace=True)

    # Initializing variables
    twofield = [0,0]
    onefield = [0,0]
    zerofield = [0,0]
    nocall = [0,0]

    d5 = ['Run', 'A', 'B']
    d6 = ['Run', 'C']
    
    all_gs_alleles = [] # to hold all alleles in gold standard
    all_pre_alleles = [] # to hold all alleles in pre
    miscalled_gs_alleles = [] # to hold only miscalled alleles in gold standard

    accession_numbers = gs["Run"].values.tolist()
    genes = gs.columns.values.tolist()

    # iterate through accesssion numbers to check accuracy one sample at a time
    for number in accession_numbers:
        
        # get the gs and pre alleles for the sample
        pre_row = pre.loc[pre['ERR'] == number]
        gs_row = gs.loc[gs['Run'] == number]

        # if we are working with d5 or d6, the monoallelic datasets
        if (gs.columns.tolist() == d5 or gs.columns.tolist() == d6):
            
            for i in range(1,len(genes)):
                
                # get values for gs, pre, and gs_locus
                gs_val = gs_row[genes[i]].astype(str).values[0]
                
                # addresses edge case in D6 where some gs_val are na. 
                if gs_val == 'nan':
                    # we ignore these alleles entirely. 'continue' to avoid tallying in any category
                    continue
                    
                gs_locus = gs_val[0] # stores first digit of gs_locus to differentiate class I vs II
                all_gs_alleles.append( gs_val ) 
                    
                pre_val1 = pre_row[genes[i]].astype(str).values[0]
                pre_val2 = pre_row[genes[i]+".1"].astype(str).values[0]
                
                all_pre_alleles.append( pre_val1 )
                all_pre_alleles.append( pre_val2 )

                
                # compute the accuracy resolution
                ans1 = compute_resolution(gs_val,pre_val1)
                ans2 = compute_resolution(gs_val,pre_val2)
                
                # update the results variables
                zerofield,onefield,twofield=single_allele_tally(ans1, ans2, zerofield, onefield, twofield, gs_locus)
                if (max(ans1,ans2) == 0):
                     miscalled_gs_alleles.append( gs_val )

        # if we are working with d1-d4, the biallelic datasets
        else:
            
            # check accuracy, one locus (pair of alleles) at a time
            for i in range(1,len(genes),2):
                
                # get values for gs, pre, and gs_locus
                gs_val1 = gs_row[genes[i]].astype(str).values[0]
                gs_val2 = gs_row[genes[i+1]].astype(str).values[0]
                gs_locus = reformat( gs_val1.split("/")[0])[0] # stores first digit of gs_locus to differentiate class I vs II
                all_gs_alleles.append( reformat( gs_val1.split("/")[0]) ) 
                all_gs_alleles.append( reformat( gs_val2.split("/")[0]) ) 
                
                 #addresses edge case in only read length gs where one or more loci is nan
                if gs_val1 == 'nan' and gs_val2 == 'nan':
                    # we ignore these loci entirely. 'continue' to avoid tallying in any category
                    continue
                
                if gs_val1 == 'nan' or gs_val2 == 'nan':
                    print ("entered gs_val1 xor gs_val2 does not exist case. TODO")
                    # TODO: not sure if this happens yet, if it does, I will implement it
                
                try:
                    # try to get the 2 pre alleles at the current iteration loci
                    pre_val1 = pre_row[genes[i]].astype(str).values[0]                    
                    pre_val2 = pre_row[genes[i+1]].astype(str).values[0]
                    
                    all_pre_alleles.append( pre_val1 )
                    all_pre_alleles.append( pre_val2 )
                    
                    
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
                    ans2 = compute_resolution(gs_val2,pre_val2)
                    
                    if (max(ans1,ans2) == 0):
                        if ans1 == 0:
                            miscalled_gs_alleles.append(reformat( gs_val1.split("/")[0]))
                        else: 
                            miscalled_gs_alleles.append(reformat( gs_val2.split("/")[0]))
                    zerofield,onefield,twofield=single_allele_tally(ans1,ans2,zerofield,onefield,twofield, gs_locus)

                elif no_val2:
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val2,pre_val1)
                    
                    if (max(ans1,ans2) == 0):
                        if ans1 == 0:
                            miscalled_gs_alleles.append(reformat( gs_val1.split("/")[0]))
                        else: 
                            miscalled_gs_alleles.append(reformat( gs_val2.split("/")[0]))
                    zerofield,onefield,twofield=single_allele_tally(ans1,ans2,zerofield,onefield,twofield, gs_locus)

                else: # most typical case -- both calls are valid alleles

                    # assuming no swapping 
                    ans1 = compute_resolution(gs_val1,pre_val1)
                    ans2 = compute_resolution(gs_val2,pre_val2)

                    # assuming swapping
                    ans3 = compute_resolution(gs_val1,pre_val2)
                    ans4 = compute_resolution(gs_val2,pre_val1)

                    
                    if (ans1+ans2 > ans3+ans4):
                        if (ans1 == 0):
                            miscalled_gs_alleles.append(reformat( gs_val1.split("/")[0]))
                        if (ans2 == 0):
                            miscalled_gs_alleles.append(reformat( gs_val2.split("/")[0]))
                    else:
                        if (ans3 == 0):
                            miscalled_gs_alleles.append(reformat( gs_val2.split("/")[0]))
                        if (ans4 == 0):
                            miscalled_gs_alleles.append(reformat( gs_val1.split("/")[0]))

                    zerofield,onefield,twofield=pair_allele_tally(ans1,ans2,ans3,ans4,zerofield,onefield,twofield, gs_locus)
 

    return zerofield,onefield,twofield,nocall,miscalled_gs_alleles,all_gs_alleles,all_pre_alleles


def get_miscalled_alleles_only(pre,gs):
    ret = master_accuracy_function(pre,gs)
    return ret[4] 


def get_accuracy_counts(pre,gs):
    ret = master_accuracy_function(pre,gs)
    return ret[:4] 

def get_miscalled_and_all_alleles(pre,gs):
    '''returns miscalled_gs_alleles,all_gs_alleles'''
    ret = master_accuracy_function(pre,gs)
    return ret[4:6]

# split predictions into europe and yoruba df, for ancestry analysis
def split_csv_by_ancestry():
    ''' returns europe_df, yoruba_df '''
    # can modify filepaths as necessary 
    groupscsv = "../datasets/SraRunTableD1.txt"
    goldstandard = "../datasets/1_gs.csv"
    
    gs = pd.read_csv(goldstandard)
    groups = pd.read_csv(groupscsv)

    dfs = []

    for group, df_by_group in groups.groupby('Population'):
        accession_numbers = df_by_group['Run'].values.tolist()
        gs_final = gs[gs['Run'].isin(accession_numbers)] #gs_final is a df containing the gold standard samples per population group
        dfs.append(gs_final)

    europe_df = pd.concat([dfs[0],dfs[1],dfs[2],dfs[3]])
    yoruba_df =dfs[4]
    
    return europe_df, yoruba_df