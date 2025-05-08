import pandas as pd
import numpy as np
import random
from collections import Counter
import matplotlib.pyplot as plt

def get_Fragments_DF(DF,label):
    #Create new DatFrame
    Fragments=DF[['length',label]]
    #Fragments['length']=Fragments['length'].astype('int')
    Fragments.set_index('length',inplace=True)

    #Delete zeros
    Fragments = Fragments[(Fragments != 0).all(axis=1)]

    #Sort Fragments
    Fragments=Fragments.sort_index()
    return Fragments

def get_median(DF,label):
    DF=DF.sort_index()
    yy=DF[label]/np.sum(DF[label])
    cum_prob=np.cumsum(yy)
    index=np.argmax(cum_prob>=0.5)
    return DF.index[index]

def get_median0(DF,label):
    DF=DF.sort_values(by='length')
    yy=DF[label]/np.sum(DF[label])
    cum_prob=np.cumsum(yy)
    index=np.argmax(cum_prob>=0.5)
    return DF['length'][index]

def load_data(file_path_nontreated, file_path_treated,label,Total):
    data_labels=['0','CD4T', 'hESC-H9', 'HEK293T', 'RNH2-KO-T3-8','RNH2-KO-T3-17']
    new_labels=['length']+data_labels[1::]

    # Read the CSV file and create DataFrames
    NonTreated = pd.read_csv(file_path_nontreated,usecols=data_labels).dropna()
    Treated = pd.read_csv(file_path_treated,usecols=data_labels).dropna()

    NonTreated.columns=new_labels
    Treated.columns=new_labels
    

    NonTreated['length']=NonTreated['length']
    NonTreated[label]=np.round(NonTreated[label]*Total/NonTreated[label].sum())
    NonTreated['length']=np.array(np.round(NonTreated['length']),dtype='int')
    NonTreated=NonTreated.groupby('length').sum().reset_index()

    Treated['length']=Treated['length']
    Treated[label]=np.round(Treated[label]*Total/Treated[label].sum())
    Treated['length']=np.array(np.round(Treated['length']),dtype='int')
    Treated=Treated.groupby('length').sum().reset_index()
    
    return Treated, NonTreated


def cut_given_length0(length, min_length,n_cuts):
    minimum=0
    num_tries=0
    while minimum<=min_length and num_tries<=10:
        num_tries+=1
        #Find random positions for the cuts
        cut_pos=random.sample(range(min_length,length-min_length+1),k=n_cuts)
        cut_pos.sort()

        #Computes the resulting lengths
        lengths=np.diff([0]+cut_pos+[length])
        minimum=min(lengths)
    return list(lengths)

def cut_given_length(length, min_length,n_cuts):
    extra_space=length-min_length*(n_cuts+1)
    cuts_extra_space=n_cuts
    cuts_pos=random.choices(range(extra_space+1),k=cuts_extra_space)
    cuts_pos.sort()
    lengths_extra=list(np.diff([0]+cuts_pos+[extra_space]))
    lengths_extra_modified=np.array(lengths_extra)
    lengths=np.array(lengths_extra_modified)+min_length
    return list(lengths)

def run_simulation(Fragments_input,label, N,min_size,total_cuts=0):
    #Copy Data Frame
    Fragments=Fragments_input.copy()
    
    #Only consider Fragments big enough
    lengths=np.array(Fragments[Fragments.index>=2*min_size].index)
    counts=np.array(Fragments[Fragments.index>=2*min_size][label])
    
    #Compute weights
    weights=lengths*counts

    # Randomly choose N segments on the given weights
    random_lengths = random.choices(lengths, weights=weights, k=N)

    #Get Counter object
    counter_lengths=Counter(random_lengths)
    
    #Make the cuts
    new_lengths=[]
    for length in counter_lengths:
        cuts_count=counter_lengths[length]
        max_cuts_per_fragment=int(np.floor(length/min_size))-1
        n_fragments=int(Fragments.at[length,label])
        
        max_cuts_possible=max_cuts_per_fragment*n_fragments
        
        #Randomly choose without replacement, which fragments to cut out of (max_cuts_per_fragment)
        #many copies of each fragment.
        random_indices=random.sample(range(max_cuts_possible),min(cuts_count,max_cuts_possible))
        
        #Transform into indices of the fragments
        random_indices=[x//max_cuts_per_fragment for x in random_indices]
        
        
        #Convert into a counter object
        counter_indices=Counter(random_indices)
        

        for n_cuts in counter_indices.values():
            Fragments.at[length,label]-=1
            lengths_cut=cut_given_length(length,min_size,n_cuts)
            new_lengths+=lengths_cut
            total_cuts+=len(lengths_cut)-1
            

    #Update the lengths in Fragments
    counter_new_lengths=Counter(new_lengths)
    for length in counter_new_lengths:
        if length in Fragments.index:
            Fragments.at[length,label]+=counter_new_lengths[length]
        else:
            Fragments.at[length,label]=counter_new_lengths[length]
    Fragments=Fragments[Fragments[label]!=0]
    return Fragments,total_cuts,new_lengths

def binary_search_number_cuts(label, Treated, NonTreated, Fragments,
                              Total,min_size,n_reps,threshold_N,N_known=None):

    Goal_Median=get_median0(Treated,label)
    Initial_Median=get_median0(NonTreated,label)

    print('Initial Median ', Initial_Median)
    print('Goal Median ', Goal_Median)
    print('------------')

    min_N=0; max_N=int(np.floor(2*Total/3))
    if N_known is None:
        N=int(Total/5)
    else:
        N=N_known
        min_N=N
        max_N=N
    continue_search=True

    while(continue_search):
        medians_try=[]
        Fragments_try=[]
        total_cuts_try=[]
        new_lengths_try=[]
        print('N=', N, end=' ')
        for rep_index in range(n_reps):
            print('*', end='')
            New_Fragments=Fragments.copy()
            New_Fragments,total_cuts,new_lengths=run_simulation(New_Fragments, label,N,min_size,total_cuts=0)
            Fragments_try.append(New_Fragments.copy())
            total_cuts_try.append(total_cuts)
            new_lengths_try.append(new_lengths.copy())
            medians_try.append(get_median(New_Fragments,label))  
        
        average_median=np.mean(medians_try)
        print('. Mean:', int(average_median))
        if average_median<Goal_Median-1:
            max_N=N
            new_N=int((max_N+min_N)/2)
        elif average_median>Goal_Median+1:
            min_N=N
            new_N=int((max_N+min_N)/2)
        else:
            break
        if abs(new_N-N)/N>=threshold_N:
            continue_search=True
            N=new_N
        else:
            continue_search=False
    return N,Fragments_try, total_cuts_try, new_lengths_try,medians_try, Goal_Median, Initial_Median

def print_plot_results(label, Treated,NonTreated, Fragments, Fragments_try, 
                       medians_try, total_cuts_try, Total,Goal_Median, Initial_Median, n_reps):
    print('Final Distribution')
    print('Label: ',label)
    print('Median before ', get_median(Fragments,label))
    print('Goal Median ', Goal_Median)
    print('New medians: ', medians_try)
    print('New medians average: ', np.mean(medians_try))
    print('Number of cuts/Fragments', np.average(total_cuts_try)/Total)

    Fragments_Treated=get_Fragments_DF(Treated,label)

    #Original Distro NonTreated
    df=Fragments.copy()
    df=df.reset_index()
    df['x_group'] = ((df["length"]) // 1000)*1000 # Group x into bins of 10
    df = df.groupby('x_group')[label].sum().reset_index()

    nontreated_x=np.array(df['x_group'])
    nontreated_y=np.array(df[label])/df[label].sum()
    plt.plot(nontreated_x,nontreated_y)


    #Original Distro Treated
    df=Fragments_Treated.copy()
    df=df.reset_index()
    df['x_group'] = ((df["length"]) // 1000)*1000 # Group x into bins of 10
    df = df.groupby('x_group')[label].sum().reset_index()

    treated_x=np.array(df['x_group'])
    treated_y=np.array(df[label])/df[label].sum()
    plt.plot(treated_x,treated_y)

    #New Distribution
    for New_Fragments in Fragments_try:
        New_Fragments = New_Fragments[(New_Fragments != 0).all(axis=1)]
        df=New_Fragments.copy()
        df=df.reset_index()
        df['x_group'] = ((df["length"]) // 1000)*1000 # Group x into bins of 10
        df = df.groupby('x_group')[label].sum().reset_index()

        simulated_x=np.array(df['x_group'])
        simulated_y=np.array(df[label])/df[label].sum()

        plt.plot(simulated_x,simulated_y)
    plt.legend(['Non Treated', 'Treated']+['Simulation'+str(i+1) for i in range(n_reps)])
    plt.show()