"""
Created on Mon Mar 23 15:51:39 2020

@author: Rabi Chhantyal-Pun
email: rcpchem@gmail.com, rc13564@bristol.ac.uk
Please use for only non-commerical use 
"""

from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tkinter import filedialog
from lmfit import Model, Parameters, minimize, fit_report
import matplotlib.backends.backend_tkagg
from datetime import datetime

root=Tk()

Label(root, text="Reaction Model").grid(row=0,column=0,padx=5,pady=5)
rxn_textbox=Text(root, width=30, height=10, padx=5)
rxn_textbox.insert('1.0', 'CH2OO + CH2OO = Product1 : k1\nCH2OO = Product2  : k2')
rxn_textbox.grid(row=0,column=0,columnspan=5,padx=5,pady=5)
#rxn=rxn_textbox.get("1.0","end-1c")

Label(root, text="Reaction Species").grid(row=5,column=0,padx=5,pady=5)
species_entry=Entry(root)
species_entry.insert(0,'CH2OO:Product1:Product2')
species_entry.grid(row=5,column=1,padx=5)

Label(root, text="Initial Concentrations, C0").grid(row=6,column=0,padx=5,pady=5)
initial_C0_params_entry=Entry(root)
initial_C0_params_entry.insert(0,'1E12:0:0')
initial_C0_params_entry.grid(row=6,column=1,padx=5)

Label(root, text="Initial Rate Coefficients, k").grid(row=7,column=0,padx=5,pady=5)
initial_k_params_entry=Entry(root)
initial_k_params_entry.insert(0,'1E-10:50')
initial_k_params_entry.grid(row=7,column=1,padx=5)

Label(root, text="Time Start").grid(row=8,column=0,padx=5,pady=5)
time_start_entry=Entry(root)
time_start_entry.insert(0, '0')
time_start_entry.grid(row=8,column=1,padx=5)

Label(root, text="Time Stop").grid(row=9,column=0,padx=5,pady=5)
step_size_entry=Entry(root)
step_size_entry.insert(0, '0.01')
step_size_entry.grid(row=9,column=1,padx=5)
#
Label(root, text="Number of Steps").grid(row=10,column=0,padx=5,pady=5)
step_num_entry=Entry(root)
step_num_entry.insert(0, '100')
step_num_entry.grid(row=10,column=1,padx=5)

Label(root, text="Measured Species").grid(row=5,column=3,padx=5,pady=5)
measured_species_entry=Entry(root)
measured_species_entry.insert(0,'None')
measured_species_entry.grid(row=5,column=4,padx=5)

Label(root, text="C0 min -Fit").grid(row=6,column=3,padx=5,pady=5)
min_C0_params_entry=Entry(root)
min_C0_params_entry.insert(0,'1E12:0:0')
min_C0_params_entry.grid(row=6,column=4,padx=5)

Label(root, text="C0 max -Fit").grid(row=7,column=3,padx=5,pady=5)
max_C0_params_entry=Entry(root)
max_C0_params_entry.insert(0,'1E13:0:0')
max_C0_params_entry.grid(row=7,column=4,padx=5)

Label(root, text="k min -Fit").grid(row=8,column=3,padx=5,pady=5)
min_k_params_entry=Entry(root)
min_k_params_entry.insert(0,'1E-11:10')
min_k_params_entry.grid(row=8,column=4,padx=5)

Label(root, text="k max -Fit").grid(row=9,column=3,padx=5,pady=5)
max_k_params_entry=Entry(root)
max_k_params_entry.insert(0,'1E-9:1000')
max_k_params_entry.grid(row=9,column=4,padx=5)

Label(root, text="Species to Plot").grid(row=10,column=3,padx=5,pady=5)
plot_species_entry=Entry(root)
plot_species_entry.insert(0,'CH2OO:Product1:Product2')
plot_species_entry.grid(row=10,column=4,padx=5)

    
def model_exp_callback():

    initial_C0_params=initial_C0_params_entry.get()
    initial_C0_params=[float(x) for x in initial_C0_params.split(':')]

    t=np.linspace(float(time_start_entry.get()),float(step_size_entry.get()),int(step_num_entry.get()))
    
    rxn_list=rxn_textbox.get("1.0","end-1c")
    rxn_list_split=rxn_list.split('\n')
    RRX=[]*len(rxn_list_split)
    for i in range(len(rxn_list_split)):
        R_i=rxn_list_split[i].split()
        RRX.append(R_i)
    
    species_list=species_entry.get()
    S=species_list.split(':')
    
    
    def rxn(Conc,t):        
        #Dictionary for assigning species strings to concentrations 
        conc_dict={}
        for i in range(len(S)):
            conc_dict[S[i]]=Conc[i]
        for key,val in conc_dict.items():
            exec(key + '=val')
        
        #Dictionary for assigning parameter strings to parameter values         
        initial_k_params=initial_k_params_entry.get()
        initial_k_params=initial_k_params.split(':')
        initial_k_params=[float(x) for x in initial_k_params]
        params=[]*len(RRX)
        for j in range(len(RRX)):
            params.append(RRX[j][len(RRX[j])-1])
        params_dict={}
        for k in range(len(params)):
            params_dict[params[k]]=initial_k_params[k]
        for key,val in params_dict.items():
            exec(key + '=val')


#        rate expression for reactants for all reactions
#        unimolecular and bimolecular reactions    
#        R_A=[]*len(R) 
#        for i in range(len(R)):
#            if R[i].index('=') == 3:
#                AA=R[i][len(R[i])-1]+str('*')+R[i][0]+str('*')+R[i][2]
#            elif R[i].index('=') == 1:
#                AA=R[i][len(R[i])-1]+str('*')+R[i][0]
#            R_A.append(AA)
        #Reactions of all molecularity
        R_A=[]*len(RRX)
        for i in range(len(RRX)):
            R_A_i=[]*(len(RRX[i][0:RRX[i].index('=')])-RRX[i][0:RRX[i].index('=')].count('+')+1) #list for reaction species and k
            R_A_i.append(RRX[i][len(RRX[i])-1])
            for j in range(len(RRX[i][0:RRX[i].index('=')])-RRX[i][0:RRX[i].index('=')].count('+')):
                R_A_i.append('*'+RRX[i][j+j])
            R_A.append(''.join(R_A_i))
 
        #Combination of rate expression to give ode for all species           
        Sdt=[]*len(S)     
        for j in range(len(S)):
            Sdt_R=[]*len(RRX)
            for k in range(len(RRX)):
                if S[j] not in RRX[k]:
                    continue
                for l in range(len(RRX[k])):
                    if RRX[k][l]==S[j] and l < RRX[k].index('='):
                        Sdt_R.append('-'+R_A[k])
                    elif RRX[k][l]==S[j] and l > RRX[k].index('='):
                        Sdt_R.append('+'+R_A[k])
            Sdt.append(Sdt_R)                               
        Sdt = [''.join(x) for x in Sdt]
        Sdt_eval=[]*len(Sdt)
        for x in range(len(Sdt)):
            Sdt_eval.append(eval(Sdt[x]))
        return Sdt_eval        
        
    Conc=odeint(rxn,initial_C0_params,t)
    
    plot_S=plot_species_entry.get()
    plot_S=plot_S.split(':')
    plot_S_num=[S.index(x) for x in plot_S]

    #Plotting modelled and exp data       
    exp=measured_species_entry.get()
    exp_names=exp.split(':')
    if exp == 'None':
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(len(plot_S)):
            ax.plot(t,Conc[:,plot_S_num[i]],label=plot_S[i])
            plt.ylabel('[Concentration] (Unit)',fontsize=15)
            plt.xlabel('Time (s)',fontsize=15)
            plt.legend()
    else:
        root.filename=filedialog.askopenfilename(initialdir='/Documents/Python', title='Select a file')
        data=np.genfromtxt(root.filename)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(len(plot_S)):
            ax.plot(t,Conc[:,plot_S_num[i]],label=plot_S[i])
            plt.ylabel('[Concentration] (Unit)',fontsize=15)
            plt.xlabel('Time (s)',fontsize=15)
            plt.legend()
        for j in range(len(data[0,:])-1):
            ax.plot(data[:,0],data[:,j+1],label='Exp '+exp_names[j])
            plt.ylabel('[Concentration] (Unit)',fontsize=15)
            plt.xlabel('Time (s)',fontsize=15)
            plt.legend()
    plt.show()
    print('Success!')
    
    
    #saving model results and log   
    response = messagebox.askquestion('Save?','Do you want to save file?')
    if response == 'yes':
        fname=simpledialog.askstring('Filename', 'Please enter filename')
        np.savetxt(fname+'.txt', np.column_stack((t,Conc)), delimiter="\t", fmt='%s')
        
        log_file='Successful run at '+str(datetime.now())+'\n'
        log_file=log_file+'\n'+'Reaction Model\n'+rxn_list+'\n'
        log_file=log_file+'Reactions Species\n'+species_entry.get()+'\n'
        log_file=log_file+'Initial Concentrations\n'+initial_C0_params_entry.get()+'\n'
        log_file=log_file+'Rate Coefficient Values\n'+initial_k_params_entry.get()
        f = open(fname+'_log.txt', 'w')
        f.write(log_file) 
    
def fit_callback():
    
    exp=measured_species_entry.get()
    S_m=exp.split(':')
    
    if exp == 'None':
        print('Give Measured Species labels')
        return
    del exp
    
    #Experimental data input
    root.filename=filedialog.askopenfilename(initialdir='/Documents/Python', title='Select a file')
    data=np.genfromtxt(root.filename)
    
    #Experimental time measurement
    t=data[:,0]


    #dictionatory for measured signals
    sig_dic={}    
    for i in range(len(S_m)):
        sig_dic[S_m[i]]=data[:,i+1]

    
    #Reaction Model and Species from user input    
    rxn_list=rxn_textbox.get("1.0","end-1c")
    rxn_list_split=rxn_list.split('\n')
    RRX=[]*len(rxn_list_split)
    for i in range(len(rxn_list_split)):
        R_i=rxn_list_split[i].split()
        RRX.append(R_i)
    
    #List of species from user input
    species_list=species_entry.get()
    S=species_list.split(':')

    #Packing Parameters for lmfit module
    initial_k_params=initial_k_params_entry.get()
    initial_k_params=initial_k_params.split(':')
    initial_k_params=[float(x) for x in initial_k_params] 
    min_k_params=min_k_params_entry.get()
    min_k_params=min_k_params.split(':')
    min_k_params=[float(x) for x in min_k_params]
    max_k_params=max_k_params_entry.get()
    max_k_params=max_k_params.split(':')
    max_k_params=[float(x) for x in max_k_params]
    
    initial_C0_params=initial_C0_params_entry.get()
    initial_C0_params=initial_C0_params.split(':')
    initial_C0_params=[float(x) for x in initial_C0_params]
    min_C0_params=min_C0_params_entry.get()
    min_C0_params=min_C0_params.split(':')
    min_C0_params=[float(x) for x in min_C0_params]
    max_C0_params=max_C0_params_entry.get()
    max_C0_params=max_C0_params.split(':')
    max_C0_params=[float(x) for x in max_C0_params]
    
    
    #vary or fix parameters?
    kk=[]*len(max_k_params)
    for tt in range(len(min_k_params)):
        if min_k_params[tt]==max_k_params[tt]:
            kk.append('false')
        else:
            kk.append('true')
    del tt
    
    cc=[]*len(max_C0_params)
    for tt in range(len(min_C0_params)):
        if min_C0_params[tt]==max_C0_params[tt]:
            cc.append('false')
        else:
            cc.append('true')
    del tt
   
    #Collection of rate and concentration parameter from reaction model input for fitting
    params = Parameters()
    params_names=[]*len(RRX)
    for l in range(len(RRX)):
        params_names.append(RRX[l][len(RRX[l])-1])
    del l

    for p in range(len(params_names)):
        if kk[p]=='true':
            params.add(str(params_names[p]), value=initial_k_params[p], min=min_k_params[p], max=max_k_params[p])
        elif kk[p]=='false':
            params.add(str(params_names[p]), value=initial_k_params[p], vary=False)
    del p
    
    species_names=species_entry.get()
    species_names=species_names.split(':')
    species_names=['Initial_' + s for s in species_names]
    
    for p in range(len(species_names)):
        if cc[p]=='true':
            params.add(species_names[p], value=initial_C0_params[p], min=min_C0_params[p], max=max_C0_params[p])
        elif cc[p]=='false':
            params.add(species_names[p], value=initial_C0_params[p], vary=False)
    del p
    
    def rxn_fit(params,t,sig_dic):
        
        
        def rxn(Conc,t):  
            
            #Relating params values to rate coefficient labels in the model   
            params_dict={}
            for k in range(len(params_names)):
                params_dict[params_names[k]]=params[params_names[k]]
            for key,val in params_dict.items():
                exec(key + '=val')
            
            #Relating Conc values to S labels in the model        
            conc_dict={}
            for i in range(len(S)):
                conc_dict[S[i]]=Conc[i]
            for key,val in conc_dict.items():
                exec(key + '=val')
    
            #Rate expression for reactants for all reactions
            R_A=[]*len(RRX)
            for i in range(len(RRX)):
                R_A_i=[]*(len(RRX[i][0:RRX[i].index('=')])-RRX[i][0:RRX[i].index('=')].count('+')+1) #list for reaction species and k
                R_A_i.append(RRX[i][len(RRX[i])-1])
                for j in range(len(RRX[i][0:RRX[i].index('=')])-RRX[i][0:RRX[i].index('=')].count('+')):
                    R_A_i.append('*'+RRX[i][j+j])
                R_A.append(''.join(R_A_i))
     
            #Combination of rate expression to give ode for all species           
            Sdt=[]*len(S)     
            for j in range(len(S)):
                Sdt_R=[]*len(RRX)
                for k in range(len(RRX)):
                    if S[j] not in RRX[k]:
                        continue
                    for l in range(len(RRX[k])):
                        if RRX[k][l]==S[j] and l < RRX[k].index('='):
                            Sdt_R.append('-'+R_A[k])
                        elif RRX[k][l]==S[j] and l > RRX[k].index('='):
                            Sdt_R.append('+'+R_A[k])
                Sdt.append(Sdt_R)                               
            Sdt = [''.join(x) for x in Sdt]
            Sdt_eval=[]*len(Sdt)
            for x in range(len(Sdt)):
                Sdt_eval.append(eval(Sdt[x]))
            return Sdt_eval
        
        #list of initial concentration values for all the reaction species 
        C0=[]*len(species_names)
        for p in range(len(species_names)):
            conc=params[species_names[p]]
            C0.append(conc)
            
        Conc=odeint(rxn,C0,t)
    
        residual=[]*len(S_m)
        for i in range(len(S_m)):
            residual.append(Conc[:,S.index(S_m[i])]-sig_dic[S_m[i]])
        return residual
    

    out = minimize(rxn_fit, params, args=(t, sig_dic), method='leastsq')
    print('Success!')
    print(fit_report(out.params)) 

    #Modelled concentration with optimized parameters
    def rxn_final(Conc,t):  
        
        #Relating Conc values to S labels in the model        
        conc_dict={}
        for i in range(len(S)):
            conc_dict[S[i]]=Conc[i]
        for key,val in conc_dict.items():
            exec(key + '=val')
            
     
        #Relating params values to rate coefficient labels in the model    
        out_params_dict={}
        for k in range(len(params_names)):
            out_params_dict[params_names[k]]=out.params[params_names[k]].value
        for key,val in out_params_dict.items():
            exec(key + '=val')
    
        #Rate expression for reactants for all reactions
        R_A=[]*len(RRX)
        for i in range(len(RRX)):
            if RRX[i].index('=') == 3:
                AA=RRX[i][len(RRX[i])-1]+str('*')+RRX[i][0]+str('*')+RRX[i][2]
            elif RRX[i].index('=') == 1:
                AA=RRX[i][len(RRX[i])-1]+str('*')+RRX[i][0]
            R_A.append(AA)
     
        #Combination of rate expression to give ode for all species           
        Sdt=[]*len(S)     
        for j in range(len(S)):
            Sdt_R=[]*len(RRX)
            for k in range(len(RRX)):
                if S[j] not in RRX[k]:
                    continue
                for l in range(len(RRX[k])):
                    if RRX[k][l]==S[j] and l < RRX[k].index('='):
                        Sdt_R.append('-'+R_A[k])
                    elif RRX[k][l]==S[j] and l > RRX[k].index('='):
                        Sdt_R.append('+'+R_A[k])
            Sdt.append(Sdt_R)                               
        Sdt = [''.join(x) for x in Sdt]
        Sdt_eval=[]*len(Sdt)
        for x in range(len(Sdt)):
            Sdt_eval.append(eval(Sdt[x]))
        return Sdt_eval   
     
    #list of initial concentration values for all the reaction species 
    C0=[]*len(species_names)
    for p in range(len(species_names)):
        conc=out.params[species_names[p]].value 
        C0.append(conc) 
    


    t1=np.linspace(np.min(t),np.max(t),np.int(step_num_entry.get())) 
    
    Conc=odeint(rxn_final,C0,t1)
    

    #Plotting modelled and exp data       
    plot_S=plot_species_entry.get()
    plot_S=plot_S.split(':')
    plot_S_num=[S.index(x) for x in plot_S]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(plot_S)):
        ax.plot(t1,Conc[:,plot_S_num[i]],label='Fit Model '+plot_S[i])
        plt.ylabel('[Concentration] (Unit)',fontsize=15)
        plt.xlabel('Time (s)',fontsize=15)
        plt.legend()
    del i    
    for i in range(len(S_m)):
        ax.plot(t,sig_dic[S_m[i]],label='Exp '+S_m[i])
        plt.ylabel('[Concentration] (Unit)',fontsize=15)
        plt.xlabel('Time (s)',fontsize=15)
        plt.legend()
    del i 
    plt.show()
    
    
    #saving model results and log   
    response = messagebox.askquestion('Save?','Do you want to save file?')
    if response == 'yes':
        fname=simpledialog.askstring('Filename', 'Please enter filename')
        np.savetxt(fname+'.txt', np.column_stack((t1,Conc)), delimiter="\t", fmt='%s')
        
        log_file='Successful run at '+str(datetime.now())+'\n'
        log_file=log_file+'\n'+'Reaction Model\n'+rxn_list+'\n'
        log_file=log_file+'Reactions Species\n'+species_entry.get()+'\n'
        log_file=log_file+'Initial Concentrations\n'+initial_C0_params_entry.get()+'\n'
        log_file=log_file+'Inital Rate Coefficient Values\n'+initial_k_params_entry.get()+'\n'
        log_file=log_file+'Min Initial Concentration Values\n'+min_C0_params_entry.get()+'\n'
        log_file=log_file+'Max Initial Concentration Values\n'+max_C0_params_entry.get()+'\n'
        log_file=log_file+'Min Rate Coefficient Values\n'+min_k_params_entry.get()+'\n'
        log_file=log_file+'Max Rate Coefficient Values\n'+max_k_params_entry.get()+'\n'
        log_file=log_file+'\n'+'Fit Result'+'\n'+str(fit_report(out.params))
        f = open(fname+'_log.txt', 'w')
        f.write(log_file) 
    
file_button=Button(root,text='Fit',command=fit_callback)
file_button.config(height=3,width=15)
file_button.grid(row=11,column=4)

model_exp_button=Button(root,text='Model',command=model_exp_callback)
model_exp_button.config(height=3,width=15)
model_exp_button.grid(row=11,column=3)

quit_button=Button(root, text="Quit", command=root.destroy)
quit_button.config(height=3,width=15)
quit_button.grid(row=11,column=0)

root.mainloop()


    