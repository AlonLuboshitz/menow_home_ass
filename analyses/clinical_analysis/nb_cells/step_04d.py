ps = dex[~dex['pred'].isin(['No predisposition','Unknown'])]
ps = ps.groupby('pred').agg(Count=('PATIENT_ID','nunique'),MedAge=('AGE','median'),
    PctF=('SEX',lambda x: (x=='Female').sum()/len(x)*100 if len(x)>0 else 0)).reset_index().sort_values('Count',ascending=False)
ps['%Patients'] = (ps['Count']/total_p*100).round(1)
ps.columns = ['Predisposition','Count','Median Age','% Female','% Patients']
print(ps.to_string(index=False))
