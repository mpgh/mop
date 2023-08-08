import numpy as np



def assess_anomaly(KS_test,SW_test,chi2dof):


    status = 'Non Anomalous'
    
    if (KS_test<0.025) & (SW_test<0.025):
    
        status = 'Anomalous'
        
    if (chi2dof>5):
    
        status = 'Anomalous'
        
    return status
    
    
