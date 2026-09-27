#!/usr/bin/env python3
"""E6 per ADDENDUM_5 + round-03 weakness-12: post-hoc recalibration
(DISCLOSED post-hoc) of the locked H1' OOF predictions. Isotonic and Platt
per ordinal threshold (P(y>=k), k=1..3), 5-fold CV inside the OOF set so
the recalibrated probabilities stay out-of-sample; report calibration
slope + Brier score pre/post per threshold. Paper keeps NOT-decision-grade
language regardless of improvement."""
import json
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
df=pd.read_csv('results/h1prime_oof_predictions.csv')
y=df.y.values; out={'spec':'ADDENDUM_5 E6; recalibration disclosed post-hoc; CV-internal','thresholds':{}}
def slope_brier(yy,pp):
    eps=1e-6; lp=np.log(np.clip(pp,eps,1-eps)/(1-np.clip(pp,eps,1-eps)))
    s=LogisticRegression().fit(lp.reshape(-1,1),yy)
    return float(s.coef_[0][0]), float(np.mean((yy-pp)**2))
for k,col in [(1,'oof_p_ge1'),(2,'oof_p_ge2'),(3,'oof_p_ge3')]:
    yy=(y>=k).astype(int); pp=df[col].values
    s0,b0=slope_brier(yy,pp)
    res={}
    for name,fitter in [('isotonic',lambda a,b:IsotonicRegression(out_of_bounds='clip').fit(a,b)),
                        ('platt',lambda a,b:LogisticRegression().fit(a.reshape(-1,1),b))]:
        pr=np.zeros_like(pp)
        for tr,te in KFold(5,shuffle=True,random_state=20260926).split(pp):
            m=fitter(pp[tr],yy[tr])
            pr[te]=m.predict(pp[te].reshape(-1,1)) if name=='platt' else m.predict(pp[te])
        s1,b1=slope_brier(yy,pr)
        res[name]={'slope_post':s1,'brier_post':b1}
    out['thresholds'][f'p_ge{k}']={'slope_pre':s0,'brier_pre':b0,**res,
        'base_rate':float(yy.mean())}
json.dump(out,open('results/e6_recalibration.json','w'),indent=2)
print(json.dumps(out,indent=1))
