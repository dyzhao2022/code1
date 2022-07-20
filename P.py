from pandas import Series, DataFrame

import scipy.stats as st

import operator

import os
from lifelines.datasets import load_waltons
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from lifelines.datasets import load_regression_dataset
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind




def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)

    fdr = p_vals * (len(p_vals) / ranked_p_values)
    fdr[fdr > 1] = 1
    # print(p_vals)
    # print( ranked_p_values)
    # print((len(p_vals) / ranked_p_values))
    # print(fdr)
    # print(pd.DataFrame([p_vals,fdr]))

    return fdr.tolist()


def qvalue(pvalues, method='BH'):
    """Calculate q-values by the Benjamini-Hochberg method.
    Args:
        pvalues (list): p-value list or array
        method (str): `BH` (Benjamini-Hochberg) or `BY` (Benjamini-Yekutieli)
    Returns:
        numpy.ndarray: q-value array
    """
    pvals = pvalues if isinstance(pvalues, np.ndarray) else np.array(pvalues)
    if np.sum((pvals < 0) | (pvals > 1)):
        raise ValueError('Invalid p-values')
    else:
        n = len(pvals)
        df_p = pd.DataFrame({'pval': pvals}).sort_values(
            by='pval', ascending=False
        ).reset_index()
    if method == 'BH':
        df_q = df_p.assign(
            qval=lambda d: (d['pval'] * n / (n - d.index.values)).cummin()
        )
    elif method == 'BY':
        w = np.sum(np.reciprocal(np.arange(1, n + 1, dtype='float32')))
        df_q = df_p.assign(
            qval=lambda d: (d['pval'] * n / (n - d.index.values) * w).cummin()
        )
    else:
        raise ValueError('Unimplemented method')
    return df_q.set_index(
        'index', drop=True
    ).sort_index()['qval'].clip(
        lower=0, upper=1
    ).values



def chuori(oritype,orifile):

    dfd=pd.DataFrame()

    if oritype=="Pro":


        df = pd.read_excel(dir + "/" + orifile, sheet_name=1, index_col=[0, 1, 2], header=0, usecols=range(0, 123))


        df.columns = shunxuspfall
        dfd = df.dropna(axis=0, how="all")


    elif oritype=="Rna":

        df = pd.read_excel(dir + "/" + orifile, sheet_name="gene_expression" ,index_col=[0,1,2],header=0,usecols=range(0,120))

        df.columns=[ij.replace("fpkm_","").replace("_","") for ij in df.columns]

        dfsrsm=df.loc[:,shunxusr]

        dfsrsm.columns=shunxusrfall



        dfd = dfsrsm.dropna(axis=0, how="all")


    elif oritype=="Phos":

        df = pd.read_excel(dir + "/" + orifile, sheet_name='Combine(>0.75)', index_col=[0, 1, 2, 3, 4], header=1,
                           usecols=[0, 1, 2, 3, 4] + list(range(245, 365)))

        df.columns = shunxuspfall
        dfd = df.dropna(axis=0, how="all")








    return dfd


def tcga(jian):

    wrr = []


    clioo = pd.read_csv(cfn, sep="\t", header=0, index_col=0)
    lins = clioo["Status"].values.tolist()
    tims = clioo["Time"].values.tolist()

    tim = []

    for ij, valuet in enumerate(tims):
        if "+" in valuet:
            tim.append(valuet.strip().replace("+", ""))
        else:
            tim.append(valuet.strip())

    clioo["Time2"] = tim



    tcga = pd.DataFrame(jian.values.T, index=jian.columns, columns=jian.index)
    # tcga.to_excel("tcga.xlsx")
    geall = pd.concat([tcga, clioo], axis=1, join="inner")  ####一个基因的所有信息

    cli = geall.iloc[:, [-1, -2]]
    uid = []
    for columnj in geall.iloc[:, :-10]:

        uid.append(columnj)
        d = geall[columnj].transpose()  #

        dfts = pd.concat([d, cli], axis=1, join="inner")  ####一个基因的所有信息

        dftps = dfts.dropna(axis=0, how="any")

        box1 = dftps[columnj].values.tolist()

        lq = np.quantile(box1, 0.25, interpolation="lower")

        hq = np.quantile(box1, 0.75, interpolation="higher")

        dfh = dftps[dftps[columnj] > hq]
        dfl = dftps[dftps[columnj] < lq]

        quantileguo = ["Low"] * len(dfh) + ["High"] * len(dfl)

        hl = pd.concat([dfl, dfh], axis=0)


        if hl.shape[0]<10:
            wrr.append(["NaN"] * 11)


        else:

            surt = dftps

            cph = CoxPHFitter()
            cph.fit(surt, duration_col='Time2', event_col='censorguo')

            coxhrsummary = cph.summary.iloc[0,
                           :].tolist()  ####coef,exp(coef),se(coef),coef lower 95%,coef upper 95%,exp(coef) lower 95%,exp(coef) upper 95%,z,p,-log2(p)


            hl["quantileguo"] = list(quantileguo)

            groups = hl['quantileguo']
            ix = (groups == 'High')

            results = logrank_test(hl['Time2'][ix].astype(int), hl['Time2'][~ix].astype(int),
                                   event_observed_A=hl['censorguo'][ix],
                                   event_observed_B=hl['censorguo'][~ix])





            wrr.append(coxhrsummary + [
                results.p_value])  ####coef,exp(coef),se(coef),coef lower 95%,coef upper 95%,exp(coef) lower 95%,exp(coef) upper 95%,z,p,-log2(p),logrank_testsurviveQ4p_value

    surtongji = pd.DataFrame(wrr)

    surtongji.columns = ["coef", "HRorexp(coef)", "se(coef)", "coef lower 95%", "coef upper 95%",
                         "exp(coef) lower 95%", "exp(coef) upper 95%", "z", "p", "-log2(p)",
                         "logrank_testsurviveQ4p_value"]  ##
    surtongji.index = uid



    return surtongji





dir = r"H:\zzu\TQTPST原始数据\msgyp\jjresult\0oridata"

#oritype, orifile,methfc,mratio,thrsam="Pro","Exp.All MS_identified_information.xlsx","mean",1.45,0####fc求平均为最后差异 >1.45 不限制样本

# oritype, orifile,methfc,mratio,thrsam="Rna","gene_expression.xlsx","median",2,5####fc求中位数为最后差异 >2 限制样本超过5对配对样本


# oritype, orifile,methfc,mratio,thrsam="Phos","Exp.PhosAll_MS_identified_information.xlsx","mean",1.45,0

output_file = orifile.replace(".xlsx", "") + "Q4xiuS167Ts.xlsx"




dfp = chuori(oritype, orifile)
hrdf = tcga(dfp)







dft= dfp.filter(regex='T$', axis=1)
dfn = dfp.filter(regex='P$', axis=1)
dft.columns=[j.replace("S","").replace("T","").replace("P","") for j in dft.columns.tolist()]
print(dft.columns)
print(dfn.columns)
dfn.columns = [j.replace("S", "").replace("T", "").replace("P", "") for j in dfn.columns.tolist()]
dftp= dft.div(dfn, axis=0)

if oritype=="Rna":


    dftp.replace([[0]],np.nan, inplace=True)
    dftp.replace([[np.inf, -np.inf]], np.nan, inplace=True)







samplesum=np.array(dftp.isnull().sum(axis=1))#按行统计缺失值个数

allprohang=np.repeat([dftp.shape[1]],dftp.shape[0])#统计行数

samplein=np.subtract(allprohang, samplesum)
#print(samplein)

tpmean  =dftp.mean(axis=1,skipna=True)
tpmedian = dftp.median(axis=1,skipna=True)


ptn = []
for ihang in range(dft.shape[0]):
    X1=dft.iloc[ihang].values.tolist()
    X2=dfn.iloc[ihang].values.tolist()
    print(len(X1), len(X2))
    if operator.eq(X1, X2):
        print(dft.iloc[ihang])
        ptn.append(1.0)

    else:

        ptn.append(st.mannwhitneyu(X1, X2)[1])


dftp['FCmeanTchuyiN'] = tpmean
dftp['FCmedianTchuyiN'] =tpmedian

dftp["Log2 FCmean"] = dftp['FCmeanTchuyiN'].apply(np.log2)
dftp["Log2 FCmedian"] = dftp['FCmedianTchuyiN'].apply(np.log2)


dftp["pvalue"] =ptn
dftp["fdrpvalue"] = fdr(ptn)##fdr(p_values),qvalue(p_values)
dftp["bhpvalue"] =  p_adjust_bh(ptn)
dftp["ziyouduo"] = samplein


dffinal=pd.concat([dfp, dftp, hrdf], axis=1, join="inner")
dffinal.to_excel(output_file)



dffinalpv = dffinal[(dffinal["bhpvalue"] < 0.05) & (dffinal["ziyouduo"]>thrsam)]





if oritype=="Rna":
    dffinalpvu = dffinalpv[dffinalpv["FCmedianTchuyiN"] > mratio  ]
    dffinalpvd = dffinalpv[dffinalpv["FCmedianTchuyiN"] < 1/mratio]


else:
    dffinalpvu = dffinalpv[dffinalpv["FCmeanTchuyiN"] > mratio  ]
    dffinalpvd = dffinalpv[dffinalpv["FCmeanTchuyiN"] < 1/mratio]

dffinalpvu.to_excel(
         output_file.replace(".xlsx","")+oritype+"UP.xlsx")

dffinalpvd.to_excel(
    output_file.replace(".xlsx", "") + oritype + "DOWN.xlsx")

