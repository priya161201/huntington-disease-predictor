import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import joblib
import io
import re
import requests
import networkx as nx

from Bio import SeqIO
from sklearn.metrics import roc_curve, auc
from sklearn.decomposition import PCA
from lifelines import KaplanMeierFitter
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

st.set_page_config(layout="wide")

st.title("🧬 CRAZY Genomic Clinical Intelligence Platform")

model = joblib.load("huntington_ml_model.pkl")

# ================= DATABASE =================

conn = sqlite3.connect("patients.db")
cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS patients
(name TEXT, gene TEXT, repeat INT, result TEXT)
""")
conn.commit()

# ================= FUNCTIONS =================

def make_pdf(name,gene,repeat,result):

    file="patient_report.pdf"

    c=canvas.Canvas(file,pagesize=letter)

    c.drawString(100,750,"Genomic Risk Report")

    c.drawString(100,700,f"Patient: {name}")
    c.drawString(100,650,f"Gene: {gene}")
    c.drawString(100,600,f"Repeat: {repeat}")
    c.drawString(100,550,f"Result: {result}")

    c.save()

    return file


def history_pdf(df):

    file="history.pdf"

    c=canvas.Canvas(file,pagesize=letter)

    y=750

    for i,row in df.iterrows():

        c.drawString(50,y,
        f"{row['name']} | {row['gene']} | {row['repeat']} | {row['result']}")

        y-=20

        if y<100:
            c.showPage()
            y=750

    c.save()

    return file


def annotate(chrom,pos):

    try:

        r=requests.get(
        f"https://rest.ensembl.org/overlap/region/human/{chrom}:{pos}-{pos}",
        headers={"Content-Type":"application/json"})

        data=r.json()

        return data[0]["external_name"] if data else "Intergenic"

    except:
        return "NA"

drug_map={
"HTT":"Tetrabenazine",
"BRCA1":"Olaparib",
"CFTR":"Ivacaftor"
}

# ================= TABS =================

tabs=st.tabs([
"Prediction",
"History",
"Batch",
"FASTA",
"VCF",
"Analytics",
"ROC"
])

# =====================================================
# PREDICTION
# =====================================================

with tabs[0]:

    gene=st.selectbox("Gene",["HTT","BRCA1","CFTR"])

    repeat=st.slider("Repeat",5,60,20)

    name=st.text_input("Patient")

    thr={"HTT":36,"BRCA1":10,"CFTR":15}[gene]

    if st.button("Predict"):

        result="Risk" if repeat>=thr else "Normal"

        st.write(result)

        cursor.execute(
        "INSERT INTO patients VALUES (?,?,?,?)",
        (name,gene,repeat,result)
        )
        conn.commit()

        pdf=make_pdf(name,gene,repeat,result)

        with open(pdf,"rb") as f:
            st.download_button("Download Report",f)

# =====================================================
# HISTORY
# =====================================================

with tabs[1]:

    df=pd.read_sql("SELECT * FROM patients",conn)

    st.dataframe(df)

    if st.button("Download History PDF"):

        file=history_pdf(df)

        with open(file,"rb") as f:
            st.download_button("Download",f)

# =====================================================
# BATCH
# =====================================================

with tabs[2]:

    csv=st.file_uploader("CSV",type=["csv"])

    if csv:

        df=pd.read_csv(csv)

        df["Prediction"]=(df["CAG_Repeats"]>=36).astype(int)

        st.dataframe(df)

# =====================================================
# FASTA
# =====================================================

with tabs[3]:

    fasta=st.file_uploader("FASTA",type=["fa","fasta"])

    if fasta:

        txt=io.StringIO(fasta.getvalue().decode())

        record=SeqIO.read(txt,"fasta")

        seq=str(record.seq)

        m=list(re.finditer(r"(?:CAG){2,}",seq))

        pos=[(x.start(),(x.end()-x.start())//3) for x in m]

        if pos:

            df=pd.DataFrame(pos,columns=["Position","Repeat"])

            st.dataframe(df)

            fig,ax=plt.subplots(figsize=(12,2))

            for p in df["Position"][:40]:
                ax.plot([p,p],[0,1])

            st.pyplot(fig)

# =====================================================
# VCF
# =====================================================

with tabs[4]:

    vcf=st.file_uploader("VCF",type=["vcf"])

    if vcf:

        lines=vcf.getvalue().decode().split("\n")

        vars=[]

        for l in lines:

            if not l.startswith("#") and l.strip():

                c=l.split("\t")

                vars.append([c[0],int(c[1])])

        df=pd.DataFrame(vars,columns=["Chrom","Pos"])

        df["Gene"]=df.apply(lambda x:
        annotate(x["Chrom"],x["Pos"]),axis=1)

        df["Drug"]=df["Gene"].map(drug_map)

        prs=df["Gene"].map({"HTT":2.5,"BRCA1":3,"CFTR":1.8}).sum()

        st.write("PRS:",prs)

        st.dataframe(df)

# =====================================================
# ANALYTICS
# =====================================================

with tabs[5]:

    st.subheader("CNV")

    cnv=pd.DataFrame({
    "Region":range(20),
    "Type":np.random.choice(["Del","Amp","Norm"],20)
    })

    st.dataframe(cnv)

    st.subheader("Expression Heatmap")

    expr=np.random.randn(10,5)

    fig,ax=plt.subplots()

    ax.imshow(expr)

    st.pyplot(fig)

    st.subheader("Network")

    G=nx.Graph()

    for i in range(10):
        G.add_edge(f"V{i}",f"G{i%3}")

    fig2=plt.figure()

    nx.draw(G,with_labels=True)

    st.pyplot(fig2)

    st.subheader("PCA")

    X=np.random.randn(8,20)

    comp=PCA(2).fit_transform(X)

    fig3,ax3=plt.subplots()

    ax3.scatter(comp[:,0],comp[:,1])

    st.pyplot(fig3)

    st.subheader("Survival")

    km=KaplanMeierFitter()

    t=np.random.exponential(10,50)
    e=np.random.choice([0,1],50)

    km.fit(t,e)

    fig4,ax4=plt.subplots()

    km.plot(ax=ax4)

    st.pyplot(fig4)

# =====================================================
# ROC
# =====================================================

with tabs[6]:

    x=np.linspace(15,55,200)

    p=model.predict_proba(
    pd.DataFrame({"CAG_Repeats":x}))[:,1]

    y=(x>=36).astype(int)

    fpr,tpr,_=roc_curve(y,p)

    fig,ax=plt.subplots()

    ax.plot(fpr,tpr)
    ax.plot([0,1],[0,1],'--')

    st.pyplot(fig)
