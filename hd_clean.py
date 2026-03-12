import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
import re
import io
import sqlite3
import requests
from Bio import SeqIO
from sklearn.metrics import roc_curve, auc
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

st.set_page_config(page_title="Genomic AI Platform", layout="wide")

st.title("🧬 Genomic Variant Intelligence Platform")

model = joblib.load("huntington_ml_model.pkl")

# DATABASE
conn = sqlite3.connect("patients.db")
cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS patients
(name TEXT, gene TEXT, repeat INT, result TEXT)
""")
conn.commit()

tab1, tab2, tab3, tab4, tab5 = st.tabs([
"Single Prediction",
"Batch Prediction",
"FASTA Analysis",
"VCF Variant Analysis",
"ROC Dashboard"
])

# ---------------- SINGLE ---------------- #

with tab1:

    gene = st.selectbox("Gene",["HTT","BRCA1","CFTR"])

    repeat = st.slider("Repeat Count",5,60,20)

    name = st.text_input("Patient Name")

    threshold = {"HTT":36,"BRCA1":10,"CFTR":15}[gene]

    if st.button("Predict"):

        result = "Risk" if repeat>=threshold else "Normal"

        st.write(result)

        cursor.execute(
        "INSERT INTO patients VALUES (?,?,?,?)",
        (name,gene,repeat,result)
        )

        conn.commit()

        fig,ax=plt.subplots()
        ax.bar(["Repeat"],[repeat])
        ax.axhline(threshold,color="red")
        st.pyplot(fig)

        pdf="report.pdf"
        c=canvas.Canvas(pdf,pagesize=letter)
        c.drawString(100,700,f"Patient:{name}")
        c.drawString(100,650,f"Gene:{gene}")
        c.drawString(100,600,f"Repeat:{repeat}")
        c.drawString(100,550,f"Result:{result}")
        c.save()

        with open(pdf,"rb") as f:
            st.download_button("Download Report",f,"report.pdf")

# ---------------- BATCH ---------------- #

with tab2:

    csv=st.file_uploader("Upload CSV",type=["csv"])

    if csv:

        df=pd.read_csv(csv)

        df["Prediction"]=(df["CAG_Repeats"]>=36).astype(int)

        st.dataframe(df)

# ---------------- FASTA ---------------- #

with tab3:

    fasta=st.file_uploader("Upload FASTA",type=["fasta","fa"])

    if fasta:

        text_stream=io.StringIO(fasta.getvalue().decode())

        record=SeqIO.read(text_stream,"fasta")

        seq=str(record.seq)

        matches=list(re.finditer(r"(?:CAG){2,}",seq))

        pos=[(m.start(),(m.end()-m.start())//3) for m in matches]

        if pos:

            df=pd.DataFrame(pos,columns=["Position","Repeat"])

            st.dataframe(df)

            fig,ax=plt.subplots()
            ax.scatter(df["Position"],df["Repeat"])
            st.pyplot(fig)

# ---------------- VCF ---------------- #

def annotate(chrom,pos):

    server="https://rest.ensembl.org"
    ext=f"/overlap/region/human/{chrom}:{pos}-{pos}"

    r=requests.get(server+ext,
                   headers={"Content-Type":"application/json"})

    if not r.ok:
        return "NA"

    data=r.json()

    return data[0]["external_name"] if data else "Intergenic"

def drug(gene):

    m={"HTT":"Tetrabenazine",
       "BRCA1":"Olaparib",
       "CFTR":"Ivacaftor"}

    return m.get(gene,"Consult Doctor")

with tab4:

    vcf=st.file_uploader("Upload VCF",type=["vcf"])

    if vcf:

        lines=vcf.getvalue().decode().split("\n")

        variants=[]

        for l in lines:

            if not l.startswith("#") and l.strip():

                c=l.split("\t")

                variants.append([c[0],c[1]])

        df=pd.DataFrame(variants,
                        columns=["Chrom","Pos"])

        df["Gene"]=df.apply(
        lambda x: annotate(x["Chrom"],x["Pos"]),axis=1)

        df["Drug"]=df["Gene"].apply(drug)

        st.dataframe(df)

# ---------------- ROC ---------------- #

with tab5:

    x=np.linspace(15,55,200)

    probs=model.predict_proba(
    pd.DataFrame({"CAG_Repeats":x}))[:,1]

    y=(x>=36).astype(int)

    fpr,tpr,_=roc_curve(y,probs)

    auc_score=auc(fpr,tpr)

    fig,ax=plt.subplots()

    ax.plot(fpr,tpr,label=f"AUC={auc_score:.2f}")
    ax.plot([0,1],[0,1],'--')

    ax.legend()

    st.pyplot(fig)
