import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
import re
import io
import sqlite3
import requests
import networkx as nx

from Bio import SeqIO
from sklearn.metrics import roc_curve, auc
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from lifelines import KaplanMeierFitter

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

st.set_page_config(page_title="Genomic AI Platform", layout="wide")

st.title("🧬 Integrated Genomic Clinical Intelligence Platform")

model = joblib.load("huntington_ml_model.pkl")

# ---------------- DATABASE ---------------- #

conn = sqlite3.connect("patients.db")
cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS patients
(name TEXT, gene TEXT, repeat INT, result TEXT)
""")
conn.commit()

# ---------------- PATHOGENICITY MODEL ---------------- #

path_model = RandomForestClassifier()

X_train = pd.DataFrame({
"pos":[100,200,300,400,500],
"impact":[1,0,1,1,0]
})

y_train = [1,0,1,1,0]

path_model.fit(X_train,y_train)

# ---------------- DRUG MAP ---------------- #

drug_map = {
"HTT":"Tetrabenazine",
"BRCA1":"Olaparib",
"CFTR":"Ivacaftor"
}

# ---------------- TABS ---------------- #

tabs = st.tabs([
"Prediction",
"Batch",
"FASTA",
"VCF",
"Advanced Analytics",
"ROC"
])

# =========================================================
# PREDICTION
# =========================================================

with tabs[0]:

    gene = st.selectbox("Gene",["HTT","BRCA1","CFTR"])

    repeat = st.slider("Repeat Count",5,60,20)

    name = st.text_input("Patient Name")

    threshold = {"HTT":36,"BRCA1":10,"CFTR":15}[gene]

    if st.button("Predict Risk"):

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
            st.download_button("Download Report",f)

# =========================================================
# BATCH
# =========================================================

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
# =========================================================
# FASTA
# =========================================================

with tabs[2]:

    fasta = st.file_uploader("Upload FASTA", type=["fa","fasta"])

    if fasta:

        text_stream = io.StringIO(fasta.getvalue().decode())

        record = SeqIO.read(text_stream,"fasta")

        seq = str(record.seq)

        matches = list(re.finditer(r"(?:CAG){2,}", seq))

        pos = [(m.start(),(m.end()-m.start())//3) for m in matches]

        if pos:

            df = pd.DataFrame(pos, columns=["Position","Repeat"])

            st.dataframe(df)

            fig,ax = plt.subplots(figsize=(12,2))

            for p in df["Position"][:30]:
                ax.plot([p,p],[0,1],linewidth=2)

            ax.set_title("Genome Track")
            st.pyplot(fig)

# =========================================================
# VCF
# =========================================================

def annotate(chrom,pos):

    server="https://rest.ensembl.org"
    ext=f"/overlap/region/human/{chrom}:{pos}-{pos}"

    r=requests.get(server+ext,
                   headers={"Content-Type":"application/json"})

    if not r.ok:
        return "NA"

    data=r.json()

    return data[0]["external_name"] if data else "Intergenic"

with tabs[3]:

    vcf = st.file_uploader("Upload VCF", type=["vcf"])

    if vcf:

        lines = vcf.getvalue().decode().split("\n")

        samples=[]
        variants=[]

        for l in lines:

            if l.startswith("#CHROM"):
                samples=l.split("\t")[9:]

            if not l.startswith("#") and l.strip():

                c=l.split("\t")

                chrom=c[0]
                pos=int(c[1])

                variants.append([chrom,pos])

        df = pd.DataFrame(variants,columns=["Chrom","Pos"])

        df["Gene"]=df.apply(lambda x:
        annotate(x["Chrom"],x["Pos"]),axis=1)

        df["Drug"]=df["Gene"].map(drug_map)

        df["Pathogenicity"]=path_model.predict(
        pd.DataFrame({
        "pos":df["Pos"],
        "impact":[1]*len(df)
        }))

        prs = df["Gene"].map({"HTT":2.5,"BRCA1":3,"CFTR":1.8}).sum()

        st.write("Polygenic Risk Score:", prs)

        st.dataframe(df)

# =========================================================
# ADVANCED ANALYTICS
# =========================================================

with tabs[4]:

    st.subheader("CNV Simulation")

    cnv = pd.DataFrame({
    "Region":range(20),
    "CNV":np.random.choice(
    ["Deletion","Normal","Amplification"],20)
    })

    st.dataframe(cnv)

    st.subheader("Expression Heatmap")

    expr = np.random.randn(10,5)

    fig,ax=plt.subplots()

    ax.imshow(expr,cmap="coolwarm")

    st.pyplot(fig)

    st.subheader("Variant Network")

    G = nx.Graph()

    for i in range(10):
        G.add_edge(f"V{i}",f"Gene{i%3}")

    fig2 = plt.figure()

    nx.draw(G,with_labels=True)

    st.pyplot(fig2)

    st.subheader("PCA Clustering")

    X = np.random.randn(8,20)

    pca = PCA(n_components=2)

    comp = pca.fit_transform(X)

    fig3,ax3=plt.subplots()

    ax3.scatter(comp[:,0],comp[:,1])

    st.pyplot(fig3)

    st.subheader("Kaplan-Meier Survival")

    kmf = KaplanMeierFitter()

    t = np.random.exponential(10,50)
    e = np.random.choice([0,1],50)

    kmf.fit(t,e)

    fig4,ax4=plt.subplots()

    kmf.plot(ax=ax4)

    st.pyplot(fig4)

# =========================================================
# ROC
# =========================================================

with tabs[5]:

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

