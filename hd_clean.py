import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, roc_curve, auc
import plotly.graph_objects as go
from fpdf import FPDF

st.set_page_config(page_title="Huntington Platform", layout="wide")

# ================= DATABASE =================
@st.cache_resource
def get_connection():
    return sqlite3.connect("patients.db", check_same_thread=False)

conn = get_connection()
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS patients(
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT,
    age INTEGER,
    cag INTEGER,
    risk TEXT
)
""")
conn.commit()

# ================= LOGIN =================
if "login" not in st.session_state:
    st.session_state.login = False

with st.sidebar:
    st.title("Login")
    user = st.text_input("Username")
    pw = st.text_input("Password", type="password")

    if st.button("Login"):
        if user == "doctor" and pw == "123":
            st.session_state.login = True
            st.rerun()
        else:
            st.error("Invalid Credentials")

if not st.session_state.login:
    st.stop()

# ================= RISK FUNCTION =================
def risk_predict(cag):
    if cag < 20:
        return "No Risk"
    elif cag <= 37:
        return "Risk"
    else:
        return "High Risk"

# ================= UI =================
st.title("🧬 Huntington Disease Risk Prediction Platform")

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Prediction", "Patient Database", "Genome Visualization", "Model Analysis", "Medical Report"
])

# ================= TAB 1: PREDICTION =================
with tab1:
    gene_file = st.file_uploader("Upload Gene Expression File", type=["csv"])
    clinical_file = st.file_uploader("Upload Clinical Data File", type=["csv"])

    if gene_file and clinical_file:
        try:
            gene_df = pd.read_csv(gene_file)
            clinical_df = pd.read_csv(clinical_file)

            # Simulated feature extraction
            gene_df["CAG"] = np.random.randint(15, 55, len(gene_df))

            results = gene_df.copy()
            results["Risk"] = results["CAG"].apply(risk_predict)

            # 🔥 FIXED SORTING ISSUE
            results["Sample"] = results.index.astype(str)
            results["Sample_num"] = results["Sample"].str.extract(r'(\d+)')
            results["Sample_num"] = results["Sample_num"].fillna(0).astype(int)

            results = results.sort_values(by="Sample_num")\
                             .drop(columns=["Sample_num"])\
                             .reset_index(drop=True)

            st.session_state["results"] = results

            st.success("Prediction Completed")
            st.dataframe(results)

            # Save to DB (clear old → avoid duplicates)
            c.execute("DELETE FROM patients")
            for _, row in results.iterrows():
                c.execute(
                    "INSERT INTO patients(name, age, cag, risk) VALUES (?, ?, ?, ?)",
                    ("Sample", 30, int(row["CAG"]), row["Risk"])
                )
            conn.commit()

        except Exception as e:
            st.error(str(e))

# ================= TAB 2: DATABASE =================
with tab2:
    st.subheader("Patient Records")
    try:
        df = pd.read_sql("SELECT * FROM patients ORDER BY id DESC", conn)
        st.dataframe(df)
    except:
        st.warning("No data available")

# ================= TAB 3: GENOME VIS =================
with tab3:
    if "results" in st.session_state and not st.session_state["results"].empty:
        df = st.session_state["results"]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df.index,
            y=df["CAG"],
            mode='markers',
            marker=dict(color=df["CAG"], colorscale='Viridis')
        ))

        fig.update_layout(title="Genome CAG Distribution")
        st.plotly_chart(fig)
    else:
        st.info("Run prediction first")

# ================= TAB 4: MODEL ANALYSIS =================
with tab4:
    if "results" in st.session_state and not st.session_state["results"].empty:
        df = st.session_state["results"]

        X = df[["CAG"]]
        y = df["Risk"].map({"No Risk": 0, "Risk": 1, "High Risk": 1})

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

        model = LogisticRegression()
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        cm = confusion_matrix(y_test, y_pred)

        st.write("Confusion Matrix")
        st.write(cm)

        y_prob = model.predict_proba(X_test)[:,1]
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        roc_auc = auc(fpr, tpr)

        st.write(f"ROC AUC: {roc_auc:.2f}")

    else:
        st.info("Run prediction first")

# ================= TAB 5: PDF REPORT =================
with tab5:
    if st.button("Generate PDF Report"):
        try:
            df = pd.read_sql("SELECT * FROM patients ORDER BY id DESC LIMIT 1", conn)

            if len(df) == 0:
                st.warning("No data")
            else:
                row = df.iloc[0]

                pdf = FPDF()
                pdf.add_page()
                pdf.set_font("Arial", size=12)

                pdf.cell(200, 10, txt="Huntington Report", ln=True)
                pdf.cell(200, 10, txt=f"Name: {row['name']}", ln=True)
                pdf.cell(200, 10, txt=f"CAG: {row['cag']}", ln=True)
                pdf.cell(200, 10, txt=f"Risk: {row['risk']}", ln=True)

                pdf.output("report.pdf")

                with open("report.pdf", "rb") as f:
                    st.download_button("Download PDF", f, file_name="report.pdf")

        except Exception as e:
            st.error(str(e))
