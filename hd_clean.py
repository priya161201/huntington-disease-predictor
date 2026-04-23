import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import hashlib
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
import plotly.express as px
from fpdf import FPDF
import tempfile

# -----------------------------
# DATABASE
# -----------------------------
conn = sqlite3.connect("app.db", check_same_thread=False)
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS users (
    username TEXT PRIMARY KEY,
    password TEXT
)
""")

c.execute("""
CREATE TABLE IF NOT EXISTS patients (
    sample TEXT,
    prediction TEXT,
    risk_score REAL,
    risk_level TEXT
)
""")

conn.commit()

# -----------------------------
# AUTH
# -----------------------------
def hash_password(password):
    return hashlib.sha256(password.encode()).hexdigest()

def signup():
    st.subheader("Sign Up")
    user = st.text_input("Username")
    pwd = st.text_input("Password", type="password")

    if st.button("Create Account"):
        try:
            c.execute("INSERT INTO users VALUES (?,?)", (user, hash_password(pwd)))
            conn.commit()
            st.success("Account created! Please login.")
        except:
            st.error("User already exists")

def login():
    st.subheader("Login")
    user = st.text_input("Username")
    pwd = st.text_input("Password", type="password")

    if st.button("Login"):
        c.execute("SELECT * FROM users WHERE username=? AND password=?",
                  (user, hash_password(pwd)))
        if c.fetchone():
            st.session_state["login"] = True
            st.session_state["user"] = user
        else:
            st.error("Invalid credentials")

def logout():
    st.session_state["login"] = False

# -----------------------------
# MAIN APP
# -----------------------------
def main_app():

    st.title("🧬 Huntington Disease Risk Prediction Platform")

    st.sidebar.success(f"Logged in as {st.session_state['user']}")
    if st.sidebar.button("Logout"):
        logout()

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Prediction",
        "Patient Database",
        "Genome Visualization",
        "Model Analysis",
        "Medical Report"
    ])

    # =============================
    # TAB 1: PREDICTION
    # =============================
    with tab1:

        gene_file = st.file_uploader("Upload Gene Expression File", type=["csv"])
        clinical_file = st.file_uploader("Upload Clinical Data File", type=["csv"])

        if gene_file is not None:

            try:
                gene_df = pd.read_csv(gene_file)
                gene_df.set_index(gene_df.columns[0], inplace=True)

                log_df = np.log2(gene_df.T + 1)
                log_df.index = log_df.index.astype(str).str.strip()

                if clinical_file is not None:
                    clinical = pd.read_csv(clinical_file)
                    clinical.columns = clinical.columns.str.strip()

                    if len(clinical.columns) == 1:
                        clinical = clinical.iloc[:, 0].str.split(",", expand=True)
                        clinical.columns = [
                            "Sample","Group","Age","Gender","PMI","RIN","Disease_Stage"
                        ]

                    if "Group" not in clinical.columns:
                        st.error("Missing 'Group' column")
                        st.stop()

                    clinical["Sample"] = clinical["Sample"].astype(str).str.strip()
                    clinical = clinical[clinical["Sample"].isin(log_df.index)]
                    clinical = clinical.set_index("Sample")
                    clinical = clinical.reindex(log_df.index)

                    labels = clinical["Group"]
                    labels = labels.dropna()

                    log_df = log_df.loc[labels.index]
                    labels = labels.astype(str).values

                else:
                    n = len(log_df)
                    labels = np.array(['Control']*(n//2)+['Disease']*(n-n//2))

                if len(labels) != len(log_df):
                    st.error("Mismatch between data")
                    st.stop()

                X = log_df.values

                X_train, X_test, y_train, y_test = train_test_split(
                    X, labels, test_size=0.2, random_state=42
                )

                model = LogisticRegression(max_iter=1000)
                model.fit(X_train, y_train)

                preds = model.predict(X)
                probs = model.predict_proba(X)[:,1]

                risk = ["Low" if r<0.3 else "Medium" if r<0.7 else "High" for r in probs]

                results = pd.DataFrame({
                    "Sample": log_df.index,
                    "Prediction": preds,
                    "Risk Score": probs,
                    "Risk Level": risk
                })

                results = results.sort_values(
                    by="Sample",
                    key=lambda x: x.str.extract(r'(\d+)').astype(int)[0]
                ).reset_index(drop=True)

                st.dataframe(results)

                if st.button("Save Results"):
                    for _, row in results.iterrows():
                        c.execute("INSERT INTO patients VALUES (?,?,?,?)",
                                  (row["Sample"], row["Prediction"],
                                   float(row["Risk Score"]), row["Risk Level"]))
                    conn.commit()
                    st.success("Saved!")

                st.session_state["results"] = results

            except Exception as e:
                st.error(e)

    # =============================
    # TAB 2: DATABASE
    # =============================
    with tab2:
        df = pd.read_sql("SELECT * FROM patients", conn)
        st.dataframe(df)

    # =============================
    # TAB 3: VISUALIZATION
    # =============================
    with tab3:
        if "results" in st.session_state:
            fig = px.bar(st.session_state["results"],
                         x="Sample", y="Risk Score")
            st.plotly_chart(fig)

    # =============================
    # TAB 4: MODEL ANALYSIS
    # =============================
    with tab4:
        if "results" in st.session_state:
            df = st.session_state["results"]
            st.metric("Total Samples", len(df))
            st.metric("High Risk", sum(df["Risk Level"]=="High"))

    # =============================
    # TAB 5: REPORT
    # =============================
    with tab5:
        if "results" in st.session_state:

            df = st.session_state["results"]
            sample = st.selectbox("Select Sample", df["Sample"])
            p = df[df["Sample"] == sample].iloc[0]

            report = f"""
Huntington Disease Report

Sample: {sample}
Prediction: {p['Prediction']}
Risk Score: {round(p['Risk Score'],3)}
Risk Level: {p['Risk Level']}
"""

            st.text(report)

            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Arial", size=12)

            for line in report.split("\n"):
                pdf.cell(200, 10, txt=line, ln=True)

            temp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdf")
            pdf.output(temp.name)

            with open(temp.name, "rb") as f:
                st.download_button("Download PDF", f,
                                   file_name=f"{sample}_report.pdf")

# -----------------------------
# FLOW
# -----------------------------
if "login" not in st.session_state:
    st.session_state["login"] = False

menu = ["Login", "Sign Up"]
choice = st.sidebar.selectbox("Menu", menu)

if not st.session_state["login"]:
    if choice == "Login":
        login()
    else:
        signup()
else:
    main_app()
