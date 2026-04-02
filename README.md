# Stress Concentration & Safety Factor Estimator

A Streamlit web app for estimating stress concentrations and safety factors
in load-bearing mechanical components.

## Files

| File | Purpose |
|------|---------|
| `app.py` | Streamlit web app (UI + charts) |
| `stress_analysis.py` | Core analysis engine (importable module) |
| `requirements.txt` | Python dependencies |

---

## Run Locally

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Launch the app
streamlit run app.py
```

The app opens automatically at http://localhost:8501

---

## Deploy for Free (Streamlit Community Cloud)

1. Push both files to a **public GitHub repository**
2. Go to https://share.streamlit.io and sign in with GitHub
3. Click **New app** → select your repo → set `app.py` as the main file
4. Click **Deploy** — you'll get a shareable URL in ~60 seconds

That's it. The URL is permanent and free.

---

## Features

- **5 stress raiser types** — circular hole, shoulder fillet, U-groove, keyway, press fit
- **7 built-in materials** — steels, aluminium alloys, titanium, cast iron + custom
- **4 load types** — axial, bending, torsion, combined
- **3 failure criteria** — Von Mises, Tresca, Max Normal Stress
- **Interactive gauge charts** for each safety factor
- **Stress breakdown bar chart** comparing nominal → peak → equiv. stress vs Sy / Su
- **Parametric sensitivity plot** — SF vs fillet radius (shoulder fillet feature)
- Automatic design warnings for low safety factors or high Kt

---

## Formulas

Stress concentration factors use Pilkey & Pilkey curve-fit formulas from:
> *Peterson's Stress Concentration Factors*, 3rd edition (2008)

For design verification only. Always consult a qualified mechanical engineer
before making structural decisions.
