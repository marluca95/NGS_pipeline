import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef, 
    precision_recall_curve,
    average_precision_socore
)

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier



import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--data_csv", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)


# -----------------------------
# 1. Load data
# -----------------------------
df = pd.read_csv(args.data_csv)
df = df[df["specificity"].isin([0, 1])].reset_index(drop=True)
cdr3_len = df["cdr3"].str.len().unique()[0]
print(f"CDR3β length: {cdr3_len}")

# -----------------------------
# 2. Position-wise encoding
# -----------------------------
def explode_cdr3(df, column="cdr3", prefix="pos"):
    return pd.DataFrame(df[column].apply(list).tolist(),
                        columns=[f"{prefix}_{i}" for i in range(df[column].str.len().iloc[0])])

X_df = explode_cdr3(df)
y = df["specificity"].values


# -------------------------
# Train / test split
# -------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X_df, y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

# -------------------------
# One-hot encoding (sparse)
# -------------------------
ohe = OneHotEncoder(
    handle_unknown="ignore",
    sparse_output=True,
    dtype=np.uint8
)

X_train_enc = ohe.fit_transform(X_train)
X_test_enc = ohe.transform(X_test)


# -------------------------
# Models (sklearn)
# -------------------------
models = {
    "SVM_linear": Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("clf", SVC(kernel="linear", probability=True, class_weight="balanced"))
    ]),
    "SVM_rbf": Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("clf", SVC(kernel="rbf", probability=True, class_weight="balanced"))
    ]),
    "LogReg": Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("clf", LogisticRegression(
            max_iter=500,
            class_weight="balanced",
            n_jobs=4
        ))
    ]),
    "RandomForest": RandomForestClassifier(
        n_estimators=300,
        max_depth=15,
        min_samples_leaf=5,
        n_jobs=4,
        class_weight="balanced",
        random_state=42
    )
}


# -------------------------
# Train + evaluate sklearn models
# -------------------------
results = []
roc_data = {}


for name, model in models.items():
    print(f"Training {name}...")
    model.fit(X_train_enc, y_train)

    y_prob = model.predict_proba(X_test_enc)[:, 1]
    y_pred = (y_prob >= 0.5).astype(int)

    results.append({
        "model": name,
        "AUROC": roc_auc_score(y_test, y_prob),
        "Accuracy": accuracy_score(y_test, y_pred),
        "Precision": precision_score(y_test, y_pred, zero_division=0),
        "Recall": recall_score(y_test, y_pred),
        "F1": f1_score(y_test, y_pred),
        "MCC": matthews_corrcoef(y_test, y_pred)
    })

    fpr, tpr, _ = roc_curve(y_test, y_prob)
    roc_data[name] = (fpr, tpr)


# -------------------------
# MLP (2 hidden layers, GPU if available)
# -------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("MLP using device:", device)

X_train_dense = torch.tensor(X_train_enc.toarray(), dtype=torch.float32)
X_test_dense = torch.tensor(X_test_enc.toarray(), dtype=torch.float32)
y_train_t = torch.tensor(y_train, dtype=torch.float32)
y_test_t = torch.tensor(y_test, dtype=torch.float32)

train_loader = DataLoader(
    TensorDataset(X_train_dense, y_train_t),
    batch_size=256,
    shuffle=True
)

class MLP(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)
        )

    def forward(self, x):
        return self.net(x).squeeze(1)


mlp = MLP(X_train_dense.shape[1]).to(device)
opt = torch.optim.Adam(mlp.parameters(), lr=1e-3)
loss_fn = nn.BCEWithLogitsLoss()

# Train
mlp.train()
for epoch in range(15):
    for xb, yb in train_loader:
        xb, yb = xb.to(device), yb.to(device)
        opt.zero_grad()
        loss = loss_fn(mlp(xb), yb)
        loss.backward()
        opt.step()

# Evaluate
mlp.eval()
with torch.no_grad():
    logits = mlp(X_test_dense.to(device))
    probs = torch.sigmoid(logits).cpu().numpy()

preds = (probs >= 0.5).astype(int)

results.append({
    "model": "MLP_2layer",
    "AUROC": roc_auc_score(y_test, probs),
    "Accuracy": accuracy_score(y_test, preds),
    "Precision": precision_score(y_test, preds, zero_division=0),
    "Recall": recall_score(y_test, preds),
    "F1": f1_score(y_test, preds),
    "MCC": matthews_corrcoef(y_test, preds)
})

fpr, tpr, _ = roc_curve(y_test, probs)
roc_data["MLP_2layer"] = (fpr, tpr)


# -------------------------
# Save results
# -------------------------
res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(args.output_dir, "metrics.csv"), index=False)


# -------------------------
# Plot: metrics
# -------------------------
plt.figure(figsize=(10, 5))
res_df.set_index("model")[["AUROC", "MCC", "Accuracy", "Precision", "Recall", "F1"]].plot(kind="bar")
plt.ylabel("Score")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "metrics_barplot.png"))
plt.close()


# -------------------------
# Plot: ROC curves
# -------------------------
plt.figure(figsize=(6, 6))
for name, (fpr, tpr) in roc_data.items():
    plt.plot(fpr, tpr, label=name)

plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "roc_curves.png"))
plt.close()

print("Analysis complete.")
