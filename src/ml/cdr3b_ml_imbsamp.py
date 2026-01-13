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
    precision_recall_curve,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef
)

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import Pipeline as ImbPipeline

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--data_csv", required=True)
parser.add_argument("--output_dir", required=True)

parser.add_argument("--oversample_pos", action="store_true",
                    help="Randomly oversample positive class in training set")
parser.add_argument("--undersample_neg", action="store_true",
                    help="Randomly undersample negative class in training set")
parser.add_argument("--target_pos_frac", type=float, default=0.5,
                    help="Target positive fraction after resampling (default: 0.5)")

args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)


# -----------------------------
# 1. Load data
# -----------------------------
df = pd.read_csv(args.data_csv)
df = df[df["specificity"].isin([0, 1])].reset_index(drop=True)

cdr3_len = df["cdr3"].str.len().unique()
assert len(cdr3_len) == 1
cdr3_len = cdr3_len[0]
print(f"CDR3β length: {cdr3_len}")

y = df["specificity"].values


# -----------------------------
# 2. Position-wise encoding
# -----------------------------
def explode_cdr3(df, column="cdr3"):
    return pd.DataFrame(
        df[column].apply(list).tolist(),
        columns=[f"pos_{i}" for i in range(len(df[column].iloc[0]))]
    )

X_df = explode_cdr3(df)


# -------------------------
# Train / test split
# -------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X_df,
    y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

print(f"Train size (original): {len(y_train)}")
print(f"Positive fraction (original): {y_train.mean():.3f}")


# -------------------------
# Optional resampling (TRAIN ONLY)
# -------------------------
if args.oversample_pos or args.undersample_neg:
    sampling_strategy = args.target_pos_frac

    if args.oversample_pos and args.undersample_neg:
        sampler = ImbPipeline([
            ("over", RandomOverSampler(
                sampling_strategy=sampling_strategy,
                random_state=42
            )),
            ("under", RandomUnderSampler(
                sampling_strategy=sampling_strategy,
                random_state=42
            ))
        ])
    elif args.oversample_pos:
        sampler = RandomOverSampler(
            sampling_strategy=sampling_strategy,
            random_state=42
        )
    elif args.undersample_neg:
        sampler = RandomUnderSampler(
            sampling_strategy=sampling_strategy,
            random_state=42
        )

    X_train, y_train = sampler.fit_resample(X_train, y_train)

    print(f"Train size (resampled): {len(y_train)}")
    print(f"Positive fraction (resampled): {y_train.mean():.3f}")


# -------------------------
# One-hot encoding
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
        ("clf", SVC(
            kernel="linear",
            probability=True,
            class_weight="balanced"
        ))
    ]),
    "SVM_rbf": Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("clf", SVC(
            kernel="rbf",
            probability=True,
            class_weight="balanced"
        ))
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
pr_data = {}

for name, model in models.items():
    print(f"Training {name}...")
    model.fit(X_train_enc, y_train)

    y_prob = model.predict_proba(X_test_enc)[:, 1]
    y_pred = (y_prob >= 0.5).astype(int)

    results.append({
        "model": name,
        "AUROC": roc_auc_score(y_test, y_prob),
        "PR_AUC": average_precision_score(y_test, y_prob),
        "Accuracy": accuracy_score(y_test, y_pred),
        "Precision": precision_score(y_test, y_pred, zero_division=0),
        "Recall": recall_score(y_test, y_pred),
        "F1": f1_score(y_test, y_pred),
        "MCC": matthews_corrcoef(y_test, y_pred)
    })

    roc_data[name] = roc_curve(y_test, y_prob)
    pr_data[name] = precision_recall_curve(y_test, y_prob)


# -------------------------
# PyTorch MLP (NO class weighting)
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

mlp.train()
for epoch in range(15):
    for xb, yb in train_loader:
        xb, yb = xb.to(device), yb.to(device)
        opt.zero_grad()
        loss = loss_fn(mlp(xb), yb)
        loss.backward()
        opt.step()

mlp.eval()
with torch.no_grad():
    logits = mlp(X_test_dense.to(device))
    probs = torch.sigmoid(logits).cpu().numpy()

preds = (probs >= 0.5).astype(int)

results.append({
    "model": "MLP_2layer",
    "AUROC": roc_auc_score(y_test, probs),
    "PR_AUC": average_precision_score(y_test, probs),
    "Accuracy": accuracy_score(y_test, preds),
    "Precision": precision_score(y_test, preds, zero_division=0),
    "Recall": recall_score(y_test, preds),
    "F1": f1_score(y_test, preds),
    "MCC": matthews_corrcoef(y_test, preds)
})

roc_data["MLP_2layer"] = roc_curve(y_test, probs)
pr_data["MLP_2layer"] = precision_recall_curve(y_test, probs)


# -------------------------
# Save metrics
# -------------------------
res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(args.output_dir, "metrics_ml_imb.csv"), index=False)


# -------------------------
# Plot metrics
# -------------------------
plt.figure(figsize=(10, 5))
res_df.set_index("model")[["AUROC", "PR_AUC", "MCC"]].plot(kind="bar")
plt.ylabel("Score")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "metrics_barplot_imb.png"))
plt.close()


# -------------------------
# Plot ROC curves
# -------------------------
plt.figure(figsize=(6, 6))
for name, (fpr, tpr, _) in roc_data.items():
    plt.plot(fpr, tpr, label=name)

plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "roc_curves_imb.png"))
plt.close()


# -------------------------
# Plot PR curves
# -------------------------
plt.figure(figsize=(6, 6))
for name, (precision, recall, _) in pr_data.items():
    ap = res_df.loc[res_df.model == name, "PR_AUC"].values[0]
    plt.plot(recall, precision, label=f"{name} (AP={ap:.2f})")

baseline = y_test.mean()
plt.hlines(baseline, 0, 1, linestyles="dashed", colors="gray", label="Random")

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "pr_curves_imb.png"))
plt.close()

print("Analysis complete.")