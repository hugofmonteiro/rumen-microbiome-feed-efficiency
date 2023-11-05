import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from joblib import Parallel, delayed
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score

data = pd.read_csv('data.csv')

X = data.drop('target', axis=1)
y = data.loc[:, ['target']]

def calculate_mse_for_candidate_feature(candidate_feature, selected_features):
    candidate_features = selected_features.union({candidate_feature})
    candidate_X = X.iloc[:, list(candidate_features)]

    mse_scores_fold = []
    y_val_fold = np.array([])
    y_pred_fold = np.array([])
    val_indices = []

    scaler = StandardScaler()

    for train_index, val_index in kfold.split(candidate_X):
        X_train, X_val = candidate_X.iloc[train_index], candidate_X.iloc[val_index]
        y_train, y_val = y[train_index], y[val_index]

        X_train_scaled = scaler.fit_transform(X_train)
        X_val_scaled = scaler.transform(X_val)

        model.fit(X_train_scaled, y_train)
        y_pred = model.predict(X_val_scaled)
        mse = mean_squared_error(y_val, y_pred)

        mse_scores_fold.append(mse)
        y_val_fold = np.concatenate((y_val_fold, y_val))
        y_pred_fold = np.concatenate((y_pred_fold, y_pred))
        val_indices.extend(val_index)

    avg_mse = np.mean(mse_scores_fold)
    return candidate_feature, avg_mse, y_val_fold, y_pred_fold, val_indices

model = Ridge(alpha=1.0)

kfold = KFold(n_splits=10, shuffle=True, random_state=0)

selected_features = set()
mse_scores = []
all_y_val = np.array([])
all_y_pred = np.array([])

best_model = None
best_features = None
lowest_mse = float('inf')

best_val_indices = []

for feature in range(X.shape[1]):
    results = Parallel(n_jobs=-1)(
        delayed(calculate_mse_for_candidate_feature)(candidate_feature, selected_features)
        for candidate_feature in range(X.shape[1])
        if candidate_feature not in selected_features
    )

    best_feature, best_mse, best_y_val, best_y_pred, best_val_indices = min(results, key=lambda x: x[1])

    selected_features.add(best_feature)
    mse_scores.append(best_mse)

    print(f"Selecting feature {best_feature} with MSE: {best_mse:.4f}")

    current_mse = mean_squared_error(best_y_val, best_y_pred)
    if current_mse < lowest_mse:
        lowest_mse = current_mse
        best_model = model
        best_features = selected_features.copy()
        best_val_indices = best_val_indices

X_best_val = X.iloc[best_val_indices, list(best_features)]
y_best_val = y[best_val_indices]

y_pred_best = np.array([])
y_true_best = np.array([])

for train_index, val_index in kfold.split(X_best_val):
    X_train, X_val = X_best_val.iloc[train_index], X_best_val.iloc[val_index]
    y_train, y_val = y_best_val[train_index], y_best_val[val_index]

    best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_val)

    y_pred_best = np.concatenate((y_pred_best, y_pred))
    y_true_best = np.concatenate((y_true_best, y_val))

r2_best = r2_score(y_true_best, y_pred_best)
mse_best = mean_squared_error(y_true_best, y_pred_best)
rmse_best = np.sqrt(mse_best)
print(f"R² based on all folds: {r2_best:.4f}")
print(f"MSE based on all folds: {mse_best:.4f}")
print(f"RMSE based on all folds: {rmse_best:.4f}")

for feature in best_features:
    print(f"Column index: {feature}, Column name: {X.columns[feature]}")

df_pred = pd.DataFrame({'observed_target': y_true_best, 'predicted_target': y_pred_best})

df_pred['Category'] = ['true positive' if (obs > 0) & (pred > 0) else
                       'true negative' if (obs < 0) & (pred < 0) else
                       'false positive' if (obs < 0) & (pred > 0) else
                       'false negative' for obs, pred in zip(df_pred['observed_target'], df_pred['predicted_target'])]

summary = df_pred['Category'].value_counts(normalize=True) * 100
summary
