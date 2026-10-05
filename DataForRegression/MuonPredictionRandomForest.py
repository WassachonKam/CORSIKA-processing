#%% import functions
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import joblib

#%% Regression 

df = pd.read_parquet('data_for_regression_filteredXmax.parquet')
# print(df.head())

# label decoder changes string to index
le_particle = LabelEncoder()
df['particle'] = le_particle.fit_transform(df['particle'])

# splitting Train set (50%) and Test set (50%)
train_df, test_df = train_test_split(df, test_size= 0.5, random_state= 42)

# Edep, Erad, and Nmu are in log10
inputs = ['costheta', 'Xmax', 'sinalpha', 'Edep', 'Erad']
# outputs = ['Nmu', 'Ne']
outputs = ['Nmu']

# store data in dictionary
data_dict = {
    "train": {
        "X": train_df[inputs],
        "y": train_df[outputs],
    },
    "test": {
        "X": test_df[inputs],
        "y": test_df[outputs],
    }
}

X_train = data_dict['train']['X']
y_train = data_dict['train']['y']
X_test = data_dict['test']['X']
y_test = data_dict['test']['y']

print('Input')
print(X_train.head())
print('Output')
print(y_train.head())

# model training
# model = RandomForestRegressor(n_estimators=200, criterion= "absolute_error", random_state=42)
model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# prediction
y_pred = model.predict(X_test)

# result evaluation
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Total events: {len(df['particle'])}/168,000")
print(f"Train events: {len(X_train)}")
print(f"Test events: {len(X_test)}")

print(f"Mean Squared Error: {mse:.4f}")
print(f"R-squared Score: {r2:.4f}")

# Feature importance
importances = pd.Series(model.feature_importances_, index=inputs)
print("\nFeature Importances:")
print(importances.sort_values(ascending=False))

# Save data 
save_dict = {
    "model": model,
    "test_data": {
        "X": X_test,
        "y": y_test,
        "energy": test_df['energy'],
        "particle_id": test_df['particle'] ,
        "particle_str": le_particle.classes_ , # particle ID sorted alphabetically from 0, 1, 2,...
        "description": {
        # "description": "Muon and EM Reconstruction using RandomForestRegressor",
        # "description": "RandomForestRegressor with absolute error",
        "description": "filtered Xmax with squared error",
        "date": "2026-04-23"
    }
    }
}
# joblib.dump(save_dict, "muon_prediction_model.pkl")
# joblib.dump(save_dict, "muon_em_prediction_model.pkl")
# joblib.dump(save_dict, "muon_prediction_model_aberror.pkl")
joblib.dump(save_dict, "muon_prediction_model_sqerror_filtered.pkl")
#%% Plotting
y_test = np.array(y_test)
y_pred = np.array(y_pred)
energy = np.array(test_df['energy'] )
mse_mu, mse_em = mean_squared_error(y_test[:,0], y_pred[:,0]), mean_squared_error(y_test[:,1], y_pred[:,1])
r2_mu, r2_em = r2_score(y_test[:,0], y_pred[:,0]), r2_score(y_test[:,1], y_pred[:,1])

plt.scatter(y_test[:,0], y_pred[:,0], s = 0.5, c=energy, cmap="viridis")
plt.xlabel(f'True log($N_{{\mu^{{\pm}}}}$)')
plt.ylabel(f'Predicted log($N_{{\mu^{{\pm}}}}$)')
plt.text(4, 6.5, f"MSE = {mse:.4f}" + "\n" + f"$R^2$ = {r2:.4f}")
plt.title(rf'$\mu^{{\pm}}$')
plt.colorbar(label=rf"log($E_{{primary}}$)")
plt.show()

plt.scatter(y_test[:,1], y_pred[:,1], s = 0.5, c=energy, cmap="viridis")
plt.xlabel(f'True log($N_{{e^{{\pm}}}}$)')
plt.ylabel(f'Predicted log($N_{{e^{{\pm}}}}$)')
plt.text(6.7, 8.5, f"MSE = {mse_em:.4f}" + "\n" + f"$R^2$ = {r2_em:.4f}")
plt.title(rf'$e^{{\pm}}$')
plt.colorbar(label=rf"log($E_{{primary}}$)")


# %% Muon component prediction

new_input = pd.DataFrame({ 
    'costheta': [0.0], 
    'Xmax': [650.0], 
    'sinalpha': [-0.95], 
    'Edep': [6.70], 
    'Erad': [4.2]
})


print(f"Input parameters: " \
    f"{new_input}")

new_input['particle'] = le_particle.transform(new_input['particle'])

prediction = model.predict(new_input)

print(f"predicted muon number: {prediction[0]:.4f}")

# %%
