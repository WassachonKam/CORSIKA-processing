#%% import functions
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import joblib

#%% Data Preparation

def load_data_as_df (p, e, sin2theta, threshR, Xmaxfiltering, NmuNorm):

    # file path
    fMuon = f'../GroundTotalParticles/{p}_{e}_{sin2theta}.npz'
    fXmax = f'../Xmax/{p}_{e}_sin2_{sin2theta}.dat'
    fRadE = f'../radEnergy/norm_sintheta2/{p}_{e}_{sin2theta}.npz'
    fEdep = f'../TotalEdepScint/{threshR}m/{p}_{e}_{sin2theta}.npz'
    
    # load data
    fileRadE = np.load(fRadE, allow_pickle=True)
    fileMuon = np.load(fMuon, allow_pickle= True)
    fileEdep = np.load(fEdep, allow_pickle=True)
    fileXmax = np.loadtxt(fXmax)

    # determine parameters
    # log10 version
    # energy = np.log10(fileEdep["primaryE"]) # logE, E is primary particle energy in GeV
    # Xmax  = fileXmax[:,6] # Xmax 
    # Nep_Xmax  = np.log10(fileXmax[:,5]) # number of total +-e at Xmax
    # Ne_ground = np.log10(fileMuon['nEP']) # number of total +-mu at ground
    # Nmu_ground = np.log10(fileMuon['nMu']) # number of total +-mu at ground
    # zenith = fileMuon['zenith'] # zenith angle in rad
    # Edep = np.log10(fileEdep["Edep_tot"])
    # RadE = np.log10(fileRadE['radE_filtered(eV)'])
    # alpha = fileRadE['alpha']

    energy = fileEdep["primaryE"] # logE, E is primary particle energy in GeV
    Xmax  = fileXmax[:,6] # Xmax 
    Nep_Xmax  = fileXmax[:,5] # number of total +-e at Xmax
    Ne_ground = fileMuon['nEP'] # number of total +-mu at ground
    Nmu_ground = fileMuon['nMu'] # number of total +-mu at ground
    zenith = fileMuon['zenith'] # zenith angle in rad
    Edep = fileEdep["Edep_tot"]/1000 # Edep in GeV
    RadE = fileRadE['radE_filtered(eV)']/1e9 #Erad in GeV
    alpha = fileRadE['alpha']

    sinalpha = np.sin(alpha)
    costheta = np.cos(zenith)

    X0 = 697.6
    Xv = X0/ costheta

    mask = Xmax >= 0 # keep only positive Xmax
    if Xmaxfiltering == True:
        mask = (Xmax >=0) & (Xmax <= Xv)

    # store all parameters in data frame
    df = pd.DataFrame({
        'particle': p,                          # primary particle
        'energy': energy[mask],                 # plog(E) primary energy
        'costheta': costheta[mask],             # sin2(zenith)
        'Xmax': Xmax[mask],                     # Xmax
        'sinalpha': sinalpha[mask],             # cos(alpha), alpha = angle between shower and magntic field
        'Edep': Edep[mask],                     # log10(Edep) deposited energy
        'Erad': RadE[mask],                     # log10(Erad)radiation energy
        'Nmu': Nmu_ground[mask],                # log10(Nmu) number of muon at ground
        'Ne' : Nep_Xmax[mask],                  # log10(Ne+-) number of electron and positron at Xmax
        'Ne_ground': Ne_ground[mask]            # log10(Ne+-) number of electron and positron at ground    
    })
    
    return df


threshR = 0
Xmaxfiltering = False
NmuNorm = False

particles = ['proton', 'iron', 'helium', 'oxygen']
energies =  [f"lgE_{x/10:.1f}" for x in range(160,181)]
sin2thetas = [x/10 for x in range(0, 10)]

all_df = []
print('starting to download files')
for p in particles:
    for e in energies:
        for z in sin2thetas:
            try:
                df_bin = load_data_as_df (p, e, z, threshR, Xmaxfiltering, NmuNorm)
                all_df.append(df_bin)
            except FileNotFoundError:
                continue 

# gather all bins into one big data frame
final_df = pd.concat(all_df, ignore_index=True)
final_df.replace([np.inf, -np.inf], np.nan, inplace=True) #replace inf values with nan
final_df.dropna(inplace=True) # remove all nans
final_df.to_parquet(f'data_for_regression_{threshR}_Xmax_{Xmaxfiltering}_NmuNorm_{NmuNorm}.parquet', index=False)
print(final_df.head())
print(f"Total rows (events): {len(final_df)}/168,000")
print('file is saved')

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
