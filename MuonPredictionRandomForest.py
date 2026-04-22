#%% import functions
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt

#%% Data Preparation

def load_data_as_df (p, e, sin2theta):

    # file path
    fMuon = f'GroundTotalParticles/{p}_{e}_{sin2theta}.npz'
    fXmax = f'Xmax/{p}_{e}_sin2_{sin2theta}.dat'
    fRadE = f'radEnergy/norm_sintheta2/{p}_{e}_{sin2theta}.npz'
    fEdep = f'TotalEdepScint/{p}_{e}_{sin2theta}.npz'
    
    # load data
    fileRadE = np.load(fRadE, allow_pickle=True)
    fileMuon = np.load(fMuon, allow_pickle= True)
    fileEdep = np.load(fEdep, allow_pickle=True)
    fileXmax = np.loadtxt(fXmax)

    # determine parameters
    energy = np.log(fileEdep["primaryE"]) # logE, E is primary particle energy in GeV
    Xmax  = fileXmax[:,6] # Xmax 
    Nep_Xmax  = np.log10(fileXmax[:,3]) # number of total +-e at Xmax
    Nmu_ground = np.log10(fileMuon['nMu']) # number of total +-mu at ground
    zenith = fileMuon['zenith'] # zenith angle in rad
    Edep = np.log10(fileEdep["Edep_tot"])
    RadE = np.log10(fileRadE['radE_filtered(eV)'])
    alpha = fileRadE['alpha']

    cosalpha = np.cos(alpha)
    costheta = np.cos(zenith)
    mask = Xmax >= 0 # keep only positive Xmax

    # store all parameters in data frame
    df = pd.DataFrame({
        'particle': p,                          # primary particle
        'energy': energy[mask],                 # plog(E) primary energy
        'costheta': costheta[mask],             # sin2(zenith)
        'Xmax': Xmax[mask],                     # Xmax
        'cosalpha': cosalpha[mask],             # cos(alpha), alpha = angle between shower and magntic field
        'Edep': Edep[mask],                     # log10(Edep) deposited energy
        'Erad': RadE[mask],                     # log10(Erad)radiation energy
        'Nmu': Nmu_ground[mask],                 # log10(Nmu) number of muon at ground
        'Ne' : Nep_Xmax[mask]                   # log10(Ne+-) number of electron and positron at Xmax
    })
    
    return df

particles = ['proton', 'iron', 'helium', 'oxygen']
energies =  [f"lgE_{x/10:.1f}" for x in range(160,181)]
sin2thetas = [x/10 for x in range(0, 10)]

all_df = []
print('starting to download files')
for p in particles:
    for e in energies:
        for z in sin2thetas:
            try:
                df_bin = load_data_as_df (p, e, z)
                all_df.append(df_bin)
            except FileNotFoundError:
                continue 

# gather all bins into one big data frame
final_df = pd.concat(all_df, ignore_index=True)
final_df.replace([np.inf, -np.inf], np.nan, inplace=True) #replace inf values with nan
final_df.dropna(inplace=True) # remove all nans
final_df.to_parquet('data_for_regression.parquet', index=False)
print(final_df.head())
print(f"Total rows (events): {len(final_df)}/168,000")
print('file is saved')

#%% Regression 

df = pd.read_parquet('data_for_regression.parquet')
# print(df.head())

# label decoder changes string to index
le_particle = LabelEncoder()
df['particle'] = le_particle.fit_transform(df['particle'])

# Edep, Erad, and Nmu are in log10
inputs = ['costheta', 'Xmax', 'cosalpha', 'Edep', 'Erad']
outputs = ['Nmu', 'Ne']
X = df[inputs]
y = df[outputs]

print('Input')
print(X[:5])
print('Output')
print(y[:5])


# splitting Train set (50%) and Test set (50%)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)

# model training
model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# prediction
y_pred = model.predict(X_test)

# result evaluation
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Total events: {len(final_df)}/168,000")
print(f"Train events: {len(X_train)}")
print(f"Test events: {len(X_test)}")

print(f"Mean Squared Error: {mse:.4f}")
print(f"R-squared Score: {r2:.4f}")

# Feature importance
importances = pd.Series(model.feature_importances_, index=inputs)
print("\nFeature Importances:")
print(importances.sort_values(ascending=False))


#%% Plotting

plt.scatter(y_test, y_pred, s = 0.5)
plt.xlabel(f'True log($N_{{\mu}}$)')
plt.ylabel(f'Predicted log($N_{{\mu}}$)')
plt.text(9, 15.5, f"MSE = {mse:.4f}" + "\n" + f"$R^2$ = {r2:.4f}")


# %% Muon component prediction

new_input = pd.DataFrame({ 
    'costheta': [0.0], 
    'Xmax': [650.0], 
    'cosalpha': [-0.95], 
    'Edep': [6.70], 
    'Erad': [4.2]
})


print(f"Input parameters: " \
    f"{new_input}")

new_input['particle'] = le_particle.transform(new_input['particle'])

prediction = model.predict(new_input)

print(f"predicted muon number: {prediction[0]:.4f}")

# %%
