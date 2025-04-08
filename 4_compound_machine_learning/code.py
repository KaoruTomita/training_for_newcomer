from typing import List, Union
import numpy.typing as npt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import root_mean_squared_error, make_scorer
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Draw, Descriptors



def draw_molecule(csvfile: str) -> None:
    # 課題 4-1
    df=pd.read_csv(csvfile)
    df["smiles_clean"] = df["SMILES"].astype(str).str.strip()
    target = df[df["Compound ID"] == "CHEMBL540227"]
    smiles = target["smiles_clean"].values[0]
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, legend="CHEMBL540227")
    img.save("CHEMBL540227.png")
    pass

def create_2d_descriptors(smiles: str) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 4-2
    mol= Chem.MolFromSmiles(smiles)
    desc_list = Descriptors.descList
    descriptor_names = [desc_name for desc_name, _ in desc_list]
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    descriptor_values = list(calculator.CalcDescriptors(mol))
    return descriptor_values

def predict_logpapp(csvfile: str) -> Union[npt.NDArray[np.float_], pd.Series, List[float]]:
    # 課題 4-3
    np.random.seed(0) # 出力を固定するためにseedを指定
    rfr = RandomForestRegressor(random_state=0) # 出力を固定するためにrandom_stateを指定
    df=pd.read_csv(csvfile)
    df["smiles_clean"] = df["SMILES"].astype(str).str.strip()
    y = df["LogP app"]
    df["2D_descriptors"] = df["smiles_clean"].apply(create_2d_descriptors)
    X = pd.DataFrame(df["2D_descriptors"].tolist(), index=df.index)
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    rfr.fit(X_train, y_train)
    predicteds = rfr.predict(X_test)
    

    return predicteds.tolist()

def grid_search(csvfile: str) -> float:
    # 課題 4-4
    # こちらも出力を固定するためにseedやrandom_stateを指定すること
    np.random.seed(0)
    rfr = RandomForestRegressor(random_state=0)
    df=pd.read_csv(csvfile)
    df["smiles_clean"] = df["SMILES"].astype(str).str.strip()
    df["2D_descriptors"] = df["smiles_clean"].apply(create_2d_descriptors)
    X = pd.DataFrame(df["2D_descriptors"].tolist(), index=df.index)
    y = df["LogP app"]
    param_grid = {"n_estimators": [100, 200, 400], "max_depth": [5, 10, 15]}
    # # Xは説明変数、yは目的変数
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    grid_search = GridSearchCV(rfr, param_grid, cv=4, scoring=make_scorer(root_mean_squared_error,greater_is_better=False))
    grid_search.fit(X_train, y_train)
    y_pred = grid_search.predict(X_test)

    return root_mean_squared_error(y_test, y_pred)

if __name__ == "__main__":
    smiles = "C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1"
    filepath = "data/fukunishi_data.csv"
    # 課題 4-1
    draw_molecule(filepath)
    # 課題 4-2
    print(create_2d_descriptors(smiles))
    # 課題 4-3
    print(predict_logpapp(filepath))
    # 課題 4-4
    print(grid_search(filepath))
