import dill 


with open("CF4_primary_data_final.pkl", "rb") as f:
    df = dill.load(f)
