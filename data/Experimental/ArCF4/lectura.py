import dill 


with open("IR_yields.pkl", "rb") as f:
    df = dill.load(f)
