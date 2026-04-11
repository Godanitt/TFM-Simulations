import dill 

# IR_yields

with open("IR_yields.pkl", "rb") as f:
    df = dill.load(f)
