import numpy as np
import pandas as pd
import sklearn
import random

output_name = "data/custom_abalone.data"
generate_size = 50000

if(len(sys.argv) != 1):
    if(len(sys.argv) != 3):
        print("usage: %s <data_location=data/custom_abalone.data> <data_size=50000>", sys.argv[0])
        return
    else:
        output_name = sys.argv[1]
        generate_size = (int)(sys.argv[2])


df = pd.read_csv("data/abalone.data", delimiter=" ", header=None)

df.iloc[0]

min_max_collection = df.agg(['min', 'max'])

def generator(df, max_size):
    for i in range(0, max_size):
        rand_num = random.random()
        single_result = []
        if(rand_num < 0.5):
            gender = 'F'
        else:
            gender = 'M'
        single_result.append(gender)
        for j in range(1, len(df.iloc[0])):
            single_result.append(random.uniform(min_max_collection[j][0] * 0.55, min_max_collection[j][1] * 1.18))
        df.loc[len(df.index)] = single_result
    return df
        
new_df = generator(df, generate_size)
new_df.to_csv(output_name, sep=' ',index = False,header = False)
