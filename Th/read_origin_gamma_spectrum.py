import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
def read_file():
    # --- Step 1: Read the text file ---
    with open("gamma_output_thread0.txt", "r") as f:
        text = f.read()

    # --- Step 2: Extract values using regex ---
    pattern = re.compile(r'^\s*(\d+)\s+([\d.]+)\s*$', re.MULTILINE)

    matches = pattern.findall(text)

    # --- Step 3: Store results in DataFrame ---
    df = pd.DataFrame(matches, columns=['Event', 'Energy/keV']) # in keV



    print(df.head(10))

    for col in df.columns:
        df[col] = df[col].astype(str).str.strip().astype(float)




    print(df.head(20))

    return df

def plot_distribution(df):
    E_list = df["Energy/keV"].to_list()
    event_num = len(df["Event"].unique())
    print("total gamma events", event_num)
    print("total gamma num", len(E_list))
    print("gamma num per event", len(E_list)/event_num)

    fig, axes = plt.subplots()


    # # Select 9 columns to plot (adjust this list as needed)
    # cols = df.columns[:9]

    counts, bins,_=axes.hist(E_list, bins=50, color='steelblue', edgecolor='black', density= True)

    plt.clf()
    fractions = counts / counts.sum()  # fraction = count / total

    energy_MeV = [i*1e-3 for i in bins]

    for i in range(len(fractions)):
        print(str(energy_MeV[i])+" "+str(fractions[i]))


    plt.bar(bins[:-1], fractions, width=np.diff(bins), align='edge')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Fraction of total count")
    plt.show()





if __name__ =="__main__":
    df = read_file()
    plot_distribution(df)