import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
def read_file():
    # --- Step 1: Read the text file ---
    with open("Bi207_Photoneutrons.txt", "r") as f:
        text = f.read()

    # --- Step 2: Extract values using regex ---
    pattern = re.compile(
        r'\b\d+\s+\w+\s+\d+\s+\d+\s+\w+\s+([\d.]+)\s+\(([^)]+)\)\s+\(([^)]+)\)'
    )
    matches = pattern.findall(text)

    # --- Step 3: Store results in DataFrame ---
    df = pd.DataFrame(matches, columns=['Kinetic', 'Pos', 'Vec']) # in keV

    # --- Step 4: Split coordinates into numeric columns ---
    df[['Pos_x', 'Pos_y', 'Pos_z']] = df['Pos'].str.split(',', expand=True)
    df[['Vec_x', 'Vec_y', 'Vec_z']] = df['Vec'].str.split(',', expand=True)

    print(df.head(10))
    df = df.drop(columns=['Pos', 'Vec'])
    # --- Step 5: Clean up and convert to floats ---
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip().astype(float)

    # --- Step 6: Drop original grouped columns ---



    print(df.head(20))

    return df

def plot_distribution(df):
    x_list = df["Pos_x"].to_list()
    y_list = df["Pos_y"].to_list()
    z_list = df["Pos_z"].to_list()
    KE = df["Kinetic"].to_list()
    vec_ux_list = df["Vec_x"].to_list()
    vec_uy_list = df["Vec_y"].to_list()
    vec_uz_list = df["Vec_z"].to_list()

    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    axes = axes.flatten()  # make it easier to loop

    # # Select 9 columns to plot (adjust this list as needed)
    # cols = df.columns[:9]

    counts_x, bins_x,_=axes[0].hist(x_list, bins=50, color='steelblue', edgecolor='black')

    axes[0].set_title("x/mm")
    axes[0].set_xlabel("position/mm")
    axes[0].set_ylabel("Count")

    counts_y, bins_y,_=axes[1].hist(y_list, bins=50, color='steelblue', edgecolor='black')
    axes[1].set_title("y/mm")
    axes[1].set_xlabel("position/mm")
    axes[1].set_ylabel("Count")

    counts_z, bins_z,_=axes[2].hist(z_list, bins=50, color='steelblue', edgecolor='black')
    axes[2].set_title("z/mm")
    axes[2].set_xlabel("position/mm")
    axes[2].set_ylabel("Count")

    axes[3].hist(vec_ux_list, bins=50, color='steelblue', edgecolor='black')
    axes[3].set_title("x unit vector")
    axes[3].set_xlabel("unit vector")
    axes[3].set_ylabel("Count")

    axes[4].hist(vec_uy_list, bins=50, color='steelblue', edgecolor='black')
    axes[4].set_title("y unit vector")
    axes[4].set_xlabel("unit vector")
    axes[4].set_ylabel("Count")

    axes[5].hist(vec_uz_list, bins=50, color='steelblue', edgecolor='black')
    axes[5].set_title("z unit vector")
    axes[5].set_xlabel("unit vector")
    axes[5].set_ylabel("Count")


    counts_ke, bins_ke,_ = axes[6].hist(KE, bins=50, color='steelblue', edgecolor='black')
    axes[6].set_title("ke energy")
    axes[6].set_xlabel("energy/keV")
    axes[6].set_ylabel("Count")


    print("x", len(x_list),counts_x, bins_x)
    print("y", counts_y, bins_y)
    print("z", counts_z, bins_z)
    print("ke",counts_ke, bins_ke)
    bins_ke_updated = (bins_ke[:-1]+bins_ke[1:])*0.001/2 # change to MeV
    for i in range(len(counts_ke)):
        print(bins_ke_updated[i],"MeV",counts_ke[i])


    # first 25 increasing

    slope_x, intercept_x = np.polyfit(bins_x[:25], counts_x[:25], 1)
    count_fited = []
    for i in range(25):
        count_fited.append(bins_x[i]*slope_x+intercept_x)
    axes[0].plot(bins_x[:25], count_fited)
    print("xy size", slope_x, intercept_x)
    z_num = 22
    counts_z_mapped = []
    for i in range(z_num):
        counts_z_mapped.append(1/(counts_z[i])**0.5)

    slope_z, intercept_z = np.polyfit(bins_z[:z_num], counts_z_mapped, 1)
    count_z_fited = []
    for i in range(z_num):
        count_z_fited.append(1/(bins_z[i]*slope_z+intercept_z)**2)
    axes[7].plot(bins_z[:z_num], counts_z_mapped)
    axes[7].set_title("z fit y = 1/(kx+b)**0.5")
    axes[7].set_xlabel("energy/keV")

    axes[2].plot(bins_z[:z_num], count_z_fited)



    print("z size", slope_z, intercept_z)
    print("z mapped", bins_z[:z_num], counts_z_mapped[:z_num])
    print("z fit", count_z_fited[:z_num])



    plt.tight_layout()
    plt.show()


if __name__ =="__main__":
    df = read_file()
    plot_distribution(df)