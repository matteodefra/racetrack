import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# fig, ax = plt.subplots( 1, 1, figsize=(15, 15))
df = pd.read_pickle("resLmovingsqrtL.pkl")

df = df[df["type"] != "Dijk"]


sns.lineplot(data=df, x='n', y='time', hue='type')#.set_title("10<=L<=25")

plt.show()
