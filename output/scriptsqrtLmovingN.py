import os, fnmatch
import pandas as pd


def findFiles (path, filter):
  for root, dirs, files in os.walk(path):
    for file in fnmatch.filter(files, filter):
      yield os.path.join(root, file)


data = [["n", "t", "L", "time", "type"]]
data2 = []

for textFile in findFiles(r"./", '*.txt'):
  # print(textFile)
  filename = textFile.split("/")[-1].split(".")[0]
  n = int(filename.split("_")[0].split("=")[1])
  t = int(filename.split("_")[1].split("=")[1])
  L = int(filename.split("_")[2].split(".")[0].split("=")[1])
  lines = []
  with open(textFile, "r") as f:
    lines = f.readlines()
  djk = int(float(lines[0].split(" ")[1].strip().replace("+","")))
  astar = int(float(lines[1].split(" ")[1].strip().replace("+","")))
  dp = int(float(lines[2].split(" ")[1].strip().replace("+","")))
  # print(n, t, L, djk, astar, dp)

  data2.append( [n, t, L, djk, "Dijk"] )
  data2.append( [n, t, L, astar, "A*"] )
  data2.append( [n, t, L, dp, "DP"] )

data.append(data2)

df = pd.DataFrame(data=data[1], columns=data[0])

df.to_pickle("resNmovingsqrtL.pkl")
print(df)
