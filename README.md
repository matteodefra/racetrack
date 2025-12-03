# racetrack

To compile, run 

```make all```

You will find the executable in the output folder.

Now, decomment lines 181-186 of fulltest.cpp to store the output files. Successively, move into output 

```cd output```

At this stage, you can now run the different bash script one by one to collect result.

Launch the python script script.py to merge all the txt result files in a pandas dataframe. Then, remove all txt files to polish

```rm -rf *.txt```

And finally, launch the script read.py to show the according plot.

You can uncomment line 107 of fulltest.cpp and comment line 111 to test the speed max with \sqrt{L} instead of S/4.


Below are all the steps detailed to reproduce each experiment:

1. Uncomment lines 181-186 of fulltest.cpp

## Experiment S/4, L moving

To reproduce this experiment

1. Launch the script ```sh LmovingS4.sh```
2. Move into output ```cd output```
3. Launch the python script ```python3 scriptS4movingL.py```
4. Remove all txt files ```rm -rf *.txt```
5. Finally, launch ```python3 readS4movingL.py```

## Experiment S/4, n moving

To reproduce this experiment

1. Launch the script ```sh NmovingS4.sh```
2. Move into output ```cd output```
3. Launch the python script ```python3 scriptS4movingN.py```
4. Remove all txt files ```rm -rf *.txt```
5. Finally, launch ```python3 readS4movingN.py```

For the experiments with sqrt{L}, uncomment line 107 of fulltest.cpp and comment line 111 to test the speed max with \sqrt{L} instead of S/4.

## Experiment sqrt{L}, L moving

To reproduce this experiment

1. Launch the script ```sh LmovingsqrtL.sh```
2. Move into output ```cd output```
3. Launch the python script ```python3 scriptsqrtLmovingL.py```
4. Remove all txt files ```rm -rf *.txt```
5. Finally, launch ```python3 readsqrtLmovingL.py```

## Experiment sqrt{L}, n moving

To reproduce this experiment

1. Launch the script ```sh NmovingsqrtL.sh```
2. Move into output ```cd output```
3. Launch the python script ```python3 scriptsqrtLmovingN.py```
4. Remove all txt files ```rm -rf *.txt```
5. Finally, launch ```python3 readsqrtLmovingN.py```
