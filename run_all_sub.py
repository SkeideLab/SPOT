import subprocess
fillretinotopy_script = "/data/p_02915/SPOT/01_dataprep/simulation_alternativebolddata/simulated_model.py"

# Iterate over the index array using a for loop
for index in range(1, 2): 

    cmd_fillretinotopy = [
        "python",
        fillretinotopy_script,
        str(index)
    ]
    #print(cmd_fillretinotopy)
    subprocess.run(cmd_fillretinotopy, check=True)