import subprocess
fillretinotopy_script = "/data/p_02915/SPOT/run_retinotopy.py"

# Iterate over the index array using a for loop
for index in range(192, 193): 

    cmd_fillretinotopy = [
        "python",
        fillretinotopy_script,
        str(index)
    ]
    #print(cmd_fillretinotopy)
    subprocess.run(cmd_fillretinotopy, check=True)