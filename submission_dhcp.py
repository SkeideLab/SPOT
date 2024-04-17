from pathlib import Path
from simple_slurm import Slurm

def submit_job(
    args_list,
    cpus=8,
    mem=32000,
    time="72:00:00",
    log_dir="logs/",
    dependency_jobs=[],
    dependency_type="afterok",
    job_name="job",
):
    """Submits a single batch job via SLURM, which can depend on other jobs.

    Parameters
    ----------
    args_list : list
        A list of shell commands and arguments. The first element will usually
        be the path of a shell script and the following elements the input
        arguments to this script.
    cpus : int, default=8
        The number of CPUs that the batch job should use.
    mem : int, default=320000
        The amount of memory (in MB) that the abtch job should use.
    time : str, default='24:00:00'
        The maximum run time (in format 'HH:MM:SS') that the batch job can use.
        Must exceed 24 hours.
    log_dir : str or Path, default='logs/'
        Directory to which the standard error and output messages of the batch
        job should be written.
    dependency_jobs : int or list, default=[]
        Other SLURM batch job IDs on which the current job depends. Can be used
        to create a pipeline of jobs that are executed after one another.
    dependency_type : str, default='afterok
        How to handle the 'dependency_jobs'. Must be one of ['after',
        'afterany', 'afternotok', 'afterok', 'singleton']. See [1] for further
        information.
    job_dir : str
        Name of the batch job for creating meaningful log file names.

    Returns
    -------
    job_id : int
        The job ID of the submitted SLURM batch job.

    Notes
    -----
    [1] https://hpc.nih.gov/docs/job_dependencies.html
    """

    # Join arguments to a single bash command
    cmd = " ".join(str(arg) for arg in args_list)

    # Create directory for output logs
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    error = f"{log_dir}/slurm-%j-{job_name}.out"
    output = f"{log_dir}/slurm-%j-{job_name}.out"

    # Prepare job scheduler
    slurm = Slurm(
        cpus_per_task=cpus,
        error=error,
        mem=mem,
        nodes=1,
        ntasks=1,
        output=output,
        time=time,
        job_name=job_name,
        tmp="300G",
        mail_type="all",
        mail_user="yoos@cbs.mpg.de",
    )
    slurm.add_cmd('check_ComputeClusterSlurm_memory-usage')

    # Make the current job depend on previous jobs
    if dependency_jobs != []:
        if isinstance(dependency_jobs, int):
            dependency_jobs = [dependency_jobs]
        dependency_str = ":".join([str(job_id) for job_id in dependency_jobs])
        dependency = {dependency_type: dependency_str}
        slurm.set_dependency(dependency)

    # Submit
    print("Submitting", cmd)
    job_id = slurm.sbatch(cmd)

    return job_id

# Define directory for SLURM batch job log files
log_dir = Path("/data/p_02915/dhcp_derivatives_SPOT/")

script = "python3 /data/p_02915/SPOT/run_all.py"


# Iterate over the index array using a for loop
for index in range(34,41): 
    subj_num = index
    args = [script, subj_num]
    job_id = submit_job(
             args,
             mem=320000,
             cpus=2,
             log_dir=log_dir,
             job_name="dhcp_retinotopy"
             )


#job_name="dhcp_preprocessing" 9-21
