import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

JOBSCRIPT_TEMPLATE = """#!/bin/bash

#SBATCH {sbatch_options}

echo ________________________________________
echo
echo SLURM Job Log
echo Start time: "$(date)"
echo
echo Job name: "$SLURM_JOB_NAME"
echo Job ID: "$SLURM_JOBID"
echo Submitted by user: "$USER"
echo User effective group ID: "$(id -ng)"
echo
echo SLURM account used: "$SLURM_ACCOUNT"
echo Hostname of submission: "$SLURM_SUBMIT_HOST"
echo Submitted to cluster: "$SLURM_CLUSTER_NAME"
echo Submitted to node: "$SLURMD_NODENAME"
echo Cores on node: "$SLURM_CPUS_ON_NODE"
echo Requested cores per task: "$SLURM_CPUS_PER_TASK"
echo Requested cores per job: "$SLURM_NTASKS"
echo Requested walltime: "$SBATCH_TIMELIMIT"
echo Nodes assigned to job: "$SLURM_JOB_NODELIST"
echo Running node index: "$SLURM_NODEID"
echo
echo Running on hostname: "$HOSTNAME"
echo Parent PID: "$PPID"
echo Process PID: "$$"
echo
echo Working directory: "$SLURM_SUBMIT_DIR"
echo ________________________________________________________
echo

cd "$SLURM_SUBMIT_DIR" || exit

source "$HOME/.bashrc"; conda activate {conda_env_name}

COMMAND='{command_to_run}'

echo "$COMMAND"
time eval "$COMMAND"
"""


@dataclass
class ConfigurableJobscript:
    jobscript_name_prefix: str
    jobscript_dir: Path
    sbatch_options: list[str]
    conda_env_name: str
    command_to_run: str
    _init_at: datetime = datetime.now()

    @property
    def path(self) -> Path:
        """Path to the generated jobscript"""
        timestamp = self._init_at.strftime("%Y%m%dT%H%M%S")
        filename = f"{self.jobscript_name_prefix}_{timestamp}.sh"
        return self.jobscript_dir / filename

    def content(self) -> str:
        """Return the content of the generated jobscript as a multiline string"""
        return JOBSCRIPT_TEMPLATE.format(
            sbatch_options="\n#SBATCH ".join(self.sbatch_options),
            conda_env_name=self.conda_env_name,
            command_to_run=self.command_to_run,
        )

    def write(self) -> Path:
        """Write the jobscript to the output directory and return the path to the file"""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("w") as f:
            f.write(self.content())

    def submit(self) -> None:
        """Submit the jobscript to the cluster with sbatch"""
        if not self.path.exists():
            self.write()
        subprocess.run(["sbatch", f"{self.path}"], check=True)
