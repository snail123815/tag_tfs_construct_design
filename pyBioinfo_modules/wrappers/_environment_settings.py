
# TODO: write a script change this file
from typing import Literal
from pathlib import Path


SHELL: Literal['bash', 'zsh'] = 'zsh'
CONDAEXE: Literal['conda', 'mamba', 'micromamba'] = 'micromamba'

ANTISMASH_ENV: Path | None = Path.home() / 'genvs/quasan'
BUSCO_ENV: Path | None = Path.home() / 'genvs/quasan'
PROKKA_ENV: Path | None = Path.home() / 'genvs/quasan'
MASH_ENV: Path | None = Path.home() / 'genvs/phylophlan'
BIGSCAPE_ENV: Path | None = Path.home() / 'genvs/bigscape'
HMMER_ENV: Path | None = Path.home() / '.micromamba'
PFAM_DB: Path | None = None
try:
    PFAM_DB = sorted([p for p
                      in (Path.home() / 'dbMisc/antismash_databases/pfam'
                          ).iterdir() if p.is_dir()])[-1]
except FileNotFoundError:
    pass


def withActivateEnvCmd(cmd: str, condaEnv: Path | None = None,
                       condaExe=CONDAEXE, shell=SHELL) -> str:
    '''cmd = withActivateEnvCmd(cmd, condaEnv, condaExe, shell)
    Join commands that activate your environment and your command,
    return joined command:
    eval "$(micromamba shell hook --shell={shell} &&
    micromamba activate {condaEnv} && {cmd}
    '''
    if condaEnv is not None:
        if condaExe == 'micromamba':
            activateEnvCmd = f'eval "$(micromamba shell hook --shell={shell})"'
        else:
            activateEnvCmd = f'eval "$({condaExe} shell.{shell} hook)"'
        activateEnvCmd += f' && micromamba activate {condaEnv}'
        cmd = activateEnvCmd + " && " + cmd
    return cmd
