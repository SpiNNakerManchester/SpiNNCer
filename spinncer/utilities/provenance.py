"""
Utility script to retrieve the HASH of the current git commit

Running command in bash:
https://stackoverflow.com/questions/4256107/running-bash-commands-in-python

Retrieving current commit HASH:
https://stackoverflow.com/questions/949314/how-to-retrieve-the-hash-for-the-current-commit-in-git
"""


def retrieve_git_commit():
    import subprocess
    from subprocess import PIPE
    bash_command = "git rev-parse HEAD"

    try:
        proc = subprocess.run(bash_command.split(),
                              stdout=PIPE, stderr=PIPE, shell=False)
        return proc.stdout
    except subprocess.CalledProcessError as e:
        print("Failed to retrieve git commit HASH-", str(e))
        return "CalledProcessError"
    except Exception as e:
        print("Failed to retrieve git commit HASH more seriously-", str(e))
        return "GeneralError"
