# Activate your Python environment first. All arguments are forwarded unchanged.
& python (Join-Path $PSScriptRoot "sem_cd_measure_200k_batch_V1_16.py") @args
exit $LASTEXITCODE
