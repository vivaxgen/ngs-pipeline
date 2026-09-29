shell_min_delay = max(0, int(config.get("shell_min_delay", 0)))

if shell_min_delay > 0:
    shell.prefix(f"set -euo pipefail; sleep {shell_min_delay}; ")
