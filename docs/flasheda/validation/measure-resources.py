import json
from pathlib import Path
import resource
import subprocess
import sys

output = Path(sys.argv[1])
result = subprocess.run(sys.argv[2:])
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
output.write_text(json.dumps(dict(max_rss_kb=usage.ru_maxrss,
                                  user_seconds=usage.ru_utime,
                                  system_seconds=usage.ru_stime),indent=2)+'\n')
sys.exit(result.returncode)
