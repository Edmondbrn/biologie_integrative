import os
import subprocess
os.chdir(os.path.dirname(__file__))
print(os.listdir())
file= "test.R"
Rscript = "C:\\Program Files\\R\\R-4.4.0\\bin\\Rscript.exe"
if os.path.isfile(file):
    res = subprocess.run([Rscript,file], capture_output=True, text=True)
    print(res.stdout)