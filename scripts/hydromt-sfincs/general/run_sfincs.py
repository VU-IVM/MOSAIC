#!/bin/sh
import shutil
import subprocess

source = r'run.sh'
destination = r'xynthia_sfincs_pyscript2/run.sh'

#shutil.copyfile(source, destination)
#chmod u+rwx xynthia_sfincs_pyscript2/run.bat
subprocess.call([r'xynthia_sfincs_pyscript2/run.sh'])

#os.system(r'xynthia_sfincs_pyscript2\run.bat')
