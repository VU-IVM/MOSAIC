import subprocess
import sys
import configparser

# setting up the ini file
config = configparser.ConfigParser()
config.read('/gpfs/home4/benitoli/papers/paper1/scripts/hydromt-sfincs/general/sfincs_coastal.ini')
config.set('setup_config', 'tref', sys.argv[2])
config.set('setup_config', 'tstart', sys.argv[2])
config.set('setup_config', 'tstop', sys.argv[3])
config.set('setup_basemaps', 'basemaps_fn', sys.argv[5])
config.set('setup_h_forcing', 'geodataset_fn', sys.argv[6])
config.set('setup_gauges', 'gauges_fn', '/gpfs/home4/benitoli/papers/paper1/scripts/hydromt-sfincs/general/' + sys.argv[7])
with open('sfincs_coastal.ini', 'w') as configfile:
    config.write(configfile)

# executing hydromt
cmd = ["hydromt", "build", "sfincs", sys.argv[1], sys.argv[4], "-r", "150", "-i", "sfincs_coastal.ini", "-vv"]
subprocess.call(cmd)