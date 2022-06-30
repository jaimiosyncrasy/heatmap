import sys

import run_RHP
import run_evalConfig
import run_createPowerGrid
import run_CPP

print('im in driver')

# setup power grid
grid_vars_obj=run_createPowerGrid.createPowerGrid_func() # grid_vars_obj holds loose grid specs

# handle run-config parms and run the right experiment
print("Argument List:", str(sys.argv))  # first arg is always the name of the file being run
exp_num=sys.argv[1]  # extract first parm from the run-config
if exp_num=='1.a':
    run_evalConfig.exp1a(grid_vars_obj)
elif exp_num=='1.b':
    run_evalConfig.exp1b(grid_vars_obj)
elif exp_num=='2':
    run_evalConfig.exp2(grid_vars_obj)
elif exp_num == '4':
    run_CPP.exp4a(grid_vars_obj)
else:
    raise Exception("unrecognized run-config parm for experiment number")
print('------ run-config complete!-------------')