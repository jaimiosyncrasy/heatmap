import sys

import run_RHP
import run_evalConfig
import run_createPowerGrid
import run_CPP
import run_OCPP
import matplotlib.pyplot as plt

print('im in driver')

# setup power grid
grid_vars_obj=run_createPowerGrid.createPowerGrid_func() # grid_vars_obj holds loose grid specs

# handle run-config parms and run the right experiment
print("Argument List:", str(sys.argv))  # first arg is always the name of the file being run
exp_num=sys.argv[1]  # extract first parm from the run-config

if exp_num=='1eval':
    run_evalConfig.exp1(grid_vars_obj)
elif exp_num=='1view':
    run_evalConfig.exp1_markGraph(grid_vars_obj)
if exp_num == '2':
    run_evalConfig.exp2(grid_vars_obj)
elif exp_num == '2view':
    run_evalConfig.exp2_markGraph(grid_vars_obj)
elif exp_num == '3eval': #new
    run_evalConfig.exp3eval(grid_vars_obj)
elif exp_num == '3':
    run_RHP.exp3(grid_vars_obj)
elif exp_num == '3view':
    run_RHP.exp3_viewResults(grid_vars_obj)
elif exp_num == '4':  # new
    run_evalConfig.exp4(grid_vars_obj)
elif exp_num == '4.a': # bad branch
    run_CPP.exp4a(grid_vars_obj)
elif exp_num == '4.b': # good branch
    run_CPP.exp4b(grid_vars_obj)
elif exp_num == '4view':
    run_CPP.make_CPP_heatmaps(grid_vars_obj)
elif exp_num == '5':
    run_OCPP.exp5(grid_vars_obj)
elif exp_num == '6': #new
    run_evalConfig.exp6(grid_vars_obj)
elif exp_num=='1234':
    run_evalConfig.exp1(grid_vars_obj)
    run_evalConfig.exp2(grid_vars_obj)
    run_evalConfig.exp3eval(grid_vars_obj)
    run_evalConfig.exp4(grid_vars_obj)
else:
    raise Exception("unrecognized run-config parm for experiment number")
print('------ Complete! close any plotting windows to finish the run-config -------------')

plt.show(block=True)