def exp2a():
    # ------- run runHeatMapProcess (RHP) ---------------
    # 6.18.22
    # run RHP on neighborhood config for 1st to 7th step:
    parmObj = hm.configParms()
    parmObj.set_version(1)  # 1 for PBC

    # set_act_locs = []
    # set_perf = []
    # addon_act_nodes = ['bus_82','bus_87','bus_76','bus_49','bus_41','bus_46','bus_60']
    # addon_perf_nodes = ['bus_77','bus_77','bus_77','bus_44','bus_44','bus_44','bus_66']

    # run RHP on neighborhood config for 7th step only:
    set_act_locs = ['bus_82', 'bus_87', 'bus_76', 'bus_49', 'bus_41', 'bus_46']
    set_perf = ['bus_77', 'bus_77', 'bus_77', 'bus_44', 'bus_44', 'bus_44']
    addon_act_nodes = ['bus_60']
    addon_perf_nodes = ['bus_66']

    parmObj.set_ctrlTypes(['PBC'] * len(set_act_locs + addon_act_locs))

    t = time.time()
    lst_feas_configs, lzn_error_run_sum, heatMapNames = prc.runHeatMapProcess(parmObj, fin_feeder, set_act_locs,
                                                                              set_perf, addon_act_nodes,
                                                                              addon_perf_nodes, node_index_map,
                                                                              substation_name, depths, file_name,
                                                                              Vbase_ll, Sbase)
    # lst_feas_configs is list of dictionaries
    elapsed = time.time() - t
    print('Time elapsed=', elapsed)
