for plates=5:5
    
   
    number_of_plates_in_chain = plates

    run('superloop_generator.m');

    message = 'starting capacitance calculation'
    
    collection_of_chains_capacitance_calculator

    numb = num2str(plates);
    filename = [numb,'_hex_chain_ultraloop_run_ar_0_25'];
    save(filename);
    clear;

end