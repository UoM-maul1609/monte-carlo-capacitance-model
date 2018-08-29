% make sure comment number_of_plates_in_chain in collection_of_chain_generator_with_intersection_beta
for number_of_plates_in_chain=2:20
    collection_of_chain_generator_with_intersection_beta;
    run2(number_of_plates_in_chain-1).collection_of_chains=collection_of_chains;
    run2(number_of_plates_in_chain-1).D_maxes=D_maxes;                          
    run2(number_of_plates_in_chain-1).aggregation_indices=aggregation_indices;  
    run2(number_of_plates_in_chain-1).number_of_plates_in_chain=number_of_plates_in_chain;  
    number_of_plates_in_chain
end

for number_of_plates_in_chain=2:20
    collection_of_chains=run2(number_of_plates_in_chain-1).collection_of_chains;   
    for k=1:10
        rotate_about_cog_pjc;
    end
    run2(number_of_plates_in_chain-1).D_m1=D_m1;
    run2(number_of_plates_in_chain-1).vt=vt;    
    run2(number_of_plates_in_chain-1).v1=v1;
    run2(number_of_plates_in_chain-1).Ar=Ar;
end
for number_of_plates_in_chain=2:20
    collection_of_chains=run1(number_of_plates_in_chain-1).collection_of_chains;   
    for k=1:10
        rotate_about_cog_pjc;
    end
    run1(number_of_plates_in_chain-1).D_m1=D_m1;
    run1(number_of_plates_in_chain-1).vt=vt;    
    run1(number_of_plates_in_chain-1).v1=v1;
    run1(number_of_plates_in_chain-1).Ar=Ar;
end
