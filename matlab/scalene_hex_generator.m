a = 1;
L = 1;
d = 1;


number_of_plates_in_chain = 1;
number_of_chains_to_generate = 20;
d_array = zeros(1,number_of_chains_to_generate);

basic_hex = [-a/2 -a*d*sqrt(3)/2 L/2; a/2 -a*d*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*d*sqrt(3)/2 -L/2; a/2 -a*d*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];

collection_of_chains = zeros(12, 3, number_of_plates_in_chain, number_of_chains_to_generate);

for current_chain = 1:number_of_chains_to_generate
    
    basic_hex = [-a/2 -a*d*sqrt(3)/2 L/2; a/2 -a*d*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*d*sqrt(3)/2 -L/2; a/2 -a*d*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];
    
    collection_of_chains(:,:,1, current_chain) = basic_hex;
    d_array(current_chain) = d;
    d = d - (1/number_of_chains_to_generate);
    
end