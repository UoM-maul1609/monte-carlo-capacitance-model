a = 1;
L = 38;
number_of_plates_in_chain = 1;

number_of_chains_to_generate = 15;

collection_of_chains = zeros(12, 3, number_of_plates_in_chain, number_of_chains_to_generate);
aspect_rations = zeros(1, number_of_chains_to_generate);



for i=1:number_of_chains_to_generate
    L= L +i*0.1;
    
    basic_hex = [-a/2 -a*sqrt(3)/2 L/2; a/2 -a*sqrt(3)/2 L/2; a 0 L/2; a/2 a*sqrt(3)/2 L/2; -a/2 a*sqrt(3)/2 L/2;  -a 0 L/2; -a/2 -a*sqrt(3)/2 -L/2; a/2 -a*sqrt(3)/2 -L/2; a 0 -L/2; a/2 a*sqrt(3)/2 -L/2; -a/2 a*sqrt(3)/2 -L/2;  -a 0 -L/2];
    
    collection_of_chains(:,:,1, i) = basic_hex;
    
    aspect_rations(i) = L/(2*a);
    

end
