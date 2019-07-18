% steady state test for w'c' according to Foken and Wichura (1996)
% see Chapter 9 in X. Lee et al. (eds.), Handbook of Micrometeorology (2004)
function sst = steady_state_test(w_prime, c_prime)

    subset_size = floor(length(w_prime)/5);
    wc_subset = zeros(1,5);
    for (j = 1:5)
        wc_subset(j) = xcov(w_prime(((j-1)*subset_size+1):(j*subset_size)), c_prime(((j-1)*subset_size+1):(j*subset_size)), 0)/subset_size;
    end
    wc_total = xcov(w_prime, c_prime, 0)/length(w_prime);
    sst = abs((mean(wc_subset) - wc_total)/wc_total);

end
