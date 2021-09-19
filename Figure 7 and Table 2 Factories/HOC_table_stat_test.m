load('HOC_table_Fib_1000trials')

% Locations in table for respective orders
real_vals_fib = [379165 11261 7782 2986 1191 679;
                 181554 3254 2519 1214 276 153;
                 98272 1021 831 473 63 35;
                 142575 544 477 341 24 13]

p_vals = zeros(4,6);
for i = 1:4
    for j = 3:6
        temp = zeros(1000,1);
        for p = 1:1000
            temp(p) = HOC_table_Fib{p}(i,j);
        end
        p_vals(i,j) = sum(temp>real_vals_fib(i,j))/1000;
    end
end



figure
hist(temp)
%%
real_vals_gm = [240477 8384 7384 4157 2006 1536;
                227352 4354 3972 2686 822 606;
                196423 1996 1881 1434 277 193;
                1000231 1802 1727 1419 109 67]

p_vals = zeros(4,6);
for i = 1:4
    for j = 3:6
        temp = zeros(1000,1);
        for p = 1:1000
            temp(p) = HOC_table_GM{p}(i,j);
        end
        p_vals(i,j) = sum(temp>real_vals_gm(i,j))/1000;
    end
end