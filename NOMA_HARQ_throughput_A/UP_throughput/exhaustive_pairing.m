function [exhaustive_pairing] = exhaustive_pairing(user_distance)
    mother_set = 1:10;
    exhaustive_paring = zeros(5,2,945);
    ind=0;

    for j=1:9
        for p=1:7
            for q=1:5
                for k=1:3
                    ind= ind+1;
                    exhaustive_paring(1,1,ind) = 1;
                    exhaustive_paring(1,2,ind) = j+1;

                    set1=mother_set;
                    set1(set1 == 1)=[];
                    set1(set1 == (j+1))=[];
                    min1 = min(set1);
                    set1(set1 == min1) = [];
                    exhaustive_paring(2,1,ind) = min1;
                    exhaustive_paring(2,2,ind) = set1(p);


                    set1(set1 == set1(p)) = [];
                    min2 = min(set1);
                    set1(set1 == min2) = [];
                    exhaustive_paring(3,1,ind) = min2;
                    exhaustive_paring(3,2,ind) = set1(q);


                    set1(set1 == set1(q)) = [];
                    min3 = min(set1);
                    set1(set1 == min3) = [];
                    exhaustive_paring(4,1,ind) = min3;
                    exhaustive_paring(4,2,ind) = set1(k);

                    set1(set1 == set1(k)) = [];
                    exhaustive_paring(5,1,ind) = set1(1);
                    exhaustive_paring(5,2,ind) = set1(2);

                end
            end
        end
    end
    exhaustive_pairing = user_distance(exhaustive_paring);
end


