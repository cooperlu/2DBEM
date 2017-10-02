function [ H,G ] = applybc( mu,ne,bct,H,G )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:ne
    for j = 1:6 % 6 unknowns each element
        if bct(6*i-6+j) == 0 % 0 is displacement, 1 is traction
            if (i == ne) && j > 4 % 1st code in the last element
                if bct(j-4) > 0
                    ch = H(:,j-4);
                    H(:,j-4) =  - G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = - ch;
                else
                    H(:,j-4) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));

                end
            else
                if (i == 1) || (j>2) || (bct(6*i-8+j) ==1)
                    ch = H(:,4*i-4+j);
                    H(:,4*i-4+j) =  - G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = - ch;
                else
                    H(:,4*i-4+j) = H(:,4*i-4+j)- G(:,6*i-6+j)*mu;
                    G(:,6*i-6+j) = zeros(size(G(:,6*i-6+j)));
                end
            end
        end
    end
end
end

