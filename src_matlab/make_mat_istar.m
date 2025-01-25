function Istar = make_mat_istar(K);

n_el=prod(1:K)/prod(1:K-2);

Imat=perms([1 2 zeros(1,K-2)]);
j_unrest=size(Imat,1);
j=1;
cnt=2;
Istar=Imat(1,:);

while j<n_el;
    if any(all(Istar==permute(Imat(cnt,:,ones(size(Istar,1),1)),[3,2,1]),2))
        cnt=cnt+1;
    else
        Istar=[Istar; Imat(cnt,:)];
        cnt=cnt+1;
        j=j+1;
    end
end
 