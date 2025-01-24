function IMSf = forecast_MS(IMSf,etaMS);

% IMS=IMSf+1;
nst=size(etaMS,1);

if nst>2
IMSenc=[[0;0] [0;1] [1;0] [1;1]];

f= etaMS'*all(kron(ones(1,4),IMSf)==IMSenc,1)';

c=rand(sum(f>0),1);
Iind=find(f);
IMSf = IMSenc(:,Iind(sum(cumsum(f(f>0)) < c) + 1));          

elseif nst==2
    f=etaMS'*([1:nst]'==IMSf+1);
    c=rand(nst,1);
    IMSf = sum(cumsum(f)<c);
end


