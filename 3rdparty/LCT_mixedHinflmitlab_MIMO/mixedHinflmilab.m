function [ K, sol_info ] = mixedHinflmilab( sys, ny, nu, alpha, channel)

gen_sys = extract_generalized_plant(sys,nu,ny);

n_pspecs = length(channel); %number of performance specs
A = gen_sys.A; nx = gen_sys.nx; Bu = gen_sys.Bu; Cy = gen_sys.Cy; 
[Bw,Cz,Dzw,Dzu,Dyw] = deal(cell(1,n_pspecs));
[nw,nz] = deal(zeros(1,n_pspecs));

for j = 1:n_pspecs
    Bw{j}  = gen_sys.Bw(:,channel(j).In);
    Cz{j}  = gen_sys.Cz(channel(j).Out,:);   
    Dzw{j} = gen_sys.Dzw(channel(j).Out,channel(j).In); 
    Dzu{j} = gen_sys.Dzu(channel(j).Out,:); 
    Dyw{j} = gen_sys.Dyw(:,channel(j).In); 
    nw(j)  = size(Bw{j},2);
    nz(j)  = size(Cz{j},1);
end

setlmis([]);

X = lmivar(1,[nx,1]);
Y = lmivar(1,[nx,1]);

Ac_hat = lmivar(2,[nx nx]);
Bc_hat = lmivar(2,[nx ny]);
Cc_hat = lmivar(2,[nu nx]);

for j = 1:n_pspecs
    gam2{j} = lmivar(2,[1 1]);
end

Q = newlmi;
lmiterm([Q,1,1,X],1,1);
lmiterm([Q,1,2,0],eye(nx));
lmiterm([Q,2,2,Y],1,1);

for j = 1:n_pspecs
Term{j} = newlmi;
lmiterm([Term{j},1,1,X],A,1,'s');
lmiterm([Term{j},1,1,Cc_hat],Bu,1,'s');
lmiterm([Term{j},1,2,Ac_hat],1,1);
lmiterm([Term{j},1,2,0],A);
lmiterm([Term{j},2,2,Y],1,A,'s');
lmiterm([Term{j},2,2,Bc_hat],1,Cy,'s');
lmiterm([Term{j},1,3,0],Bw{j});
lmiterm([Term{j},3,3,0],-eye(nw(j)));
lmiterm([Term{j},1,4,X],1,Cz{j}');
%lmiterm([Term{j},1,4,Cc_hat],Dzu{j}',1);
lmiterm([Term{j},2,4,0],Cz{j}');
lmiterm([Term{j},3,4,0],Dzw{j});
lmiterm([Term{j},4,4,gam2{j}],-1,eye(nz(j)));

end
LMIs = getlmis;
c = zeros(decnbr(LMIs),1);

for j = 1:n_pspecs
   G{j} = decinfo(LMIs,gam2{j}); 
   f = G{j};
   c(f(f~=0)) = alpha((f~=0));
end
[fopt,xopt] = mincx(LMIs,c);
LMIsopt = evallmi(LMIs,xopt);

K = 0;


end

