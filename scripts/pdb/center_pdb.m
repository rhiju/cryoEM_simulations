function pdbstruct = center_pdb( pdbstruct )
% pdbstruct = center_pdb( pdbstruct )
%
%

% Let's apply an offset
for i = 1:length(pdbstruct.Model(1).Atom)
    all_r(:,i) = [pdbstruct.Model(1).Atom(i).X,pdbstruct.Model(1).Atom(i).Y,pdbstruct.Model(1).Atom(i).Z];
end

mean_r = mean(all_r');
mean_x = mean_r(1);
mean_y = mean_r(2);
mean_z = mean_r(3);

for i = 1:length(pdbstruct.Model(1).Atom)
    pdbstruct.Model(1).Atom(i).X = pdbstruct.Model(1).Atom(i).X - mean_x;
    pdbstruct.Model(1).Atom(i).Y = pdbstruct.Model(1).Atom(i).Y - mean_y;
    pdbstruct.Model(1).Atom(i).Z = pdbstruct.Model(1).Atom(i).Z - mean_z;
end

for i = 1:length(pdbstruct.Model(1).Atom)
    all_r(:,i) = [pdbstruct.Model(1).Atom(i).X,pdbstruct.Model(1).Atom(i).Y,pdbstruct.Model(1).Atom(i).Z];
end
