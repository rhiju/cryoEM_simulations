function pdbstruct = center_pdb( pdbstruct )
% pdbstruct = center_pdb( pdbstruct )
%
%

% Let's apply an offset
for i = 1:length(pdbstruct.Model.Atom)
    all_r(:,i) = [pdbstruct.Model.Atom(i).X,pdbstruct.Model.Atom(i).Y,pdbstruct.Model.Atom(i).Z];
end

mean_r = mean(all_r');
mean_x = mean_r(1);
mean_y = mean_r(2);
mean_z = mean_r(3);

for i = 1:length(pdbstruct.Model.Atom)
    pdbstruct.Model.Atom(i).X = pdbstruct.Model.Atom(i).X - mean_x;
    pdbstruct.Model.Atom(i).Y = pdbstruct.Model.Atom(i).Y - mean_y;
    pdbstruct.Model.Atom(i).Z = pdbstruct.Model.Atom(i).Z - mean_z;
end

for i = 1:length(pdbstruct.Model.Atom)
    all_r(:,i) = [pdbstruct.Model.Atom(i).X,pdbstruct.Model.Atom(i).Y,pdbstruct.Model.Atom(i).Z];
end
