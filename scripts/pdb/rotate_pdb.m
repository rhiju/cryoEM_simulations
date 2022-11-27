function pdbstruct = rotate_pdb( pdbstruct, M )
% pdbstruct = rotate_pdb( pdbstruct, M )
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
    rotate_r(:,i) = M * all_r(:,i);
    pdbstruct.Model(1).Atom(i).X = rotate_r(1,i);
    pdbstruct.Model(1).Atom(i).Y = rotate_r(2,i);
    pdbstruct.Model(1).Atom(i).Z = rotate_r(3,i);
end
