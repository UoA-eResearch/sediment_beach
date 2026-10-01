function nd = compare_outputs(file_a, file_b)
% COMPARE_OUTPUTS  Compare the O structure of two ShorelineS output.mat files
%   Prints the maximum absolute difference of every numeric field that is not
%   identical, and returns the number of differing fields.
A = load(file_a);
B = load(file_b);
f = fieldnames(A.O);
nd = 0;
for k = 1:numel(f)
    x = A.O.(f{k});
    y = B.O.(f{k});
    if isnumeric(x) && ~isequaln(x,y)
        nd = nd + 1;
        if ~isequal(size(x),size(y))
            fprintf('O.%s: size differs\n', f{k});
            continue
        end
        d = abs(double(x(:))-double(y(:)));
        fprintf('O.%s: max abs diff = %g (same NaN pattern = %d)\n', f{k}, max(d(~isnan(d))), isequal(isnan(x),isnan(y)));
    end
end
fprintf('%d of %d output fields differ\n', nd, numel(f));
end
