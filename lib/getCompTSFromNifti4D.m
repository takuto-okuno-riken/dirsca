%%
% get Component time-series (matrix) from NIfTI 4D volume.
% returns Component time-series (X)
% input:
%  V            nifti 4D volume (X x Y x Z x frames)
%  compV        nifti 3D atlas (X x Y x Z x component)

function X = getCompTSFromNifti4D(V, compV)
    
    X = single(zeros(size(compV,4), size(V,4)));
    A = reshape(V,[],size(V,4));
    for i=1:size(compV,4)
        C = compV(:,:,:,i);
        m = C(:)' * A;
        X(i,:) = m;
    end
end

