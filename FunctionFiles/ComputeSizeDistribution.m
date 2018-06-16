function size_distribution = ComputeSizeDistribution(cell_loss_mat,LOSS_RATIO)
%*****************************************************************************
%This function uses "cell_loss_mat", result of simulations indicating
%positions of cell loss, "LOSS_RATIO", measured ratio of cell loss. The
%output is "size_distribution", frequency distribution of cluster size.
%*****************************************************************************

%% Compute frequency distribution of cluster size from simulated data
N = size(cell_loss_mat);
LostCell = floor(N(1)*N(2)*LOSS_RATIO);
Record = zeros(N(3),LostCell);

for i = 1:N(3)
clusterNo = bwlabeln(cell_loss_mat(:,:,i),8);
    for k = 1:max(max(clusterNo))
    Kmat = clusterNo == k;
    cellNo = sum(sum(Kmat));
    Record(i, k) = cellNo;
    end
end

RecordVec = reshape(Record, N(3)*LostCell, 1);
Index = RecordVec > 0;
RecordVec2 = RecordVec(Index);
size_distribution = tabulate(RecordVec2);
