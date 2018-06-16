function cellLossMat = SimulateNeighborEffect(COCHLEAR_SIZE, LOSS_RATIO ...
    , neighbor_effect, TRIAL_NO)
%*****************************************************************************
%This function uses "COCHLEAR_SIZE", numbers of rows and columns of virtual
%outer hair cells, "LOSS_RATIO", measured ratio of cell loss, "neighbor_effect",
%effect size, and "TRIAL_NO". The output is "cellLossMat", result of
%simulations indicating positions of cell loss.
%*****************************************************************************

% Simulate cell loss on specified size of "Neighborhood Effect".
cellLossMat = zeros(COCHLEAR_SIZE(1), COCHLEAR_SIZE(2), TRIAL_NO);
LostCell = floor(COCHLEAR_SIZE(1)*COCHLEAR_SIZE(2)*LOSS_RATIO);

for k = 1:TRIAL_NO
    cellMat = zeros(COCHLEAR_SIZE(1), COCHLEAR_SIZE(2));
    probMat = zeros(COCHLEAR_SIZE(1)+2, COCHLEAR_SIZE(2)+2)+1;
    
    LostNo = 0;
    judge = 0;
    while judge ==0
        judge2 = 0;
        while judge2 ==0
            linearprobMat = reshape(probMat(2:COCHLEAR_SIZE(1)+1, 2:COCHLEAR_SIZE(2)+1)'...
                , 1, COCHLEAR_SIZE(1)*COCHLEAR_SIZE(2));
            sumlinearprobMat = cumsum(linearprobMat);
            randomposition = ceil(max(sumlinearprobMat)*rand);
            linearposition = find(sumlinearprobMat > randomposition, 1);
            positionY = ceil(linearposition/COCHLEAR_SIZE(2));
            positionX = rem(linearposition,COCHLEAR_SIZE(2));
            if positionX ==0
                positionX = COCHLEAR_SIZE(2);
            end
            
            if cellMat(positionY, positionX) == 0
                cellMat(positionY, positionX) = 1;
                judge2 = 1;
                probMat(positionY:positionY+2, positionX:positionX+2) ...
                    = probMat(positionY:positionY+2, positionX:positionX+2) +neighbor_effect;
                LostNo = LostNo +1;
            end
            
        end
        if LostNo == LostCell
            judge =1;
        end
    end
    cellLossMat(:,:,k) = cellMat;
end

