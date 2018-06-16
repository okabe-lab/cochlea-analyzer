function cellLossMat = SimulateCombinationEffect(COCHLEAR_SIZE, LOSS_RATIO ...
    , neighborEffect, positionEffect, TRIAL_NO)
%*****************************************************************************
%This function uses "COCHLEAR_SIZE", numbers of rows and columns of virtual
%outer hair cells, "LOSS_RATIO", measured ratio of cell loss, "neighborEffect",
%effect size of neighborhood effect, "positionEffect", effect size of
%positional effect, and "TRIAL_NO". The output is "cellLossMat", result of
%simulations indicating positions of cell loss.
%*****************************************************************************

powerInd = 0.01;
gaussInd = 3;
positionEffect = 1-1/(exp(positionEffect));

cellLossMat = zeros(COCHLEAR_SIZE(1), COCHLEAR_SIZE(2), TRIAL_NO);
LostCell = floor(COCHLEAR_SIZE(1)*COCHLEAR_SIZE(2)*LOSS_RATIO);

for i=1:TRIAL_NO
    Power= 0.1*powerInd*rand(COCHLEAR_SIZE(2),1).^-(1-0.1*powerInd);
    Power2 = Power*positionEffect+ones(size(Power))*(1-positionEffect);
    GPower =imgaussfilt(Power2,gaussInd);
    GPower = (GPower/sum(GPower)/3)';
    
    probMat = cat(1,GPower,GPower,GPower);
    
    SumGPower = cumsum(probMat(:));
    judge = 0;
    Ncount = 0;
    while judge == 0
        positionXY = find(SumGPower > rand, 1);
        [positionX, positionY] = ind2sub(size(probMat),positionXY);
        
        if cellLossMat(positionX, positionY, i) ==0
            cellLossMat(positionX, positionY, i) =1;
            probMat = padarray(probMat,[1,1]);
            probMat(positionX:positionX+2, positionY:positionY+2) =...
                probMat(positionX:positionX+2, positionY:positionY+2) ...
                + (neighborEffect*(1/numel(SumGPower)));
            probMat = probMat(2:end-1,2:end-1);
            probMat(positionX,positionY) = 0;
            probMat = probMat/sum(probMat(:));
            SumGPower = cumsum(probMat(:));
            Ncount = Ncount +1;
        end
        if Ncount == LostCell
            judge= 1;
        end
    end
end